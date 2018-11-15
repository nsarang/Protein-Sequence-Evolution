#include "alignment.h"


// MULTI-THREADING
std::mutex al_mtx_1, al_mtx_2, al_mtx_print;


// ALIGNMENT REPOSITORY
std::vector<std::vector<std::vector<std::string>>> matAMR;


// SCORE PROFILES
int TEMPLATE_COUNT;
std::array<std::array<double, 7>, MAX_PROT_LEN> alignment_potential;


// PREPROCESS ARRAYS
std::array<std::array<double, 7>, MAX_PROT_LEN> alignment_position_score;
std::array<long long, MAX_PROT_LEN> alignment_position_sum;
std::array<long long, 7> amino_acid_template_frequency;




void StructureAlignment(std::string target, std::string fPath, bool prcss_frgmnts) {
    std::string stdout, line, sequence_1, sequence_2;
    while (system_call_err(ex_TMALIGN + " " + target + "  " + fPath + " -a 2>&1", stdout) != 0);

    std::stringstream ret(stdout);


    double score;
    getline(ret, line);
    getline(ret, line);
    ret >> line >> score;
    ret >> sequence_1 >> sequence_2;

    assert(sequence_1.size() == sequence_2.size());
    int align_len = sequence_1.size();

    std::vector<double> dist_sc(align_len + 1),
        accum_dist(align_len + 1);
    accum_dist[0] = 0;
    for (int i = 1; i <= align_len; ++i)
    {
        ret >> dist_sc[i];
        accum_dist[i] = (dist_sc[i] == GAP_SCORE ? GAP_PENALTY : dist_sc[i])
                        + accum_dist[i - 1];
    }


    if (score >= ALIGN_SC_CUTTOFF) {
        std::lock_guard<std::mutex> lock(al_mtx_1);
        TEMPLATE_COUNT++;
        int pos = 0;
        for (int i = 0; i < align_len; ++i)
        {
            if (sequence_1[i] == '-')
                continue;
            pos++;

            if (sequence_2[i] == '-' || sequence_2[i] == 'X')
                continue;
            assert(AA_classification.count(sequence_2[i]) > 0);
            int class_num = AA_classification[sequence_2[i]];

            alignment_position_score[pos - 1][class_num] += 10 - std::min(10., dist_sc[i + 1]);
            alignment_position_sum[pos - 1]++;
            amino_acid_template_frequency[class_num]++;
        }
    }

    if (prcss_frgmnts) {
        int current_p = 0;
        for (int i = 0; i < align_len; ++i) {
            for (int l = FRAG_MIN_LEN; i + l <= align_len; ++l) {
                if (dist_sc[i + 1] == 0 || dist_sc[i + l] == 0)
                    continue;
                double sc = accum_dist[i + l] - accum_dist[i];
                assert(sc > 0);
                if (sc <= l * DIST_CUTOFF)
                {
                    std::string fragment;
                    int real_l = 0;
                    for (int j = i; j < i + l; ++j) {
                        if (sequence_1[j] == '-')
                            continue;
                        fragment += sequence_2[j];
                        real_l++;
                    }
                    std::lock_guard<std::mutex> lock(al_mtx_2);
                    matAMR[current_p][current_p + real_l - 1].push_back(fragment);
                }
            }
            if (sequence_1[i] != '-')
                current_p++;
        }
    }
}


void Alignment_Eval(std::string target) {
    std::ofstream outFile(db_Align + File_md5(target));

    long long sum_AA = 0;
    for (int i = 0; i < 7; ++i)
        sum_AA += amino_acid_template_frequency[i];

    outFile << "#ALIGNMENT_POTENTIAL\n\n"
            << "PROT_LEN " << PROT_LEN << "\n"
            << "TEMPLATE_COUNT " << TEMPLATE_COUNT << "\n\n";
    for (int i = 0; i < PROT_LEN; ++i)
    {
        double score_sum = std::accumulate(alignment_position_score[i].begin(), alignment_position_score[i].end(), 0);

        outFile << std::left << std::setw(3) << i + 1;
        for (int j = 0; j < 7; ++j)
        {
            if (score_sum != 0) {
                alignment_potential[i][j] = (alignment_position_score[i][j] * alignment_position_sum[i] / score_sum +
                                             (double)amino_acid_template_frequency[j] * std::sqrt(alignment_position_sum[i]) / sum_AA)
                                            / (alignment_position_sum[i] + std::sqrt(alignment_position_sum[i]));
            }
            outFile << " " << std::fixed << std::setprecision(10)
                    << alignment_potential[i][j];
        }
        outFile << "\n";
    }

    outFile << "\n";
    for (int i = 0; i < PROT_LEN; ++i)
        for (int j = i + FRAG_MIN_LEN - 1; j < PROT_LEN; ++j)
            for (auto& fragment : matAMR[i][j])
                outFile << i << " " << j << " " << fragment << "\n";
}


void AlignmentThread(std::string &target, bool prcss_frgmnts, int s, int t)
{
    static int count_now = 0;
    for (int i = s; i < t; ++i)
    {
        StructureAlignment(target, db_protein_list[i], prcss_frgmnts);
        std::lock_guard<std::mutex> lock(al_mtx_print);
        progress_indicator("Aligning against protein database", ++count_now, db_protein_list.size());
    }
}


void DB_AlignMultiThrd(std::string target, bool prcss_frgmnts, int thread_num) {
    alignment_position_score.fill(std::array<double, 7> {});
    amino_acid_template_frequency.fill(0);
    alignment_position_sum.fill(0);

    DB_ListFiles();

    std::vector<std::thread> vecThreads;
    int size = db_protein_list.size();

    progress_indicator("Aligning against protein database", 0, 1);
    for (int i = 0; i < thread_num; ++i)
        vecThreads.push_back(std::thread(AlignmentThread, std::ref(target), prcss_frgmnts, i * size / thread_num, (i + 1) * (size) / thread_num));


    for (int i = 0; i < thread_num; ++i)
        vecThreads[i].join();

    Alignment_Eval(target);
}