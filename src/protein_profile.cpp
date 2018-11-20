#include "protein_profile.h"
using namespace std::placeholders;


// DATABASE FILE LIST
std::vector<std::string> db_protein_list;


// MULTI-THREADING
std::mutex db_mtx_1, db_mtx_2, db_mtx_print;


// SCORE PROFILES
std::array<std::array<double, 7>, 20> burial_potential,
    structure_potential;


// POT STATISTICS
std::array<double, 20> pot_bar, pot_stdev;
std::vector<double> vec_pot_score_db[20];


// AA FREQUENCY
std::array<long long, 20> AA_frequency;
std::array<double, 20> AA_freq_mean, AA_freq_sd;
std::vector<double> vec_AA_freq_db[20];


// PREPROCESS ARRAYS
std::array<std::array<long long, 7>, 20> amino_acid_burial_count,
    amino_acid_structure_count;



void ProteinProfile::CalculateProfiles(bool bAlgn, bool bSolvent, bool bPot, bool bSS,
                                       bool bSave_Frags, int potS_Param) {
    std::vector<std::function<void(Protein&)> > vecProcFuncs;
    std::vector<std::function<void()> > vecCalcFuncs;

    if (bAlgn) {
        vecCalcFuncs.push_back([ = ] { Calculate_Alignment_Profile(bSave_Frags); });
    }
    if (bSolvent) {
        vecProcFuncs.push_back(std::bind(&ProteinProfile::Process_Solvent, this, _1));
        vecCalcFuncs.push_back(std::bind(&ProteinProfile::Calculate_Solvent_Profile, this));
    }
    if (bPot) {
        vecProcFuncs.push_back(std::bind(&ProteinProfile::Process_Pot_AAFreq, this, _1));
        vecCalcFuncs.push_back(std::bind(&ProteinProfile::Calculate_Pot_AAFreq_Profile, this));
    }
    if (bSS) {
        vecProcFuncs.push_back(std::bind(&ProteinProfile::Process_SS, this, _1));
        vecCalcFuncs.push_back(std::bind(&ProteinProfile::Calculate_SS_Profile, this));
    }

    Thread_Manager(vecProcFuncs, _vecHomologous_Proteins);
    for (auto& CalcFunction : vecCalcFuncs)
        CalcFunction();
}


void ProteinProfile::Find_Homologous_Proteins(std::vector<std::string> vecDB,
        double dCutOff, double bVerbose)
{
    _vecHomologous_Proteins.clear();
    _vecTupleAlignments.clear();

    _dScore_CutOff = dCutOff;

    std::vector<std::function<void(std::string&)> > vecProcs;
    vecProcs.push_back(std::bind(&ProteinProfile::Process_IsHomologue, this, _1));
    Thread_Manager(vecProcs, vecDB, bVerbose, "Finding homologous proteins");
}


template<class FuncType, class Args>
void ProteinProfile::Thread_Manager(std::vector<std::function<FuncType> > vecFuncs,
                                    std::vector<Args>& vecDB, bool bVerbose,
                                    std::string sMsg, int nThreads)
{
    std::vector<std::thread> vecThreads;
    int size = vecDB.size(),
        nCount_Now = 0;

    for (int i = 0; i < nThreads; ++i) {
        std::vector vecBatch(vecDB.begin() + i * size / nThreads,
                             vecDB.begin() + (i + 1) * (size) / nThreads);
        vecThreads.push_back(std::thread([&] { Processing_Thread(vecFuncs, vecBatch, nCount_Now); }));
    }

    if (bVerbose) {
        Progress_Indicator(sMsg, 0, size);
        while (true) {
            Progress_Indicator(sMsg, nCount_Now, size);
            if (nCount_Now == size)
                break;
            usleep(250);
        }
    }

    for (int i = 0; i < nThreads; ++i)
        vecThreads[i].join();
}


template<class FuncType, class Args>
void ProteinProfile::Processing_Thread(std::vector<std::function<FuncType> > vecFuncs,
                                       std::vector<Args> vecBatch, int& nCount_Now)
{
    for (auto data : vecBatch) {
        for (auto ProcFunction : vecFuncs)
            ProcFunction(data);
        std::lock_guard<std::mutex> lock(_mtx_count);
        nCount_Now++;
    }
}


void ProteinProfile::Process_IsHomologue(std::string fPath) {
    std::string stdout, line, algn_line1, algn_line2;

    while (system_call_err(ex_TMALIGN + " " + _refProtein.fPath + "  " + fPath + " -a 2>&1", stdout) != 0);
    assert(stdout != "");
    std::stringstream ret(stdout);

    double score;
    getline(ret, line); getline(ret, line);

    ret >> line >> score;
    if (score > _dScore_CutOff)
        return;

    ret >> algn_line1 >> algn_line2;
    assert(algn_line1.size() == algn_line2.size());
    int align_len = algn_line1.size();

    std::vector<double> atom_dist( align_len + 1 );
    for (int i = 1; i <= align_len; ++i)
        ret >> atom_dist[i];

    _vecTupleAlignments.emplace_back(score, algn_line1, algn_line2, atom_dist);
    _vecHomologous_Proteins.emplace_back(Protein(fPath));
}


void ProteinProfile::Process_Solvent(Protein& target) {
    int n = target.length();
    std::vector<int> neighbours( n );

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
        {
            if (target.CA_Atom_Distance(i, j) <= 14) {
                neighbours[i]++;
                neighbours[j]++;
            }
        }

    for (int i = 0; i < n; ++i) {
        auto& amino_acid = target.vecAmino_Acid[i];
        if (!IsStandardAA(amino_acid.name))
            continue;

        int class_num = std::upper_bound(solvent_classes.begin(), solvent_classes.end(),
                                         neighbours[i]) - solvent_classes.begin();

        std::lock_guard<std::mutex> lock(_mtx_solvent);
        int idx = sym_to_idx[amino_acid.symbol];
        _aSolvent_AA_Count[idx][class_num]++;
        _aAA_Total_Count[idx]++;
    }
}


void ProteinProfile::Process_Pot_AAFreq(Protein& target) {
    std::array<double, 20> AA_freq{}, pot_r, pot0, pot1{}, var0;
    int n = target.length();

    for (auto& amino_acid : target.vecAmino_Acid) {
        if (IsStandardAA(amino_acid.name))
            AA_freq[sym_to_idx[amino_acid.symbol]]++;
    }

    for (int i = 0; i < 20; ++i) {
        if (AA_freq[i] < 2)
            continue;
        pot_r[i] = std::exp(-std::sqrt(AA_freq[i]) / (n * _potS_Param));
        pot0[i] = (AA_freq[i] * (AA_freq[i] - 1) * pot_r[i] * (n - 1 / (1 - pot_r[i])))
                  / (n * (n - 1) * (1 - pot_r[i]));
        var0[i] = std::sqrt(std::pow(pot_r[i] * AA_freq[i] * (1 - AA_freq[i] / n), 2)
                            / ((1 - pot_r[i] * pot_r[i]) * n));
    }

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
        {
            if (!IsStandardAA(target.vecAmino_Acid[i].name))
                continue;
            int idx = sym_to_idx[target.vecAmino_Acid[i].symbol];
            pot1[idx] += (target.vecAmino_Acid[i].symbol == target.vecAmino_Acid[j].symbol)
                         * std::pow(pot_r[idx], j - i);
        }

    std::lock_guard<std::mutex> lock(_mtx_pot);
    for (int i = 0; i < 20; ++i)
    {
        if (AA_freq[i] < 2)
            continue;
        _vecPotScores[i].push_back((pot1[i] - pot0[i]) / var0[i]);
        _vecAA_Freqs[i].push_back(AA_freq[i] / n);
    }
}


void ProteinProfile::Process_SS(Protein& target) {
    std::string stdout, line, prev_line;

    while (system_call_err(ex_STRIDE + " -o " + target.fPath + " 2>&1", stdout) != 0);
    assert(stdout != "");
    std::stringstream ret( stdout );

    while (getline(ret, line)) {
        if (line.substr(0, 3) == "STR")
        {
            std::lock_guard<std::mutex> lock(_mtx_SS);
            for (int i = 10; i < 60; ++i)
            {
                if (isspace(prev_line[i]) || prev_line[i] == 'X')
                    continue;
                if (isspace(line[i]))
                    line[i] = 'C';
                if (line[i] == 'b')
                    line[i] = 'B';

                assert(sec_classes.count(line[i]) != 0);
                _aSec_AA_Count[ sym_to_idx[prev_line[i]] ][ sec_classes[line[i]] ]++;
            }
        }
        if (line.substr(0, 3) == "LOC")
            return;

        prev_line = line;
    }
}


void ProteinProfile::Calculate_Alignment_Profile(bool bSave_Frags, int nGap_Score, int nDist_CutOff, int nMin_Frag) {
    for (auto&[score, seq1, seq2, atom_dist] : _vecTupleAlignments) {

        int align_len = seq1.length();
        std::vector<double> sum_dist( align_len + 1 );

        for (int i = 1; i <= align_len; ++i)
            sum_dist[i] = (atom_dist[i] == nGap_Score ? GAP_PENALTY : atom_dist[i])
                         + sum_dist[i - 1];


        for (int i = 0, pos = 0; i < align_len; ++i) {
            if (seq1[i] == '-')
                continue;

            pos++;
            if (seq2[i] == '-' || seq2[i] == 'X')
                continue;

            assert(algn_classes.count(seq2[i]) != 0);
            int nClass = algn_classes[seq2[i]];

            _aAlgn_Profile[pos - 1][nClass] += 10 - std::min(10., atom_dist[i + 1]); // Convert to similarity
            _aAlgn_Position_Count[pos - 1]++;
            _aAlgn_Class_Count[nClass]++;
        }

        if (bSave_Frags) {
            int pos = 0;
            for (int i = 0; i < align_len; ++i) {
                for (int l = nMin_Frag; i + l <= align_len; ++l) {
                    if (atom_dist[i + 1] == 0 || atom_dist[i + l] == 0)
                        continue;
                    double sc = sum_dist[i + l] - sum_dist[i];
                    assert(sc > 0);
                    if (sc <= l * nDist_CutOff) // hit
                    {
                        std::string fragment;
                        int nFragLen = 0;
                        for (int j = i; j < i + l; ++j) {
                            if (seq1[j] == '-') // skip gaps on seq1
                                continue;
                            fragment += seq2[j];
                            nFragLen++;
                        }
                        _matFragments[pos][pos + nFragLen - 1].push_back(fragment);
                    }
                }
                if (seq1[i] != '-')
                    pos++;
            }
        }
    }
}


/*
void BurSecPot_Eval() {
    std::ofstream outFile(db_Scores + "scores.txt"),
        outFile2(db_Scores + "counts.txt");

    std::array<long long, 7> burial_count{},
        structure_count{};

    long long AA_sum = std::accumulate(AA_frequency.begin(), AA_frequency.end(), 0);

    for (int i = 0; i < 7; ++i)
    {
        for (int j = 0; j < 20; ++j)
        {
            burial_count[i] += amino_acid_burial_count[j][i];
            structure_count[i] += amino_acid_structure_count[j][i];
        }
    }


    outFile << "#BURIAL_POTENTIAL\n\n";
    for (int i = 0; i < 20; ++i)
    {
        outFile << index_to_symbol[i];
        for (int j = 0; j < 7; ++j)
        {
            burial_potential[i][j] = (double)amino_acid_burial_count[i][j] / (burial_count[j]
                                     * (AA_frequency[i] / (double)AA_sum));
            outFile << " " << std::fixed << std::setprecision(10)
                    << burial_potential[i][j];
        }
        outFile << "\n";
    }


    outFile << "\n\n\n#STRUCTURE_POTENTIAL\n\n";
    for (int i = 0; i < 20; ++i)
    {
        outFile << index_to_symbol[i];
        for (int j = 0; j < 7; ++j)
        {
            structure_potential[i][j] = (double)amino_acid_structure_count[i][j] / (structure_count[j]
                                        * (AA_frequency[i] / (double)AA_sum));
            outFile << " " << std::fixed << std::setprecision(10)
                    << structure_potential[i][j];
        }
        outFile << "\n";
    }


    outFile << "\n\n\n#POT_STATISTIC\n\n";
    for (int i = 0; i < 20; ++i)
    {
        vec_pot_score_db[i] = OutlierElimination_IQR(vec_pot_score_db[i]);
        // std::ofstream outRes(std::string() + "PotScore-" + index_to_symbol[i] + ".txt");
        // for (auto s : pot_db_score[i])
        //     outRes << s << "\n";
        Mean_SD(vec_pot_score_db[i], pot_bar[i], pot_stdev[i]);
        outFile << index_to_symbol[i]
                << " " << std::fixed << std::setprecision(10)
                << pot_bar[i] << " " << pot_stdev[i] << "\n";
    }


    outFile << "\n\n\n#AA_FREQUENCY\n\n";
    for (int i = 0; i < 20; ++i)
    {
        vec_AA_freq_db[i] = OutlierElimination_IQR(vec_AA_freq_db[i]);
        // std::ofstream outRes(std::string() + "AA_FREQ-" + index_to_symbol[i] + ".txt");
        // for (auto s : AA_freq_db[i])
        //     outRes << s << "\n";
        Mean_SD(vec_AA_freq_db[i], AA_freq_mean[i], AA_freq_sd[i]);
        outFile << index_to_symbol[i]
                << " " << std::fixed << std::setprecision(10)
                << AA_freq_mean[i] << " " << AA_freq_sd[i] << "\n";
    }


    outFile2 << "#BURIAL_COUNT\n\n";
    for (int i = 0; i < 20; ++i)
    {
        outFile2 << index_to_symbol[i];
        for (int j = 0; j < 7; ++j)
            outFile2 << " " << std::setw(10) << amino_acid_burial_count[i][j];
        outFile2 << "\n";
    }


    outFile2 << "\n\n#STRUCTURE_COUNT\n\n";
    for (int i = 0; i < 20; ++i)
    {
        outFile2 << index_to_symbol[i];
        for (int j = 0; j < 7; ++j)
            outFile2 << " " << std::setw(10) << amino_acid_structure_count[i][j];
        outFile2 << "\n";
    }

    outFile2 << "\n\n#AA_COUNT\n\n";
    for (int i = 0; i < 20; ++i)
    {
        outFile2 << index_to_symbol[i]
                 << " " << AA_frequency[i] << "\n";
    }
}


void DB_ListFiles(std::string path) {
    if (db_protein_list.empty()) {
        DIR *hDir;
        dirent *hFile;
        assert(hDir = opendir(path.c_str()));
        int total = stoi(system_call("ls " + path + " | wc -l"));

        db_protein_list.reserve(total);
        while ((hFile = readdir(hDir))) {
            std::string fName = hFile->d_name;
            if (fName.size() != 7 || !isdigit(fName[0]))
                continue;

            db_protein_list.push_back(path + fName);
        }
    }
}


void PreprocessThread(int s, int t)
{
    static int count_now = 0;
    for (int i = s; i < t; ++i)
    {
        Burial_PotStatistic(db_protein_list[i]);
        SecondaryStructure(db_protein_list[i]);

        std::lock_guard<std::mutex> lock(db_mtx_print);
        progress_indicator("Calculating Burial & Pot & SS scores", ++count_now, db_protein_list.size());
    }
}


void DB_PreprocMultiThrd(int thread_num) {
    for (int i = 0; i < 20; ++i) {
        vec_pot_score_db[i].clear();
        vec_AA_freq_db[i].clear();
    }
    amino_acid_burial_count.fill(std::array<long long, 7> {});
    amino_acid_structure_count.fill(std::array<long long, 7> {});
    AA_frequency.fill(0);


    DB_ListFiles();

    std::vector<std::thread> vecThreads;
    int size = db_protein_list.size();

    progress_indicator("Calculating Burial & Pot & SS scores", 0, 1);
    for (int i = 0; i < thread_num; ++i)
        vecThreads.push_back(std::thread(PreprocessThread, i * size / thread_num, (i + 1) * (size) / thread_num));


    for (int i = 0; i < thread_num; ++i)
        vecThreads[i].join();

    BurSecPot_Eval();
}
*/
