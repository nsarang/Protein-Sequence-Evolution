#include "db_process.h"


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




void Burial_PotStatistic(std::string fPath) {
    std::vector<AminoAcid> vecAminoAcid;
    std::array<double, 20> pot_frequency{}, pot_r, pot0, pot1{}, var0;

    ParsePDB(fPath, vecAminoAcid);
    int n = vecAminoAcid.size();

    for (auto& amino_acid : vecAminoAcid) {
        if (IsStandardAA(amino_acid.name))
            pot_frequency[symbol_map[residue_map[amino_acid.name]]]++;
    }

    for (int i = 0; i < 20; ++i)
    {
        if (pot_frequency[i] < 2)
            continue;
        pot_r[i] = std::exp(-std::sqrt(pot_frequency[i]) / (n * Pot_S_Constant));
        pot0[i] = (pot_frequency[i] * (pot_frequency[i] - 1) * pot_r[i] * (n - 1 / (1 - pot_r[i])))
                  / (n * (n - 1) * (1 - pot_r[i]));
        var0[i] = std::sqrt(std::pow(pot_r[i] * pot_frequency[i] * (1 - pot_frequency[i] / n), 2)
                            / ((1 - pot_r[i] * pot_r[i]) * n));
    }

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
        {
            if (dist(vecAminoAcid[i].cords, vecAminoAcid[j].cords) <= 14) {
                vecAminoAcid[i].neighbour_count++;
                vecAminoAcid[j].neighbour_count++;
            }

            if (!IsStandardAA(vecAminoAcid[i].name))
                continue;
            int idx = symbol_map[residue_map[vecAminoAcid[i].name]];
            pot1[idx] += (vecAminoAcid[i].name == vecAminoAcid[j].name)
                         * std::pow(pot_r[idx], j - i);
        }

    std::lock_guard<std::mutex> lock(db_mtx_1);
    for (auto& amino_acid : vecAminoAcid)
    {
        if (!IsStandardAA(amino_acid.name))
            continue;

        int class_num = std::upper_bound(burial_classification.begin(), burial_classification.end(),
                                         amino_acid.neighbour_count)
                        - burial_classification.begin();

        int idx = symbol_map[residue_map[amino_acid.name]];
        amino_acid_burial_count[idx][class_num]++;
        AA_frequency[idx]++;
    }

    for (int i = 0; i < 20; ++i)
    {
        if (pot_frequency[i] < 2)
            continue;
        vec_pot_score_db[i].push_back((pot1[i] - pot0[i]) / var0[i]);
        vec_AA_freq_db[i].push_back(pot_frequency[i] / n);
    }
}


void SecondaryStructure(std::string fPath) {
    std::string stdout, line, prev_line;
    while (system_call_err(ex_STRIDE + " -o " + fPath + " 2>&1", stdout) != 0);
    std::stringstream ret(stdout);

    while (getline(ret, line)) {
        if (line.substr(0, 3) == "STR")
        {
            std::lock_guard<std::mutex> lock(db_mtx_2);
            for (int i = 10; i < 60; ++i)
            {
                if (isspace(prev_line[i]) || prev_line[i] == 'X')
                    continue;
                if (isspace(line[i]))
                    line[i] = 'C';
                if (line[i] == 'b')
                    line[i] = 'B';

                assert(sec_structure_classification.count(line[i]) > 0);
                amino_acid_structure_count[symbol_map[prev_line[i]]][sec_structure_classification[line[i]]]++;
            }
        }
        if (line.substr(0, 3) == "LOC")
            return;

        prev_line = line;
    }
}


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
