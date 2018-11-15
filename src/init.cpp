#include "init.h"


// TARGET PROTEIN PROFILE
int PROT_LEN;

std::array<int, MAX_PROT_LEN> target_secondary_structure,
    target_solvent_accessibility;
std::array<std::array<double, MAX_PROT_LEN>, MAX_PROT_LEN> target_atom_distance;
std::string __target_sequence;
std::array<double, 20> target_AA_frequency;



// WEIGHTS
std::array<std::array<double, 3>, WGT_NUM> sc_weight;
std::array<std::array<double, 2>, 20> freq_bounds;




void Initialize(std::string target_path) {
    TargetProteinProfile(target_path); // MUST BE THE FIRST CALL

    dDFIRE_ReadLib(db_Scores + "dDFIRE.lib");

    GetScores();

    GetWeights();

    GetAlignments(target_path, true);
}


void TargetProteinProfile(std::string fName) {
    if (!FileExists(fName))
        throw std::runtime_error("ERR: Target PDB file not found!");

    progress_indicator("Target protein profile", 0, 1);

    // Solvent Accessibility
    std::vector<AminoAcid> amino_acids_vec;
    std::ifstream inFile(fName);
    std::string line;

    while (getline(inFile, line)) {
        if (line.size() >= 3 && line.substr(0, 3) == "TER")
            break;
        else if (line.size() < 47 || line.substr(0, 4) != "ATOM")
            continue;

        std::string atom_name = trim(line.substr(12, 4));
        if (atom_name != "CA" || !(line[16] == 'A' || line[16] == ' ' || line[16] == '1'))
            continue;

        amino_acids_vec.emplace_back(AminoAcid(
                                         line.substr(17, 3),
                                         std::stod(trim(line.substr(30, 8))),
                                         std::stod(trim(line.substr(38, 8))),
                                         std::stod(trim(line.substr(46, 8)))
                                     ));
        assert(residue_map.find(amino_acids_vec.back().name) !=
               residue_map.end());
        char sym = residue_map[amino_acids_vec.back().name];
        __target_sequence += sym;

        if (symbol_map.count(sym) > 0)
            target_AA_frequency[symbol_map[sym]]++;
    }
    PROT_LEN = amino_acids_vec.size();

    for (int i = 0; i < 20; ++i)
        target_AA_frequency[i] /= PROT_LEN;


    for (int i = 0; i < PROT_LEN; ++i)
        for (int j = i + 1; j < PROT_LEN; ++j)
        {
            double d = dist(amino_acids_vec[i].cords, amino_acids_vec[j].cords);
            target_atom_distance[i][j] = target_atom_distance[j][i] = d;

            if (d <= 14) {
                amino_acids_vec[i].neighbour_count++;
                amino_acids_vec[j].neighbour_count++;
            }
        }

    for (int i = 0; i < PROT_LEN; ++i)
        target_solvent_accessibility[i] = std::upper_bound(burial_classification.begin(),
                                          burial_classification.end(),
                                          amino_acids_vec[i].neighbour_count)
                                          - burial_classification.begin();

    // Secondary Structure
    std::string stdout, prev_line;
    while (system_call_err(ex_STRIDE + " -o " + fName + " 2>&1", stdout) != 0);
    //std::cerr << ex_STRIDE + " -o " + fName + " 2>&1" << "\n";
    //assert(stdout != "");
    std::stringstream ret(stdout);
    int POS = 0;

    while (getline(ret, line)) {
        if (line.substr(0, 3) == "STR")
        {
            for (int i = 10; i < 60; ++i)
            {
                if (isspace(prev_line[i]))
                    continue;
                if (isspace(line[i]))
                    line[i] = 'C';
                if (line[i] == 'b')
                    line[i] = 'B';

                assert(sec_structure_classification.find(line[i]) !=
                       sec_structure_classification.end());
                target_secondary_structure[POS++] = sec_structure_classification[line[i]];
            }
        }
        if (line.substr(0, 3) == "LOC")
            break;
        prev_line = line;
    }

    // assert(PROT_LEN == POS);
    matAMR.resize(PROT_LEN, std::vector<std::vector<std::string>>(PROT_LEN)); // alignment repo
    progress_indicator("Target protein profile", 1, 1);
}


void GetScores() {
    if (FileExists(db_Scores + "scores.txt"))
    {
        std::ifstream inFile(db_Scores + "scores.txt");
        std::string HEAD;
        char TYPE;

        inFile >> HEAD;
        for (int i = 0; i < 20; ++i)
        {
            progress_indicator("Burial & Pot & SS scores", i, 60);
            inFile >> TYPE;
            for (int j = 0; j < 7; ++j)
                inFile >> burial_potential[symbol_map[TYPE]][j];
        }

        inFile >> HEAD;
        for (int i = 0; i < 20; ++i)
        {
            progress_indicator("Burial & Pot & SS scores", i + 20, 80);
            inFile >> TYPE;
            for (int j = 0; j < 7; ++j)
                inFile >> structure_potential[symbol_map[TYPE]][j];
        }

        inFile >> HEAD;
        for (int i = 0; i < 20; ++i)
        {
            progress_indicator("Burial & Pot & SS scores", i + 40, 80);
            inFile >> TYPE >> pot_bar[symbol_map[TYPE]] >> pot_stdev[symbol_map[TYPE]];
        }

        inFile >> HEAD;
        for (int i = 0; i < 20; ++i)
        {
            progress_indicator("Burial & Pot & SS scores", i + 60, 80);
            inFile >> TYPE >> AA_freq_mean[symbol_map[TYPE]] >> AA_freq_sd[symbol_map[TYPE]];
        }

        progress_indicator("Burial & Pot & SS scores", 1, 1);
    } else {
        std::cout << "Burial & Pot & SS scores: " << BOLDRED << "NOT FOUND"
                  << RESET << std::endl;
        DB_PreprocMultiThrd();
    }
}


void GetWeights() {
    if (!FileExists(db_Scores + "weights.txt"))
        throw std::runtime_error("ERROR: Weights file not found!");

    std::ifstream inFile(db_Scores + "weights.txt");
    std::string HEAD;
    inFile >> HEAD;
    for (int i = 0; i < WGT_NUM; ++i)
    {
        inFile >> HEAD;
        for (int j = 0; j < 3; ++j)
            inFile >> sc_weight[i][j];
    }

    for (int i = 0; i < 20; ++i)
    {
        inFile >> HEAD;
        assert(symbol_map.count(HEAD[4]) > 0);
        int idx = symbol_map[HEAD[4]];
        inFile >> freq_bounds[idx][0] >> freq_bounds[idx][1];
    }
}


void GetAlignments(std::string target, bool prcss_frgmnts) {
    if (FileExists(db_Align + File_md5(target)))
    {
        progress_indicator("Laoding alignments", 0, 1);

        std::ifstream inFile(db_Align + File_md5(target));
        std::string HEAD;

        size_t size = FileSize(inFile);

        inFile >> HEAD;
        inFile >> HEAD >> PROT_LEN;
        inFile >> HEAD >> TEMPLATE_COUNT;

        for (int i = 0, p; i < PROT_LEN; ++i)
        {
            inFile >> p;
            for (int j = 0; j < 7; ++j)
                inFile >> alignment_potential[p - 1][j];

            progress_indicator("Laoding alignments", inFile.tellg(), size);
        }

        if (prcss_frgmnts) {
            for (int i, j; inFile >> i >> j;) {
                std::string fragment; inFile >> fragment;
                matAMR[i][j].emplace_back(fragment);
                progress_indicator("Laoding alignments", inFile.tellg(), size);
            }
        }
        progress_indicator("Laoding alignments", 1, 1);
    } else {
        std::cout << "Loading alignments: " << BOLDRED << "NOT FOUND"
                  << RESET << std::endl;
        DB_AlignMultiThrd(target, prcss_frgmnts);
    }
}