#include "protein.h"



Protein::Protein(std::string fPath, bool bDist, bool bSolv, bool bSec)
    : fPath( fPath ), md5( File_md5(fPath) ), bDist_Rdy( bDist ), bSolv_Rdy( bSolv ), bSS_Rdy( bSec )
{
    Parse_PDB(fPath, vecAmino_Acid);
    for (auto& amino_acid : vecAmino_Acid)
        sequence += amino_acid.symbol;
    
    if (bDist)
        Calculate_Distances();
    if (bSolv)
        Calculate_Solvent();
    if (bSec)
        Calculate_SS();
}


void Protein::Parse_PDB(std::string fPath, std::vector<AminoAcid> &retVec) {
    retVec.clear();
    std::ifstream inFile(fPath);
    std::string line;

    while (getline(inFile, line)) {
        if (line.size() < 47 || line.substr(0, 4) != "ATOM")
            continue;

        std::string atom_name = trim(line.substr(12, 4));
        if (atom_name != "CA" || !(line[16] == 'A' || line[16] == ' ' || line[16] == '1'))
            continue;

        retVec.emplace_back(AminoAcid(
                                resName_to_sym.count(line.substr(17, 3)) ?  line.substr(17, 3) : "UNK",
                                ' ',
                                std::stod(trim(line.substr(30, 8))),
                                std::stod(trim(line.substr(38, 8))),
                                std::stod(trim(line.substr(46, 8)))
                            ));

        retVec.back().symbol = resName_to_sym[retVec.back().name];
    }
}


int Protein::length() {
    return sequence.length();
}


double Protein::CA_Atom_Distance(int i, int j) {
    return vecAtom_Distance[i][j];
}


char Protein::operator[](int i) {
    return sequence[i];
}


void Protein::Calculate_Distances() {
    for (int i = 0; i < vecAmino_Acid.size(); ++i) {
        for (int j = i + 1; j < vecAmino_Acid.size(); ++j)
        {
            double d = dist(vecAmino_Acid[i].cords, vecAmino_Acid[j].cords);
            vecAtom_Distance[i][j] = vecAtom_Distance[j][i] = d;
        }
    }
}


void Protein::Calculate_Solvent() {
    if (bDist_Rdy == false) {
        Calculate_Distances();
        bDist_Rdy = true;
    }
    bSolv_Rdy = true;
    int n = sequence.length();
    std::vector<int> neighbours( n );

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (CA_Atom_Distance(i, j) <= 14) {
                neighbours[i]++;
                neighbours[j]++;
            }

    for (int i = 0; i < n; ++i) {
        int class_num = std::upper_bound(solvent_classes.begin(), solvent_classes.end(),
                                         neighbours[i]) - solvent_classes.begin();

        aSolvent_Accessibility[i] = class_num;
    }
}


void Protein::Calculate_SS() {
    std::string stdout, line, prev_line;

    while (system_call_err(ex_STRIDE + " -o " + fPath + " 2>&1", stdout) != 0);
    assert(stdout != "");
    std::stringstream ret( stdout );

    while (getline(ret, line)) {
        if (line.substr(0, 3) == "STR") {
            for (int i = 10; i < 60; ++i) {
                if (isspace(prev_line[i]) || prev_line[i] == 'X')
                    continue;
                if (isspace(line[i]))
                    line[i] = 'C';
                if (line[i] == 'b')
                    line[i] = 'B';

                assert(sec_classes.count(line[i]) != 0);
                sSecondary_Structure[i] = sec_classes[ line[i] ];
            }
        }
        if (line.substr(0, 3) == "LOC")
            return;

        prev_line = line;
    }
}
