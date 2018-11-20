#include "protein.h"

Protein::Protein(std::string fPath) : fPath( fPath )
{
    Parse_PDB(fPath, vecAmino_Acid);
    for (auto& amino_acid : vecAmino_Acid)
        sequence += amino_acid.symbol;
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