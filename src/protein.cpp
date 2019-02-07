#include "protein.h"

namespace sp = subprocess;



Protein::Protein(std::string fPath, int nFlag, double dPotS_Param)
    : fPath( fPath )
{
    assert(utility::FileExists(fPath));
    md5 = utility::File_md5(fPath);
    Parse_PDB(fPath, vecAmino_Acid);
    for (auto& amino_acid : vecAmino_Acid)
        sequence += amino_acid.symbol;

    if (nFlag & 1)
        Calculate_Distances();
    if (nFlag & 2)
        Calculate_Solvent();
    if (nFlag & 4)
        Calculate_SS();
    if (nFlag & 8)
        Calculate_Pot(dPotS_Param);
}


void Protein::Parse_PDB(std::string fPath, std::vector<AminoAcid> &retVec) {
    retVec.clear();
    std::ifstream inFile(fPath);
    std::string line, last, current;

    while (getline(inFile, line)) {
        if (line.size() < 47 || line.substr(0, 4) != "ATOM")
            continue;

        std::string atom_name = utility::trim(line.substr(12, 4));
        if (atom_name != "CA" || (last == line.substr(22, 5)) )
            continue;

        last = line.substr(22, 5);
        retVec.emplace_back(AminoAcid(
                                resName_to_sym.count(line.substr(17, 3)) ?  line.substr(17, 3) : "UNK",
                                ' ',
                                std::stod(utility::trim(line.substr(30, 8))),
                                std::stod(utility::trim(line.substr(38, 8))),
                                std::stod(utility::trim(line.substr(46, 8)))
                            ));

        retVec.back().symbol = resName_to_sym[retVec.back().name];
    }
}


std::string Protein::Get_Sequence() {
    return sequence;
}


int Protein::length() {
    return sequence.length();
}


double Protein::CA_Atom_Distance(int i, int j) {
    if (bDist_Rdy == false)
        throw;
    return vecAtom_Distance[i][j];
}


char Protein::operator[](int i) {
    return sequence[i];
}


void Protein::Calculate_Distances(bool bForceCalc) {
    if (bDist_Rdy && !bForceCalc)
        return;

    vecAtom_Distance.resize( length(), std::vector<double>(length()) );

    for (int i = 0; i < length(); ++i) {
        for (int j = i + 1; j < length(); ++j)
        {
            double d = dist3D(vecAmino_Acid[i].cords, vecAmino_Acid[j].cords);
            vecAtom_Distance[i][j] = vecAtom_Distance[j][i] = d;
        }
    }
    bDist_Rdy = true;
}


void Protein::Calculate_Solvent(bool bForceCalc) {
    if (bSolv_Rdy && !bForceCalc)
        return;
    if (bDist_Rdy == false)
        Calculate_Distances();

    int n = length();
    aSolvent_Accessibility.resize( n );
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
    bSolv_Rdy = true;
}


void Protein::Calculate_SS(bool bForceCalc) {
    if (bSS_Rdy && !bForceCalc)
        return;
    aSecondary_Structure.resize( length() );
    std::string stdout, line, prev_line;

    auto oBuffer = subprocess::check_output({ex_STRIDE.c_str(), "-o", fPath.c_str()}, sp::error{sp::PIPE});
   // std::ofstream ott("out.txt");
   // ott << oBuffer.buf.data();
    std::stringstream ret(oBuffer.buf.data());
   // std::cerr << ret.str() << "\n";

    int pos = 0;
    while (getline(ret, line)) {
        if (line.substr(0, 3) == "STR") {
            for (int i = 10; i < 60; ++i) {
 //               std::cerr << prev_line << "\n" << line << "\n";

                if (isspace(prev_line[i]))
                    continue;
                if (isspace(line[i]))
                    line[i] = 'C';
                if (line[i] == 'b')
                    line[i] = 'B';

                assert(isalpha(line[i]));
                assert(sec_classes.count(line[i]) != 0);
                aSecondary_Structure[pos++] = sec_classes[ line[i] ];
            }
        }
        if (line.substr(0, 3) == "LOC")
            break;

        prev_line = line;
    }
    assert(pos == length());
    bSS_Rdy = true;
}


void Protein::Calculate_Pot(double dPotS_Param, bool bForceCalc) {
    if (bPot_Rdy && !bForceCalc)
        return;

    std::array<double, 20> AA_freq{0}, pot_r, pot0, pot1{}, var0;
    int n = length();

    for (auto c : sequence)
        if (IsStandardAA(c))
            AA_freq[ sym_to_idx[c] ]++;


    for (int i = 0; i < 20; ++i) {
        if (AA_freq[i] < 2)
            continue;
        pot_r[i] = std::exp(-std::sqrt(AA_freq[i]) / (n * dPotS_Param));
        pot0[i] = (AA_freq[i] * (AA_freq[i] - 1) * pot_r[i] * (n - 1 / (1 - pot_r[i])))
                  / (n * (n - 1) * (1 - pot_r[i]));
        var0[i] = std::sqrt(std::pow(pot_r[i] * AA_freq[i] * (1 - AA_freq[i] / n), 2)
                            / ((1 - pot_r[i] * pot_r[i]) * n));
    }

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            if (!IsStandardAA(sequence[i]))
                continue;
            int idx = sym_to_idx[ sequence[i] ];
            pot1[idx] += (sequence[i] == sequence[j]) * std::pow(pot_r[idx], j - i);
        }

    for (int i = 0; i < 20; ++i) {
        if (AA_freq[i] < 2)
            continue;

        aPot_Values[i] = (pot1[i] - pot0[i]) / var0[i];
        aAA_Freqs[i] = AA_freq[i] / n;
    }

    bPot_Rdy = true;
}


double Protein::dist3D(std::tuple<double, double, double> &t1, std::tuple<double, double, double> &t2) {
    double dx = std::get<0>(t1) - std::get<0>(t2),
           dy = std::get<1>(t1) - std::get<1>(t2),
           dz = std::get<2>(t1) - std::get<2>(t2);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}


bool Protein::IsStandardAA(std::string abrv) {
    return !(abrv == "ASX" || abrv == "GLX" || abrv == "SEC" || abrv == "PYL" || abrv == "UNK");
}


bool Protein::IsStandardAA(char symbol) {
    return !(symbol == 'B' || symbol == 'Z' || symbol == 'U' || symbol == 'O' || symbol == 'X');
}
