#ifndef __PROTEIN_H__
#define __PROTEIN_H__


#include <cmath>
#include <array>
#include <vector>
#include <thread>
#include <mutex>
#include <unordered_map>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <dirent.h>
#include "mapping.h"
#include "utility.h"
#include "constants.h"
#include "statistics.h"



struct AminoAcid {
	AminoAcid(std::string name, char symbol, double x = 0, double y = 0, double z = 0)
		: name( name ), symbol( symbol ), cords( std::tie(x, y, z)) {}

	std::string name;
	char symbol;
	std::tuple<double, double, double> cords;
};


class Protein {
	friend class ProteinProfile;
	friend class Evaluator;
	friend class DeepAI;

public:
	Protein(std::string fPath, int nFlag = 0, double dPotS_Param = Pot_S_Constant,
	        bool calcMD5 = true);
	static void Parse_PDB(std::string fPath, std::vector<AminoAcid>& retVec);

	double CA_Atom_Distance(int i, int j);
	std::string Get_Sequence();
	int length();
	char operator[](int i);

	void Calculate_Distances(bool bForceCalc = false);
	void Calculate_Solvent(bool bForceCalc = false);
	void Calculate_SS(bool bForceCalc = false);
	void Calculate_Pot(double dPotS_Param, bool bForceCalc = false);

private:
	static double dist3D(std::tuple<double, double, double>&, std::tuple<double, double, double> &t);
	static bool IsStandardAA(std::string abrv);
	static bool IsStandardAA(char symbol);

	std::string sequence,
	    fPath,
	    md5;

	bool bDist_Rdy{ false }, bSolv_Rdy{ false }, bSS_Rdy{ false }, bPot_Rdy{ false };

	std::vector<AminoAcid> vecAmino_Acid;
	std::vector<std::vector<double > > vecAtom_Distance;
	std::vector<int> aSolvent_Accessibility, aSecondary_Structure;
	std::array<double, 20> aPot_Values{ }, aAA_Freqs{ };
};




#endif // __PROTEIN_H__
