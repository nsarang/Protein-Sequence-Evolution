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
	
public:
	Protein(std::string fPath);
	static void Parse_PDB(std::string fPath, std::vector<AminoAcid>& retVec);
	double CA_Atom_Distance(int i, int j);
	int length();
	char operator[](int i);

private:
	void Calculate_Distances();
	std::string sequence,
	    fPath,
	    md5;
	std::vector<AminoAcid> vecAmino_Acid;
	std::vector<std::vector<double > > vecAtom_Distance;
	std::vector<int> aSolvent_Accessibility, sSecondary_Structure;
    

};




#endif // __PROTEIN_H__
