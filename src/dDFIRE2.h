#ifndef _DDFIRE2_H
#define _DDFIRE2_H

#include <string>
#include <array>
#include <unordered_map>
#include <fstream>
#include <vector>
#include "constants.h"
#include "mapping.h"
#include "utility.h"




class DFIRE2 {
public:
	DFIRE2(std::string libPath);
	void ReadLib(std::string libPath);
	double Calc_CFE(Protein& target);
private:
	int nTypes, nBins;
	bool libReady;
	std::unordered_map<std::string, int> atom_map;
	std::vector<std::vector<std::vector<double> > > edDFIRE; // 3D vector
};



#endif // _DDFIRE2_H
