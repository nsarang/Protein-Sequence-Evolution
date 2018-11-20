#ifndef __PROTEIN_PROFILE_H__
#define __PROTEIN_PROFILE_H__


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
#include <functional>
#include <dirent.h>
#include <unistd.h>
#include <tuple>
#include "mapping.h"
#include "utility.h"
#include "constants.h"
#include "statistics.h"
#include "protein.h"




class ProteinProfile {
public:
	ProteinProfile(Protein target, std::vector<std::string> vecDatabase);
	void CalculateProfiles(bool bAlgn, bool bSolvent, bool bPot,
	                       bool bSS, bool bSave_Frags, int potS_Param);
	void Find_Homologous_Proteins(std::vector<std::string> vecDB,
	                              double dCutOff, double bVerbose);

private:
	template<class FuncType, class vecType>
	void Thread_Manager(std::vector<std::function<FuncType> > vecFunctions,
	                    std::vector<vecType>& vecDatabase, bool bVerbose = false,
	                    std::string sMsg = "", int nThreads = 12);

	template<class FuncType, class vecType>
	void Processing_Thread(std::vector<std::function<FuncType> > vecFunctions,
	                       std::vector<vecType> vecBatch, int& nCount_Now);


	void Calculate_Solvent_Profile();
	void Calculate_SS_Profile();
	void Calculate_Pot_AAFreq_Profile();
	void Calculate_Alignment_Profile(bool bSave_Frags, int nGap_Score = GAP_SCORE,
	                                 int nDist_CutOff = DIST_CUTOFF, int nMin_Frag = FRAG_MIN_LEN);
	void Process_IsHomologue(std::tuple<std::string>& tArgs);
	void Process_Solvent(std::tuple<Protein>& tArgs);
	void Process_Pot_AAFreq(std::tuple<Protein>& tArgs);
	void Process_SS(std::tuple<Protein>& tArgs);


	Protein _refProtein;
	bool _dScore_CutOff;
	std::vector<std::string> _vecHomologous_Proteins;
	std::array<std::array<double, 7>, 20> _aSolvent_Profile, _aSec_Profile;
	std::vector<std::array<double, 7> > _aAlgn_Profile;
	std::vector<std::vector<std::vector<std::string> > > _matFragments;
	std::vector<std::tuple<double, std::string, std::string, std::vector<double> >> _vecTupleAlignments;
	int _potS_Param;


	std::array<std::array<long long, 7>, 20> _aSolvent_AA_Count, _aSec_AA_Count;
	std::array<std::vector<double>, 20> _vecPotScores, _vecAA_Freqs;
	std::array<long long, 20> _aAA_Total_Count;
	std::vector<long long> _aAlgn_Position_Count;
	std::array<long long, 7> _aAlgn_Class_Count;

	std::mutex _mtx_count, _mtx_SS, _mtx_solvent, _mtx_pot;
};




#endif // __PROTEIN_PROFILE_H__
