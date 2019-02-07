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
#include <unistd.h>
#include <optional>
#include <tuple>
#include "mapping.h"
#include "utility.h"
#include "constants.h"
#include "statistics.h"
#include "protein.h"




class ProteinProfile {
	friend class Evaluator;
	friend class DeepAI;

public:
	ProteinProfile(Protein target,
	               double dAlgn_Score_CutOff = ALIGN_SC_CUTTOFF,
	               double dFrag_Score_Cutoff = DIST_CUTOFF,
	               double dGap_Penalty = GAP_PENALTY,
	               double dPotS_Param = Pot_S_Constant,
	               int nMin_Frag = FRAG_MIN_LEN);

	void Find_Homologous_Proteins(std::vector<std::string> vecDB,
	                              double bVerbose = false);

	void CalculateProfiles(int nFlag = 1 + 2 + 4 + 8,
	                       bool bVerbose = false,
	                       bool bSave_Frags = true);

	bool CalculateRemainingProfiles(std::optional<std::vector<std::string>> vecDB,
	                                bool bVerbose = false);

	void Read_FromFile(std::string sDirectory = db_Profiles);
	void Write_ToFile(bool bWriteCounts = false, std::string sDirectory = db_Profiles);
	std::string QuickInfo(bool bIncludeAlignInfo = false);
	bool FragmentExists(int i, int j);
	std::string FragmentFetch(int i, int j);


private:
	template<class FuncType, class vecType>
	void Thread_Manager(std::vector<std::function<FuncType> > vecFunctions,
	                    std::vector<vecType>& vecDatabase,
	                    bool bVerbose = false, std::string sMsg = "",
	                    int nThreads = std::thread::hardware_concurrency());

	template<class FuncType, class vecType>
	void Processing_Thread(std::vector<std::function<FuncType> > vecFunctions,
	                       std::vector<vecType> vecBatch, int& nCount_Now);


	void Calculate_Solvent_Profile();
	void Calculate_SS_Profile();
	void Calculate_Pot_AAFreq_Profile();
	void Calculate_Alignment_Profile(bool bSave_Frags, double dFrag_Score_Cutoff,
	                                 double dGap_Penalty, int nMin_Frag);
	void Process_IsHomologue(std::string fPath);
	void Process_Solvent(Protein& target);
	void Process_Pot_AAFreq(Protein& target);
	void Process_SS(Protein& target);
	std::string RelativeFileName(std::string sProfile_Name);


	Protein _refProtein; // Target protein
	double _dAlgn_Score_CutOff, _dFrag_Score_Cutoff, _dGap_Penalty; // Alignment parameters
	double _dPotS_Param; // Pot parameter
	int _nMin_Frag;

	bool bAlgn_Rdy{ false }, bFrags_Rdy{ false }, bSolvent_Rdy{ false },
	     bSS_Rdy{ false }, bPot_Rdy { false };

	std::vector<Protein> _vecHomologous_Proteins; // Protein with homologous structures extracted from CATH
	std::array<std::array<double, 7>, 20> _aSolvent_Profile{ }, _aSec_Profile{ };
	std::array<double, 20> _aPot_Bar{ }, _aPot_Stdev{ }, _aAA_Freq_Mean{ }, _aAA_Freq_Stdev{ };
	std::vector<std::array<double, 7> > _aAlgn_Profile;

	std::vector<std::vector<std::vector<std::string> > > _matFragments; // Alignment fragments
	std::vector<std::tuple<double, std::string, std::string, std::vector<double> >> _vecTupleAlignments;

	// Objects used in preprocessing
	std::array<std::array<long long, 7>, 20> _aSolvent_AA_Count{ }, _aSec_AA_Count{ };
	std::array<std::vector<double>, 20> _vecPotScores, _vecAA_Freqs;
	std::array<long long, 20> _aAA_Total_Count{ };
	std::vector<long long> _aAlgn_Position_Count{ };
	std::array<long long, 7> _aAlgn_Class_Count{ };
	inline static std::mutex _mtx_count, _mtx_SS, _mtx_solvent, _mtx_pot, _mtx_algn; // Mutexs used in threads
};




#endif // __PROTEIN_PROFILE_H__