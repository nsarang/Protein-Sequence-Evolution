#ifndef __EVALUATOR_H__
#define __EVALUATOR_H__


#include <iostream>
#include "mapping.h"
#include "constants.h"
#include "protein.h"
#include "protein_profile.h"
#include "dDFIRE2.h"
#include "utility.h"



class Evaluator {
public:
	Evaluator(std::string fPath);
	Evaluator(std::vector<std::vector<double> > vWeights);

	double O_Fitna(Protein &target, ProteinProfile& profiles, DFIRE2& dDFIRE_Inst);
	double PotScore(Protein& target,
	                std::array<double, 20>& aPot_Bar,
	                std::array<double, 20>& aPot_Stdev,
	                int dS_Parameter);

	double Solvent_Score(Protein& target,
	                     std::array<std::array<double, 7>, 20>& aProfile);

	double Secondary_Struct(Protein& target,
	                        std::array<std::array<double, 7>, 20>& aProfile);

	double AlignmentScore(Protein& target, std::vector<std::array<double, 7> >& aProfile);

	std::array<double, 20> FrequencyScore(Protein& target,
	                                      std::array<double, 20>& aAA_Freq_Mean,
	                                      std::array<double, 20>& aAA_Freq_Stdev);

	double eNormalaize(double e, double lB, double uB);

	double operator()(Protein &target, ProteinProfile& profiles, DFIRE2& dDFIRE_Inst);
private:
	std::vector<std::vector<double> > _vecWeights;
};


#endif // __EVALUATOR_H__
