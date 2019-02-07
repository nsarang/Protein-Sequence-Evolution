#include "evaluator.h"


Evaluator::Evaluator(std::string fPath) {
    assert(utility::FileExists(fPath));

    std::ifstream inFile(fPath);
    int n; std::string HEAD;
    inFile >> HEAD >> n;
    assert(HEAD == "#FITNESS_WEIGHTS" && n > 0);

    _vecWeights.resize(n, std::vector<double>( 3 ));
    for (int i = 0; i < n; ++i) {
        inFile >> HEAD;
        for (int j = 0; j < 3; ++j)
            inFile >> _vecWeights[i][j];
    }
}


Evaluator::Evaluator(std::vector<std::vector<double> > vWeights) :
    _vecWeights( vWeights )
{
}


double Evaluator::O_Fitna(Protein &target, ProteinProfile& profiles, DFIRE2& dDFIRE_Inst) {
    std::vector<double> score( _vecWeights.size() );
    double dFinal_Score = 0;

    score[0] = Solvent_Score(target, profiles._aSolvent_Profile);
    score[1] = Secondary_Struct(target, profiles._aSec_Profile);
    score[2] = AlignmentScore(target, profiles._aAlgn_Profile);
    score[3] = dDFIRE_Inst.Calc_CFE(target) / (target.length() * dDFIRE_COEF);
 //   score[4] = PotScore(target, profiles._aPot_Bar, profiles._aPot_Stdev, profiles._dPotS_Param);
    auto ret = FrequencyScore(target, profiles._aAA_Freq_Mean, profiles._aAA_Freq_Stdev);
    std::copy(ret.begin(), ret.end(), score.begin() + 5);

    for (int i = 0, n = _vecWeights.size(); i < n; ++i) // Avoid gcc comparison warning
        dFinal_Score += _vecWeights[i][2] * eNormalaize(score[i], _vecWeights[i][0], _vecWeights[i][1]);

    return dFinal_Score;
}


std::array<double, 20> Evaluator::PotScore(Protein& target,
        std::array<double, 20>& aPot_Bar,
        std::array<double, 20>& aPot_Stdev,
        double dPotS_Param)
{
    target.Calculate_Pot(dPotS_Param);
    std::array<double, 20> retPot_scores{};


    int n = target.length();
    for (int i = 0; i < 20; ++i) {
        if (target.aAA_Freqs[i] * n < 2)
            continue;
        retPot_scores[i] = 0.5 * std::pow((target.aPot_Values[i] - aPot_Bar[i]) / aPot_Stdev[i], 2)
                           - log(1 / (aPot_Stdev[i] * std::sqrt(2 * M_PI)));
    }
    return retPot_scores;
}


double Evaluator::Solvent_Score(Protein& target,
                                std::array<std::array<double, 7>, 20>& aProfile) {
    target.Calculate_Solvent();

    double score = 0;
    int n = target.length(),
        excluded = 0;
    for (int i = 0; i < n; ++i) {
        if (!target.IsStandardAA(target[i])) {
            excluded++;
            continue;
        }
        int idx = sym_to_idx[target[i]];
        score += aProfile[idx][ target.aSolvent_Accessibility[i] ];
    }
    return score / (n - excluded);
}


double Evaluator::Secondary_Struct(Protein& target,
                                   std::array<std::array<double, 7>, 20>& aProfile) {
    target.Calculate_SS();

    double score = 0;
    int n = target.length(),
        excluded = 0;
    for (int i = 0; i < n; ++i) {
        if (!target.IsStandardAA(target[i])) {
            excluded++;
            continue;
        }
        int idx = sym_to_idx[target[i]];
        score += aProfile[idx][ target.aSecondary_Structure[i] ];
    }
    return score / (n - excluded);
}


double Evaluator::AlignmentScore(Protein& target,
                                 std::vector<std::array<double, 7> >& aProfile) {
    double score = 0;
    int n = target.length(),
        excluded = 0;
    for (int i = 0; i < n; ++i) {
        if (!target.IsStandardAA(target[i])) {
            excluded++;
            continue;
        }
        char sym = target[i];
        score += aProfile[i][ algn_classes[sym] ];
    }
    return score / (n - excluded);
}


std::array<double, 20> Evaluator::FrequencyScore(Protein& target,
        std::array<double, 20>& aAA_Freq_Mean,
        std::array<double, 20>& aAA_Freq_Stdev)
{
    std::array<double, 20> retAA_scores{};

    for (int i = 0; i < 20; ++i)
        retAA_scores[i] = std::abs((target.aAA_Freqs[i] - aAA_Freq_Mean[i]) / aAA_Freq_Stdev[i]);

    return retAA_scores;
}


double Evaluator::eNormalaize(double e, double lB, double uB) {
    e = (e - lB) / (uB - lB);
    return std::min(1., std::max(0., e));
}


double Evaluator::operator()(Protein &target, ProteinProfile& profiles, DFIRE2& dDFIRE_Inst) {
    return O_Fitna(target, profiles, dDFIRE_Inst);
}