#include "evaluator.h"


Evaluator::Evaluator(std::string fPath) {
    assert(FileExists(fPath));

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
    score[4] = PotScore(target, profiles._aPot_Bar, profiles._aPot_Stdev, profiles._potS_Param);
    auto ret = FrequencyScore(target, profiles._aAA_Freq_Mean, profiles._aAA_Freq_Stdev);
    std::copy(ret.begin(), ret.end(), score.begin() + 5);

    for (int i = 0; i < _vecWeights.size(); ++i)
        dFinal_Score += _vecWeights[i][2] * eNormalaize(score[i], _vecWeights[i][0], _vecWeights[i][1]);

    return dFinal_Score;
}


double Evaluator::PotScore(Protein& target,
                           std::array<double, 20>& aPot_Bar,
                           std::array<double, 20>& aPot_Stdev,
                           int dS_Parameter)
{
    std::array<double, 20> AA_freq{}, pot_r, pot0, pot1{}, var0;
    int n = target.length();

    for (int i = 0; i < n; ++i) {
        if (sym_to_idx.count(target[i]) == 0)
            continue;
        AA_freq[ sym_to_idx[target[i]] ]++;
    }

    for (int i = 0; i < 20; ++i) {
        if (AA_freq[i] < 2)
            continue;
        pot_r[i] = std::exp(-std::sqrt(AA_freq[i]) / (n * dS_Parameter));
        pot0[i] = (AA_freq[i] * (AA_freq[i] - 1) * pot_r[i] * (n - 1 / (1 - pot_r[i])))
                  / (n * (n - 1) * (1 - pot_r[i]));
        var0[i] = std::sqrt(std::pow(pot_r[i] * AA_freq[i] * (1 - AA_freq[i] / n), 2)
                            / ((1 - pot_r[i] * pot_r[i]) * n));
    }

    for (int i = 0; i < n; ++i) {
        if (sym_to_idx.count(target[i]) == 0)
            continue;
        int idx = sym_to_idx[target[i]];
        for (int j = i + 1; j < n; ++j)
            pot1[idx] += (target[i] == target[j]) * std::pow(pot_r[idx], j - i);
    }

    double avg = 0;
    int excluded = 0;
    for (int i = 0; i < 20; ++i) {
        if (AA_freq[i] < 2) {
            excluded += AA_freq[i];
            continue;
        }
        double potScore_i = (pot1[i] - pot0[i]) / var0[i];
        avg += AA_freq[i] * ( 0.5 * std::pow((potScore_i - aPot_Bar[i]) / aPot_Stdev[i], 2)
                              - log(1 / (aPot_Stdev[i] * std::sqrt(2 * M_PI))) );
    }
    return avg / (n - excluded);
}


double Evaluator::Solvent_Score(Protein& target,
                                std::array<std::array<double, 7>, 20>& aProfile) {
    double score = 0;
    int n = target.length(),
        excluded = 0;
    for (int i = 0; i < n; ++i) {
        if (sym_to_idx.count(target[i]) == 0) {
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
    double score = 0;
    int n = target.length(),
        excluded = 0;
    for (int i = 0; i < n; ++i) {
        if (sym_to_idx.count(target[i]) == 0) {
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
        if (sym_to_idx.count(target[i]) == 0) {
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
    std::array<double, 20> aAA_freq{}, aAA_score{};
    auto& target_sequence = target.sequence;

    for (auto c : target_sequence)
        if (sym_to_idx.count(c) != 0)
            aAA_freq[sym_to_idx[c]]++;

    for (int i = 0; i < 20; ++i) {
        aAA_freq[i] /= target_sequence.length();
        aAA_score[i] = std::abs((aAA_freq[i] - aAA_Freq_Mean[i]) / aAA_Freq_Stdev[i]);
    }
    return aAA_score;
}


double Evaluator::eNormalaize(double e, double lB, double uB) {
    e = (e - lB) / (uB - lB);
    return std::min(1., std::max(0., e));
}


double Evaluator::operator()(Protein &target, ProteinProfile& profiles, DFIRE2& dDFIRE_Inst) {
    return O_Fitna(target, profiles, dDFIRE_Inst);
}