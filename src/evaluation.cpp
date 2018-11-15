#include "evaluation.h"



double PotScore(std::string& target_sequence) {
    std::array<double, 20> pot_frequency{}, pot_r, pot0, pot1{}, var0;
    int n = target_sequence.size();

    for (int i = 0; i < n; ++i) {
        if (symbol_map.count(target_sequence[i]) == 0)
            continue;
        pot_frequency[symbol_map[target_sequence[i]]]++;
    }

    for (int i = 0; i < 20; ++i) {
        if (pot_frequency[i] < 2)
            continue;
        pot_r[i] = std::exp(-std::sqrt(pot_frequency[i]) / (n * Pot_S_Constant));
        pot0[i] = (pot_frequency[i] * (pot_frequency[i] - 1) * pot_r[i] * (n - 1 / (1 - pot_r[i])))
                  / (n * (n - 1) * (1 - pot_r[i]));
        var0[i] = std::sqrt(std::pow(pot_r[i] * pot_frequency[i] * (1 - pot_frequency[i] / n), 2)
                            / ((1 - pot_r[i] * pot_r[i]) * n));
    }
    for (int i = 0; i < n; ++i) {
        if (symbol_map.count(target_sequence[i]) == 0)
            continue;
        int idx = symbol_map[target_sequence[i]];
        for (int j = i + 1; j < n; ++j) {
            pot1[idx] += (target_sequence[i] == target_sequence[j])
                         * std::pow(pot_r[idx], j - i);
        }
    }

    double avg = 0;
    int excluded = 0;
    for (int i = 0; i < 20; ++i) {
        if (pot_frequency[i] < 2) {
            excluded += pot_frequency[i];
            continue;
        }

        double pot_score_i = (pot1[i] - pot0[i]) / var0[i];
        //     std::cout << "Pot " << index_to_symbol[i] << ":\t" << "raw: " << pot_score_i << "  calc:" << (0.5 * std::pow((pot_score_i - pot_bar[i]) / pot_stdev[i], 2)
        //               - log(1 / (pot_stdev[i] * std::sqrt(2 * M_PI)))) << "\n";
        avg += pot_frequency[i] * ( 0.5 * std::pow((pot_score_i - pot_bar[i]) / pot_stdev[i], 2)
                                    - log(1 / (pot_stdev[i] * std::sqrt(2 * M_PI))) );
    }
    return avg / (n - excluded);
}


double BurialScore(std::string& target_sequence) {
    double score = 0;
    int n = target_sequence.size(),
        excluded = 0;
    for (int i = 0; i < n; ++i) {
        if (symbol_map.count(target_sequence[i]) == 0) {
            excluded++;
            continue;
        }
        int idx = symbol_map[target_sequence[i]];
        score += burial_potential[idx][target_solvent_accessibility[i]];
    }
    return score / (n - excluded);
}


double SecondaryStructScore(std::string& target_sequence) {
    double score = 0;
    int n = target_sequence.size(),
        excluded = 0;
    for (int i = 0; i < n; ++i) {
        if (symbol_map.count(target_sequence[i]) == 0) {
            excluded++;
            continue;
        }
        int idx = symbol_map[target_sequence[i]];
        score += structure_potential[idx][target_secondary_structure[i]];
    }
    return score / (n - excluded);
}


double AlignmentScore(std::string& target_sequence) {
    double score = 0;
    int n = target_sequence.size(),
        excluded = 0;
    for (int i = 0; i < n; ++i) {
        if (symbol_map.count(target_sequence[i]) == 0) {
            excluded++;
            continue;
        }
        char sym = target_sequence[i];
        score += alignment_potential[i][AA_classification[sym]];
    }
    return score / (n - excluded);
}


double FrequencyScore(std::string& target_sequence) {
    std::array<double, 20> AA_freq{};
    double total_score = 0;

    for (auto c : target_sequence) {
        if (symbol_map.count(c) > 0)
            AA_freq[symbol_map[c]]++;
    }
    for (int i = 0; i < 20; ++i) {
        AA_freq[i] /= target_sequence.length();
        double sc = std::abs((AA_freq[i] - AA_freq_mean[i]) / AA_freq_sd[i]);
        // std::cerr << "FREQ " << index_to_symbol[i] << "  " << (sc) << " " << eNormalaize(sc, freq_bounds[i][0], freq_bounds[i][1]) << "\n";
        total_score += AA_freq[i] * eNormalaize(sc, freq_bounds[i][0], freq_bounds[i][1]);
    }
    return total_score;
}


double O_Fitna(std::string target_sequence, double ret[]) {
    double score[WGT_NUM], final_score = 0;
    int n = target_sequence.size();

    score[0] = BurialScore(target_sequence);
    score[1] = SecondaryStructScore(target_sequence);
    score[2] = AlignmentScore(target_sequence);
    score[3] = dDFIRE_CFE(target_sequence, target_atom_distance) / (dDFIRE_COEF * n);
    score[4] = PotScore(target_sequence);
    score[5] = FrequencyScore(target_sequence);
    // std::cerr << "POT: :\t" << score[4] << "\n";

    for (int i = 0; i < WGT_NUM; ++i) {
        double d = eNormalaize(score[i], sc_weight[i][0], sc_weight[i][1]) *
                   sc_weight[i][2];
        // std::cerr << eNormalaize(score[i], sc_weight[i][0], sc_weight[i][1]) << " "  << d << "  " << i << "\n";
        final_score += d;
        if (ret != NULL)
            ret[i] = d;
    }
    return final_score;
}


double eNormalaize(double e, double lB, double uB) {
    e = (e - lB) / (uB - lB);
    return std::min(1., std::max(0., e));
}
