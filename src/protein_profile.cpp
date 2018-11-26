#include "protein_profile.h"
using namespace std::placeholders;


ProteinProfile::ProteinProfile(Protein target)
    : _refProtein ( target )
{
    _aAlgn_Profile.resize(_refProtein.length());
    _aAlgn_Position_Count.resize(_refProtein.length());
    _matFragments.resize(_refProtein.length(),
                         std::vector<std::vector<std::string> >( _refProtein.length() ));
}



void ProteinProfile::CalculateProfiles(int nFlag,
                                       bool bVerbose,
                                       bool bSave_Frags,
                                       double dFrag_Score_Cutoff,
                                       double dGap_Penalty,
                                       double dPotS_Param,
                                       int nMin_Frag)
{
    if (_vecHomologous_Proteins.empty()) {
        std::cerr << "FatalError: No protein to process.";
        return;
    }

    std::vector<std::function<void(Protein&)> > vecProcFuncs;
    std::vector<std::function<void()> > vecCalcFuncs;
    std::string Msg = "Calculating ";


    if (nFlag & 1) {
        vecProcFuncs.push_back(std::bind(&ProteinProfile::Process_Solvent, this, _1));
        vecCalcFuncs.push_back(std::bind(&ProteinProfile::Calculate_Solvent_Profile, this));
        Msg += "Solvent, ";
    }
    if (nFlag & 2) {
        vecProcFuncs.push_back(std::bind(&ProteinProfile::Process_SS, this, _1));
        vecCalcFuncs.push_back(std::bind(&ProteinProfile::Calculate_SS_Profile, this));
        Msg += "SS, ";
    }
    if (nFlag & 4) {
        vecProcFuncs.push_back(std::bind(&ProteinProfile::Process_Pot_AAFreq, this, _1));
        vecCalcFuncs.push_back(std::bind(&ProteinProfile::Calculate_Pot_AAFreq_Profile, this));
        _dPotS_Param = dPotS_Param;
        Msg += "Pot, ";
    }
    if (nFlag & 8) {
        vecCalcFuncs.push_back([ = ] { Calculate_Alignment_Profile(bSave_Frags, dFrag_Score_Cutoff,
                                       dGap_Penalty, nMin_Frag);
                                     });
        Msg += "Algn, ";
    }

    Msg.erase(Msg.size() - 2);
    Msg += " profiles";

    Thread_Manager(vecProcFuncs, _vecHomologous_Proteins, bVerbose, Msg);
    for (auto CalcFunction : vecCalcFuncs)
        CalcFunction();
}



void ProteinProfile::CalculateRemainingProfiles(bool bVerbose) {
    int nFlag = 0;
    bool abFlags[] = {bSolvent_Rdy, bSS_Rdy, bPot_Rdy, bAlgn_Rdy};
    for (int i = 0, j = 1; i < 4; ++i, j *= 2)
        if (abFlags[i])
            nFlag += j;
    CalculateProfiles(nFlag, bVerbose);
}



void ProteinProfile::Find_Homologous_Proteins(std::vector<std::string> vecDB,
        double dAlgn_Score_CutOff, double bVerbose)
{
    _vecHomologous_Proteins.clear();
    _vecTupleAlignments.clear();

    _dAlgn_Score_CutOff = dAlgn_Score_CutOff;

    std::vector<std::function<void(std::string&)> > vecProcs;
    vecProcs.push_back(std::bind(&ProteinProfile::Process_IsHomologue, this, _1));
    Thread_Manager(vecProcs, vecDB, bVerbose, "Searching for proteins homologous to the target structure");
}



template<class FuncType, class Args>
void ProteinProfile::Thread_Manager(std::vector<std::function<FuncType> > vecFuncs,
                                    std::vector<Args>& vecDB, bool bVerbose,
                                    std::string sMsg, int nThreads)
{
    std::vector<std::thread> vecThreads;
    int size = vecDB.size(),
        nCount_Now = 0;

    nThreads = (nThreads > 0 ? nThreads : FALLBACK_THRD_NUM);

    for (int i = 0; i < nThreads; ++i) {
        std::vector vecBatch(vecDB.begin() + i * size / nThreads,
                             vecDB.begin() + (i + 1) * (size) / nThreads);
        vecThreads.push_back(std::thread([ =, &nCount_Now] { Processing_Thread(vecFuncs, vecBatch, nCount_Now); })); // BUG #1: must send by copy
    }

    if (bVerbose) {
        utility::Progress_Indicator(sMsg, 0, size);
        while (true) {
            usleep(500 * 1000);
            utility::Progress_Indicator(sMsg, nCount_Now, size);
            if (nCount_Now == size)
                break;
        }
    }

    for (int i = 0; i < nThreads; ++i)
        vecThreads[i].join();
}



template<class FuncType, class Args>
void ProteinProfile::Processing_Thread(std::vector<std::function<FuncType> > vecFuncs,
                                       std::vector<Args> vecBatch, int& nCount_Now)
{
    for (auto& data : vecBatch) {
        for (auto ProcFunction : vecFuncs)
            ProcFunction(data);

        std::lock_guard<std::mutex> lock(_mtx_count);
        nCount_Now++;
    }
}



void ProteinProfile::Process_IsHomologue(std::string fPath) {
    if (_refProtein.md5 == utility::File_md5(fPath))
        return;

    std::string stdout, line, algn_line1, algn_line2;
    auto oBuffer = subprocess::check_output({ex_TMALIGN.c_str(), _refProtein.fPath.c_str(),
                                            fPath.c_str(), "-a"
                                            });
    std::stringstream ret(oBuffer.buf.data());

    double score;
    getline(ret, line); getline(ret, line);

    ret >> line >> score;
    if (score < _dAlgn_Score_CutOff)
        return;

    ret >> algn_line1 >> algn_line2;
    assert(algn_line1.size() == algn_line2.size());
    int align_len = algn_line1.size();

    std::vector<double> atom_dist( align_len + 1 );
    for (int i = 1; i <= align_len; ++i)
        ret >> atom_dist[i];

    std::lock_guard<std::mutex> lock(_mtx_algn);
    _vecTupleAlignments.emplace_back(score, algn_line1, algn_line2, atom_dist);
    _vecHomologous_Proteins.emplace_back(Protein(fPath));
}



void ProteinProfile::Process_Solvent(Protein& target) {
    target.Calculate_Solvent();

    std::lock_guard<std::mutex> lock(_mtx_solvent);
    for (int i = 0, n = target.length(); i < n; ++i) {
        if (!target.IsStandardAA(target[i]))
            continue;

        int idx = sym_to_idx[ target[i] ];
        _aSolvent_AA_Count[idx][ target.aSolvent_Accessibility[i] ]++;
        _aAA_Total_Count[idx]++;
    }
}



void ProteinProfile::Process_SS(Protein& target) {
    target.Calculate_SS();

    std::lock_guard<std::mutex> lock(_mtx_SS);
    for (int i = 0, n = target.length(); i < n; ++i)
        if (target.IsStandardAA(target[i]))
            _aSec_AA_Count[ sym_to_idx[target[i]] ][ target.aSecondary_Structure[i] ]++;
}



void ProteinProfile::Process_Pot_AAFreq(Protein& target) {
    target.Calculate_Pot(_dPotS_Param);

    std::lock_guard<std::mutex> lock(_mtx_pot);
    for (int i = 0; i < 20; ++i) {
        _vecPotScores[i].push_back(target.aPot_Values[i]);
        _vecAA_Freqs[i].push_back(target.aAA_Freqs[i]);
    }
}



void ProteinProfile::Calculate_Alignment_Profile(bool bSave_Frags, double dFrag_Score_Cutoff,
        double dGap_Penalty, int nMin_Frag)
{
    _dFrag_Score_Cutoff = dFrag_Score_Cutoff;
    _dGap_Penalty = dGap_Penalty;
    _nMin_Frag = nMin_Frag;

    for (auto& [score, seq1, seq2, atom_dist] : _vecTupleAlignments) {
        int align_len = seq1.length();
        std::vector<double> sum_dist( align_len + 1 );

        for (int i = 1; i <= align_len; ++i)
            sum_dist[i] = (atom_dist[i] == GAP_NUMBER ? _dGap_Penalty : atom_dist[i])
                          + sum_dist[i - 1];

        for (int i = 0, pos = 0; i < align_len; ++i) {
            if (seq1[i] == '-')
                continue;
            pos++;
            if (seq2[i] == '-' || seq2[i] == 'X')
                continue;

            assert(algn_classes.count(seq2[i]) != 0);
            int nClass = algn_classes[seq2[i]];

            _aAlgn_Profile[pos - 1][nClass] += 10 - std::min(10., atom_dist[i + 1]); // Convert to similarity
            _aAlgn_Position_Count[pos - 1]++;
            _aAlgn_Class_Count[nClass]++;
        }

        if (bSave_Frags) {
            int pos = 0;
            for (int i = 0; i < align_len; ++i) {
                for (int l = nMin_Frag; i + l <= align_len; ++l) {
                    if (atom_dist[i + 1] == 0 || atom_dist[i + l] == 0)
                        continue;
                    double sc = sum_dist[i + l] - sum_dist[i];
                    assert(sc > 0);
                    if (sc <= l * dFrag_Score_Cutoff) // hit
                    {
                        std::string fragment;
                        int nFragLen = 0;
                        for (int j = i; j < i + l; ++j) {
                            if (seq1[j] == '-') // skip gaps on seq1
                                continue;
                            fragment += seq2[j];
                            nFragLen++;
                        }
                        _matFragments[pos][pos + nFragLen - 1].push_back(fragment);
                    }
                }
                if (seq1[i] != '-')
                    pos++;
            }
        }
    }

    long long nAA_Total = std::accumulate(_aAlgn_Class_Count.begin(),
                                          _aAlgn_Class_Count.end(), 0);

    for (int i = 0; i < _refProtein.length(); ++i) {
        double posScore_Sum = std::accumulate(_aAlgn_Profile[i].begin(),
                                              _aAlgn_Profile[i].end(), 0);

        for (int j = 0; j < 7; ++j)
            if (posScore_Sum != 0) {
                _aAlgn_Profile[i][j] = ( (_aAlgn_Position_Count[i] * _aAlgn_Profile[i][j] / posScore_Sum) +
                                         ((double)_aAlgn_Class_Count[j] / nAA_Total) * std::sqrt(_aAlgn_Position_Count[i])
                                       )
                                       / (_aAlgn_Position_Count[i] + std::sqrt(_aAlgn_Position_Count[i]));
            }
    }
    bAlgn_Rdy = true;
    bFrags_Rdy |= bSave_Frags;
}



void ProteinProfile::Calculate_Solvent_Profile() {
    std::array<long long, 7> solvent_class_total{};

    long long nAA_Total = std::accumulate(_aAA_Total_Count.begin(), _aAA_Total_Count.end(), 0);

    for (int i = 0; i < 7; ++i)
        for (int j = 0; j < 20; ++j)
            solvent_class_total[i] += _aSolvent_AA_Count[j][i];

    for (int i = 0; i < 20; ++i)
        for (int j = 0; j < 7; ++j)
            if (solvent_class_total[j] && _aAA_Total_Count[i])
                _aSolvent_Profile[i][j] = (double)_aSolvent_AA_Count[i][j] / (solvent_class_total[j]
                                          * (_aAA_Total_Count[i] / (double)nAA_Total));
    bSolvent_Rdy = true;
}



void ProteinProfile::Calculate_SS_Profile() {
    std::array<long long, 7> sec_class_total{};

    long long nAA_Total = std::accumulate(_aAA_Total_Count.begin(), _aAA_Total_Count.end(), 0);

    for (int i = 0; i < 7; ++i)
        for (int j = 0; j < 20; ++j)
            sec_class_total[i] += _aSec_AA_Count[j][i];

    for (int i = 0; i < 20; ++i)
        for (int j = 0; j < 7; ++j)
            if (sec_class_total[j] && _aAA_Total_Count[i])
                _aSec_Profile[i][j] = (double)_aSec_AA_Count[i][j] / (sec_class_total[j]
                                      * (_aAA_Total_Count[i] / (double)nAA_Total));

    bSS_Rdy = true;
}



void ProteinProfile::Calculate_Pot_AAFreq_Profile() {
    for (int i = 0; i < 20; ++i) {
        if (_vecPotScores[i].size() >= 2) {
            _vecPotScores[i] = OutlierElimination_IQR(_vecPotScores[i]);
            std::tie(_aPot_Bar[i], _aPot_Stdev[i]) = Mean_SD(_vecPotScores[i]);
        }

        if (_vecAA_Freqs[i].size() >= 2) {
            _vecAA_Freqs[i] = OutlierElimination_IQR(_vecAA_Freqs[i]);
            std::tie(_aAA_Freq_Mean[i], _aAA_Freq_Stdev[i]) = Mean_SD(_vecAA_Freqs[i]);
        }
    }
    bPot_Rdy = true;
}



void ProteinProfile::Write_ToFile(bool bWriteCounts, std::string sDirectory) {
    assert(utility::DirectoryExists(sDirectory));

    if (bSolvent_Rdy) {
        std::ofstream outFile(sDirectory + RelativeFileName("solvent"));

        outFile << QuickInfo() << "\n\n" << "#SOLVENT_PROFILE\n\n"
                << sFileStartSym << "\n";

        for (int i = 0; i < 20; ++i) {
            outFile << idx_to_sym[i];
            for (int j = 0; j < 7; ++j)
                outFile << " " << std::fixed << std::setprecision(10) << _aSolvent_Profile[i][j];
            outFile << "\n";
        }
    }

    if (bSS_Rdy) {
        std::ofstream outFile(sDirectory + RelativeFileName("sec"));

        outFile << QuickInfo() << "\n\n" << "#SECONDARY_STRUCTURE_PROFILE\n\n"
                << sFileStartSym << "\n";

        for (int i = 0; i < 20; ++i) {
            outFile << idx_to_sym[i];
            for (int j = 0; j < 7; ++j)
                outFile << " " << std::fixed << std::setprecision(10) << _aSec_Profile[i][j];
            outFile << "\n";
        }
    }

    if (bPot_Rdy) {
        std::ofstream outFile(sDirectory + RelativeFileName("pot"));

        outFile << QuickInfo() << "\n\n" << "#POT_STATISTIC\n\n" << sFileStartSym << "\n";

        for (int i = 0; i < 20; ++i)
            outFile << idx_to_sym[i] << " " << std::fixed << std::setprecision(10)
                    << _aPot_Bar[i] << " " << _aPot_Stdev[i] << "\n";


        outFile << "\n\n\n#AA_FREQUENCY\n\n" << sFileStartSym << "\n";
        for (int i = 0; i < 20; ++i)
            outFile << idx_to_sym[i] << " " << std::fixed << std::setprecision(10)
                    << _aAA_Freq_Mean[i] << " " << _aAA_Freq_Stdev[i] << "\n";

    }

    if (bAlgn_Rdy) {
        std::ofstream outFile(sDirectory + RelativeFileName("alignment"));
        outFile << QuickInfo(true) << "\n\n" << "#ALIGNMENT_PROFILE\n\n" << sFileStartSym << "\n";

        for (int i = 0; i < _refProtein.length(); ++i) {
            outFile << std::left << std::setw(3) << i + 1;
            for (int j = 0; j < 7; ++j)
                outFile << " " << std::fixed << std::setprecision(10) << _aAlgn_Profile[i][j];
            outFile << "\n";
        }

        if (bFrags_Rdy) {
            std::ofstream outFile(sDirectory + RelativeFileName("fragment"));
            outFile << QuickInfo(true) << "\n\n" << "#ALIGNMENT_FRAGMENTS\n\n" << sFileStartSym << "\n";
            for (int i = 0; i < _refProtein.length(); ++i)
                for (int j = i; j < _refProtein.length(); ++j)
                    for (auto& sFragment : _matFragments[i][j])
                        outFile << i << " " << j << " " << sFragment << "\n";
        }

    }

    if (bWriteCounts) {
        std::ofstream outFile(sDirectory + RelativeFileName("count"));
        outFile << "#BURIAL_COUNT\n\n";
        for (int i = 0; i < 20; ++i) {
            outFile << idx_to_sym[i];
            for (int j = 0; j < 7; ++j)
                outFile << " " << std::setw(7) << _aSolvent_AA_Count[i][j];
            outFile << "\n";
        }

        outFile << "\n\n#STRUCTURE_COUNT\n\n";
        for (int i = 0; i < 20; ++i) {
            outFile << idx_to_sym[i];
            for (int j = 0; j < 7; ++j)
                outFile << " " << std::setw(7) << _aSec_AA_Count[i][j];
            outFile << "\n";
        }

        outFile << "\n\n#AA_COUNT\n\n";
        for (int i = 0; i < 20; ++i) {
            outFile << idx_to_sym[i]
                    << " " << _aAA_Total_Count[i] << "\n";
        }
    }
}



void ProteinProfile::Read_FromFile(std::string sDirectory) {
    assert(utility::DirectoryExists(sDirectory));
    std::string HEAD;
    double PARAM;

    std::string fSolvName = sDirectory + RelativeFileName("solvent");
    if (utility::FileExists(fSolvName)) {
        std::ifstream inFile(fSolvName);

        while (HEAD != sFileStartSym)
            inFile >> HEAD;

        for (int i = 0; i < 20; ++i) {
            inFile >> HEAD;
            for (int j = 0; j < 7; ++j)
                inFile >> _aSolvent_Profile[i][j];
        }

        bSolvent_Rdy = true;
    }

    std::string fSecName = sDirectory + RelativeFileName("sec");
    if (utility::FileExists(fSolvName)) {
        std::ifstream inFile(fSecName);

        HEAD = "";
        while (HEAD != sFileStartSym)
            inFile >> HEAD;

        for (int i = 0; i < 20; ++i) {
            inFile >> HEAD;
            for (int j = 0; j < 7; ++j)
                inFile >> _aSec_Profile[i][j];
        }

        bSS_Rdy = true;
    }

    std::string fPotName = sDirectory + RelativeFileName("pot");
    if (utility::FileExists(fPotName)) {
        std::ifstream inFile(fPotName);

        HEAD = "";
        while (HEAD != sFileStartSym) {
            if (HEAD.find("POT_S_PARAMETER") == 0) {
                inFile >> PARAM;
                _dPotS_Param = PARAM;
            }
            inFile >> HEAD;
        }

        for (int i = 0; i < 20; ++i)
            inFile >> HEAD >> _aPot_Bar[i] >> _aPot_Stdev[i];


        while (HEAD != sFileStartSym)
            inFile >> HEAD;
        for (int i = 0; i < 20; ++i)
            inFile >> HEAD >> _aAA_Freq_Mean[i] >> _aAA_Freq_Stdev[i];

        bPot_Rdy = true;
    }

    std::string fAlgnName = sDirectory + RelativeFileName("alignment");
    if (utility::FileExists(fAlgnName)) {
        std::ifstream inFile(fAlgnName);

        while (HEAD != sFileStartSym) {
            if (HEAD.find("ALGN_SCORE_CUTOFF") == 0) {
                inFile >> PARAM;
                _dAlgn_Score_CutOff = PARAM;
            }
            if (HEAD.find("FRAG_DIST_CUTOFF") == 0) {
                inFile >> PARAM;
                _dFrag_Score_Cutoff = PARAM;
            }
            if (HEAD.find("GAP_PENALTY") == 0) {
                inFile >> PARAM;
                _dGap_Penalty = PARAM;
            }
            if (HEAD.find("MIN_FRAG_LENGTH") == 0) {
                inFile >> PARAM;
                _nMin_Frag = PARAM;
            }
            inFile >> HEAD;
        }

        for (int i = 0; i < _refProtein.length(); ++i) {
            inFile >> HEAD;
            for (int j = 0; j < 7; ++j)
                inFile >> _aAlgn_Profile[i][j];
        }

        bAlgn_Rdy = true;

        std::string fFragName = sDirectory + RelativeFileName("fragment");
        if (utility::FileExists(fFragName)) {
            std::ifstream inFile(fFragName);

            while (HEAD != sFileStartSym)
                inFile >> HEAD;

            for (int i, j; inFile >> i >> j >> HEAD; )
                _matFragments[i][j].emplace_back(HEAD);

            bFrags_Rdy = true;
        }
    }
}



std::string ProteinProfile::QuickInfo(bool bIncludeAlignInfo) {
    std::ostringstream output;
    output << "FILE_PATH:  " << _refProtein.fPath << "\n"
           << "MD5:  " << _refProtein.md5 << "\n"
           << "PROT_LEN:  " << _refProtein.length() << "\n"
           << "TEMPLATES_COUNT:  " << _vecHomologous_Proteins.size() << "\n"
           << "ALGN_SCORE_CUTOFF:  " << _dAlgn_Score_CutOff << "\n"
           << "POT_S_PARAMETER:  " << _dPotS_Param << "\n";
    if (bIncludeAlignInfo) {
        output << "FRAG_DIST_CUTOFF:  " << _dFrag_Score_Cutoff << "\n"
               << "GAP_PENALTY:  " << _dGap_Penalty << "\n"
               << "MIN_FRAG_LENGTH:  " << _nMin_Frag << "\n";
    }
    return output.str();
}



std::string ProteinProfile::RelativeFileName(std::string sPN) {
    std::string md5_First15 = _refProtein.md5.substr(0, 15);

    if (sPN == "Solvent" || sPN == "solvent" || sPN == "burial")
        return "Solvent_Prof_" + md5_First15;

    if (sPN == "Pot" || sPN == "pot" || sPN == "AAFreq")
        return "Pot_AA_Prof_" + md5_First15;

    if (sPN == "SS" || sPN == "Sec" || sPN == "sec" || sPN == "secondary")
        return "Sec_Prof_" + md5_First15;

    if (sPN == "algn" || sPN == "Algn" || sPN == "alignment")
        return "Algn_Prof_" + md5_First15;

    if (sPN == "frags" || sPN == "fragment" || sPN == "fragments")
        return "Algn_Frags_" + md5_First15;

    if (sPN == "Counts" || sPN == "count" || sPN == "counts")
        return "Counts_" + md5_First15;

    throw std::runtime_error("Unknown profile name.");
}



bool ProteinProfile::FragmentExists(int i, int j) {
    return !_matFragments[i][j].empty();
}



std::string ProteinProfile::FragmentFetch(int i, int j) {
    return _matFragments[i][j][rand() % _matFragments[i][j].size()];
}