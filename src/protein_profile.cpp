#include "protein_profile.h"
using namespace std::placeholders;


ProteinProfile::ProteinProfile(Protein target)
    : _refProtein ( target )
{
    _aAlgn_Profile.resize(_refProtein.length());
    _aAlgn_Position_Count.resize(_refProtein.length());
    _matFragments.resize(_refProtein.length(), std::vector<std::vector<std::string> >( _refProtein.length() ));
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 7; ++j)
        {
            assert(_aSolvent_Profile[i][j] == 0);
            assert(_aSec_Profile[i][j] == 0);
            assert(_aSolvent_AA_Count[i][j] == 0);
            assert(_aSec_AA_Count[i][j] == 0);
        }
    }
}


void ProteinProfile::CalculateProfiles(bool bAlgn, bool bSolvent, bool bPot, bool bSS,
                                       bool bVerbose, bool bSave_Frags, int potS_Param,
                                       double dFrag_Score_Cutoff, double dGap_Score, int nMin_Frag)
{
    if (_vecHomologous_Proteins.empty()) {
        std::cerr << "FatalError: No protein to process.";
        return;
    }

    bAlgn_Rdy |= bAlgn;
    bSolvent_Rdy |= bSolvent;
    bPot_Rdy |= bPot;
    bSS_Rdy |= bSS;
    bFrags_Rdy |= bSave_Frags;

    std::vector<std::function<void(Protein&)> > vecProcFuncs;
    std::vector<std::function<void()> > vecCalcFuncs;
    std::string Msg = "Calculating ";
    std::cerr << Msg << std::endl;

    if (bAlgn) {
        vecCalcFuncs.push_back([ = ] { Calculate_Alignment_Profile(bSave_Frags, dFrag_Score_Cutoff,
                                       dGap_Score, nMin_Frag);
                                     });
        Msg += "Algn, ";
    }
    if (bSolvent) {
        vecProcFuncs.push_back(std::bind(&ProteinProfile::Process_Solvent, this, _1));
        vecCalcFuncs.push_back(std::bind(&ProteinProfile::Calculate_Solvent_Profile, this));
        Msg += "Solvent, ";
    }
    if (bPot) {
        vecProcFuncs.push_back(std::bind(&ProteinProfile::Process_Pot_AAFreq, this, _1));
        vecCalcFuncs.push_back(std::bind(&ProteinProfile::Calculate_Pot_AAFreq_Profile, this));
        Msg += "Pot, ";
    }
    if (bSS) {
        vecProcFuncs.push_back(std::bind(&ProteinProfile::Process_SS, this, _1));
        vecCalcFuncs.push_back(std::bind(&ProteinProfile::Calculate_SS_Profile, this));
        Msg += "SS, ";
    }
    Msg.erase(Msg.size() - 2, 1);
    Msg += "profiles";
    std::cerr << Msg << std::endl;

    Thread_Manager(vecProcFuncs, _vecHomologous_Proteins, bVerbose, Msg);
    std::cerr << "done>1\njk???\noll";
    for (auto CalcFunction : vecCalcFuncs)
        CalcFunction();
    std::cerr << "famm\njkmmmm???";

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

    for (int i = 0; i < nThreads; ++i) {
        std::vector vecBatch(vecDB.begin() + i * size / nThreads,
                             vecDB.begin() + (i + 1) * (size) / nThreads);
        vecThreads.push_back(std::thread([ =, &nCount_Now] { Processing_Thread(vecFuncs, vecBatch, nCount_Now); })); // BUG #1: must send by copy
    }

    if (bVerbose) {
        Progress_Indicator(sMsg, 0, size);
        while (true) {
            Progress_Indicator(sMsg, nCount_Now, size);
            if (nCount_Now == size)
                break;
            usleep(500 * 1000);
        }
    }
    std::cerr << "Were out of threads\n";
    for (int i = 0; i < nThreads; ++i)
        vecThreads[i].join();
}



template<class FuncType, class Args>
void ProteinProfile::Processing_Thread(std::vector<std::function<FuncType> > vecFuncs,
                                       std::vector<Args> vecBatch, int& nCount_Now)
{
    int i = 0;
    for (auto& data : vecBatch) {
        i++;
        int j = 0;
        //   std::cerr << i << "\n";
        for (auto ProcFunction : vecFuncs) {
            ProcFunction(data);
            j++;
            //        std::cerr << "THRD " << i << " " << j << " " << vecBatch.size() << "\n";
        }
        std::lock_guard<std::mutex> lock(_mtx_count);
        nCount_Now++;
    }
}



void ProteinProfile::Process_IsHomologue(std::string fPath) {
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
    std::cerr << __PRETTY_FUNCTION__ << "\n";
    target.Calculate_Solvent();

    std::lock_guard<std::mutex> lock(_mtx_solvent);
    for (int i = 0, n = target.length(); i < n; ++i) {
        auto& amino_acid = target.vecAmino_Acid[i];
        if (!IsStandardAA(amino_acid.name))
            continue;

        int idx = sym_to_idx[amino_acid.symbol];
        _aSolvent_AA_Count[idx][ target.aSolvent_Accessibility[i] ]++;
        _aAA_Total_Count[idx]++;
    }
}



void ProteinProfile::Process_SS(Protein& target) {
    std::cerr << __PRETTY_FUNCTION__ << "\n";

    target.Calculate_SS();

    std::lock_guard<std::mutex> lock(_mtx_SS);
    std::cerr << target.length() << "\n";
    for (int i = 0, n = target.length(); i < n; ++i) {
        //    std::cerr << i << " " << target.fPath << "\n";
        //    std::cerr << target.sequence << " "  << "\n";
        //    std::cerr << sym_to_idx[ target[i]] << " " << target[i] << " " << target.aSecondary_Structure[i] << "\n";
        _aSec_AA_Count[ sym_to_idx[ target[i]] ][ target.aSecondary_Structure[i] ]++;
        //    std::cerr << "Done\n";
    }
    std::cerr << "Out of func\n";
}



void ProteinProfile::Process_Pot_AAFreq(Protein& target) {
    std::cerr << __PRETTY_FUNCTION__ << "\n";

    std::array<double, 20> AA_freq{}, pot_r, pot0, pot1{}, var0;
    int n = target.length();

    for (auto& amino_acid : target.vecAmino_Acid) {
        if (IsStandardAA(amino_acid.name))
            AA_freq[sym_to_idx[amino_acid.symbol]]++;
    }

    for (int i = 0; i < 20; ++i) {
        if (AA_freq[i] < 2)
            continue;
        pot_r[i] = std::exp(-std::sqrt(AA_freq[i]) / (n * _potS_Param));
        pot0[i] = (AA_freq[i] * (AA_freq[i] - 1) * pot_r[i] * (n - 1 / (1 - pot_r[i])))
                  / (n * (n - 1) * (1 - pot_r[i]));
        var0[i] = std::sqrt(std::pow(pot_r[i] * AA_freq[i] * (1 - AA_freq[i] / n), 2)
                            / ((1 - pot_r[i] * pot_r[i]) * n));
    }

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
        {
            if (!IsStandardAA(target.vecAmino_Acid[i].name))
                continue;
            int idx = sym_to_idx[ target[i] ];
            pot1[idx] += (target[i] == target[j]) * std::pow(pot_r[idx], j - i);
        }

    std::lock_guard<std::mutex> lock(_mtx_pot);
    for (int i = 0; i < 20; ++i)
    {
        if (AA_freq[i] < 2)
            continue;
        _vecPotScores[i].push_back((pot1[i] - pot0[i]) / var0[i]);
        _vecAA_Freqs[i].push_back(AA_freq[i] / n);
    }
}


void ProteinProfile::Calculate_Alignment_Profile(bool bSave_Frags, double dFrag_Score_Cutoff,
        double dGap_Score, int nMin_Frag)
{
    std::cerr << __PRETTY_FUNCTION__ << "\n";

    _dFrag_Score_Cutoff = dFrag_Score_Cutoff;
    _dGap_Score = dGap_Score;
    _nMin_Frag = nMin_Frag;
    std::cerr << "Rdy\n";
    for (auto& [score, seq1, seq2, atom_dist] : _vecTupleAlignments) {

        std::cerr << score << " " << seq1 << " " << seq2 << "\n";
        int align_len = seq1.length();
        std::vector<double> sum_dist( align_len + 1 );

        for (int i = 1; i <= align_len; ++i)
            sum_dist[i] = (atom_dist[i] == dGap_Score ? GAP_PENALTY : atom_dist[i])
                          + sum_dist[i - 1];


        for (int i = 0, pos = 0; i < align_len; ++i) {
            if (seq1[i] == '-')
                continue;

            pos++;
            if (seq2[i] == '-' || seq2[i] == 'X')
                continue;

            assert(algn_classes.count(seq2[i]) != 0);
            int nClass = algn_classes[seq2[i]];
            std::cerr << pos << " " << nClass << " " << _aAlgn_Profile.size() << " " << "\n";
            _aAlgn_Profile[pos - 1][nClass] += 10 - std::min(10., atom_dist[i + 1]); // Convert to similarity
            _aAlgn_Position_Count[pos - 1]++;
            _aAlgn_Class_Count[nClass]++;
        }
        std::cerr << "save frags\n";
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
                        std::cerr << "fragpush " << pos << " " << pos + nFragLen - 1 << " " << _refProtein.length() << "\n";
                        std::cerr << "poses " << i << " " << l << " " << nFragLen << " " << sc << " " << _dFrag_Score_Cutoff << "\n";
                        _matFragments[pos][pos + nFragLen - 1].push_back(fragment);
                    }
                }
                if (seq1[i] != '-')
                    pos++;
            } //
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
}



void ProteinProfile::Calculate_Solvent_Profile() {
    std::cerr << __PRETTY_FUNCTION__ << "\n";

    std::array<long long, 7> solvent_class_total{};

    long long nAA_Total = std::accumulate(_aAA_Total_Count.begin(), _aAA_Total_Count.end(), 0);

    for (int i = 0; i < 7; ++i)
        for (int j = 0; j < 20; ++j)
            std::cerr << _aSolvent_AA_Count[i][j] << " \n"[j == 19], solvent_class_total[i] += _aSolvent_AA_Count[j][i];

    for (int i = 0; i < 20; ++i)
        for (int j = 0; j < 7; ++j)
            if (solvent_class_total[j] != 0)
                _aSolvent_Profile[i][j] = (double)_aSolvent_AA_Count[i][j] / (solvent_class_total[j]
                                          * (_aAA_Total_Count[i] / (double)nAA_Total));
}



void ProteinProfile::Calculate_SS_Profile() {
    std::cerr << __PRETTY_FUNCTION__ << "\n";

    std::array<long long, 7> sec_class_total{};

    long long nAA_Total = std::accumulate(_aAA_Total_Count.begin(), _aAA_Total_Count.end(), 0);

    for (int i = 0; i < 7; ++i)
        for (int j = 0; j < 20; ++j)
            std::cerr << _aSec_AA_Count[j][i] << " \n"[j == 19], sec_class_total[i] += _aSec_AA_Count[j][i];

    for (int i = 0; i < 20; ++i)
        for (int j = 0; j < 7; ++j)
            if (sec_class_total[j] != 0)
                _aSec_Profile[i][j] = (double)_aSec_AA_Count[i][j] / (sec_class_total[j]
                                      * (_aAA_Total_Count[i] / (double)nAA_Total));
}



void ProteinProfile::Calculate_Pot_AAFreq_Profile() {
    std::cerr << __PRETTY_FUNCTION__ << "\n";

    for (int i = 0; i < 20; ++i) {
        std::cerr << "PotSize " << _vecPotScores[i].size() << "\n";
        if (_vecPotScores[i].size()) {
            _vecPotScores[i] = OutlierElimination_IQR(_vecPotScores[i]);
            std::cerr << "finished IQR\n";
            std::tie(_aPot_Bar[i], _aPot_Stdev[i]) = Mean_SD(_vecPotScores[i]);
        }

        if (_vecAA_Freqs[i].size()) {
            _vecAA_Freqs[i] = OutlierElimination_IQR(_vecAA_Freqs[i]);
            std::tie(_aAA_Freq_Mean[i], _aAA_Freq_Stdev[i]) = Mean_SD(_vecAA_Freqs[i]);
        }
    }
}


void ProteinProfile::Write_ToFile(bool bWriteCounts) {
    std::cerr << __PRETTY_FUNCTION__ << "\n";

    if (bSolvent_Rdy) {
        std::ofstream outFile(db_Profiles + RelativeFileName("solvent"));

        outFile << QuickInfo() << "\n\n" << "#SOLVENT_PROFILE\n\n"
                << sFileStartSym << "\n";

        for (int i = 0; i < 20; ++i) {
            outFile << idx_to_sym[i];
            for (int j = 0; j < 7; ++j)
                outFile << " " << std::fixed << std::setprecision(10) << _aSolvent_Profile[i][j];
            outFile << "\n";
        }
    }
    std::cerr << "CHECK1\n" << "\n";
    if (bSS_Rdy) {
        std::ofstream outFile(db_Profiles + RelativeFileName("sec"));

        outFile << QuickInfo() << "\n\n" << "#SECONDARY_STRUCTURE_PROFILE\n\n"
                << sFileStartSym << "\n";

        for (int i = 0; i < 20; ++i) {
            outFile << idx_to_sym[i];
            for (int j = 0; j < 7; ++j)
                outFile << " " << std::fixed << std::setprecision(10) << _aSec_Profile[i][j];
            outFile << "\n";
        }
    }
    std::cerr << "CHECK2\n" << "\n";

    if (bPot_Rdy) {
        std::ofstream outFile(db_Profiles + RelativeFileName("pot"));

        outFile << QuickInfo() << "\n\n" << "#POT_STATISTIC\n\n" << sFileStartSym << "\n";

        for (int i = 0; i < 20; ++i)
            outFile << idx_to_sym[i] << " " << std::fixed << std::setprecision(10)
                    << _aPot_Bar[i] << " " << _aPot_Stdev[i] << "\n";


        outFile << "\n\n\n#AA_FREQUENCY\n\n";
        for (int i = 0; i < 20; ++i)
            outFile << idx_to_sym[i] << " " << std::fixed << std::setprecision(10)
                    << _aAA_Freq_Mean[i] << " " << _aAA_Freq_Stdev[i] << "\n";

    }
    std::cerr << "CHECK3\n" << "\n";

    if (bAlgn_Rdy) {
        std::ofstream outFile(db_Profiles + RelativeFileName("alignment"));
        outFile << QuickInfo(true) << "\n\n" << "#ALIGNMENT_PROFILE\n\n" << sFileStartSym << "\n";

        for (int i = 0; i < _refProtein.length(); ++i) {
            outFile << std::left << std::setw(3) << i + 1;
            for (int j = 0; j < 7; ++j)
                outFile << " " << std::fixed << std::setprecision(10) << _aAlgn_Profile[i][j];
            outFile << "\n";
        }
        std::cerr << "CHECK4\n" << "\n";

        if (bFrags_Rdy) {
            std::ofstream outFile(db_Profiles + RelativeFileName("fragment"));
            outFile << QuickInfo(true) << "\n\n" << "#ALIGNMENT_FRAGMENTS\n\n" << sFileStartSym << "\n";

            for (int i = 0; i < _refProtein.length(); ++i)
                for (int j = i + _nMin_Frag - 1; j < _refProtein.length(); ++j)
                    for (auto& sFragment : _matFragments[i][j])
                        outFile << i << " " << j << " " << sFragment << "\n";

        }

    }

    if (bWriteCounts) {
        std::ofstream outFile(db_Profiles + RelativeFileName("count"));
        outFile << "#BURIAL_COUNT\n\n";
        for (int i = 0; i < 20; ++i) {
            outFile << idx_to_sym[i];
            for (int j = 0; j < 7; ++j)
                outFile << " " << std::setw(10) << _aSolvent_AA_Count[i][j];
            outFile << "\n";
        }

        outFile << "\n\n#STRUCTURE_COUNT\n\n";
        for (int i = 0; i < 20; ++i) {
            outFile << idx_to_sym[i];
            for (int j = 0; j < 7; ++j)
                outFile << " " << std::setw(10) << _aSec_AA_Count[i][j];
            outFile << "\n";
        }

        outFile << "\n\n#AA_COUNT\n\n";
        for (int i = 0; i < 20; ++i) {
            outFile << idx_to_sym[i]
                    << " " << _aAA_Total_Count[i] << "\n";
        }
    }
}


void ProteinProfile::Read_FromFile(std::string sParentDirectory) {
    std::string HEAD;

    std::string fSolvName = sParentDirectory + RelativeFileName("solvent");
    if (FileExists(fSolvName)) {
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

    std::string fSecName = sParentDirectory + RelativeFileName("sec");
    if (FileExists(fSolvName)) {
        std::ifstream inFile(fSecName);

        while (HEAD != sFileStartSym)
            inFile >> HEAD;

        for (int i = 0; i < 20; ++i) {
            inFile >> HEAD;
            for (int j = 0; j < 7; ++j)
                inFile >> _aSec_Profile[i][j];
        }

        bSS_Rdy = true;
    }

    std::string fPotName = sParentDirectory + RelativeFileName("pot");
    if (FileExists(fPotName)) {
        std::ifstream inFile(fPotName);

        while (HEAD != sFileStartSym)
            inFile >> HEAD;

        for (int i = 0; i < 20; ++i)
            inFile >> HEAD >> _aPot_Bar[i] >> _aPot_Stdev[i];


        inFile >> HEAD;
        for (int i = 0; i < 20; ++i)
            inFile >> HEAD >> _aAA_Freq_Mean[i] >> _aAA_Freq_Stdev[i];

        bPot_Rdy = true;
    }

    std::string fAlgnName = sParentDirectory + RelativeFileName("alignment");
    if (FileExists(fAlgnName)) {
        std::ifstream inFile(fAlgnName);

        while (HEAD != sFileStartSym)
            inFile >> HEAD;

        for (int i = 0; i < _refProtein.length(); ++i) {
            inFile >> HEAD;
            for (int j = 0; j < 7; ++j)
                inFile >> _aAlgn_Profile[i][j];
        }

        bAlgn_Rdy = true;

        std::string fFragName = sParentDirectory + RelativeFileName("fragment");
        if (FileExists(fFragName)) {
            std::ifstream inFile(fFragName);

            while (HEAD != sFileStartSym)
                inFile >> HEAD;

            for (int i, j; inFile >> i >> j >> HEAD; )
                _matFragments[i][j].emplace_back(HEAD);

            bFrags_Rdy = true;
        }
    }
}


std::string ProteinProfile::QuickInfo(bool bIncludeAlignment) {
    std::ostringstream output;
    output << "FILE_PATH:  " << _refProtein.fPath << "\n"
           << "MD5:  " << _refProtein.md5 << "\n"
           << "PROT_LEN:  " << _refProtein.length() << "\n"
           << "TEMPLATES_COUNT:  " << _vecHomologous_Proteins.size()
           << "ALGN_SCORE_CUTOFF:  " << _dAlgn_Score_CutOff;
    if (bIncludeAlignment) {
        output << "FRAG_SCORE_CUTOFF:  " << _dFrag_Score_Cutoff
               << "GAP_SCORE:  " << _dGap_Score
               << "MIN_FRAG_LENGTH" << _nMin_Frag;
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