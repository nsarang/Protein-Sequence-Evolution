#include "deep.h"


std::vector<ProteinProfile> DeepAI::PrepareProfiles(std::string sFDir)
{
	std::vector<ProteinProfile> retProfiles;

	for (auto& eEntry : std::filesystem::directory_iterator(sFDir))
	{
		auto fPath = eEntry.path().string();
		auto baseName = eEntry.path().filename().string();

		if (baseName.rfind("Family", 0) != 0)
			continue;

		std::ifstream ifsFamily(fPath);
		std::string HEAD;
		ifsFamily >> HEAD >> HEAD;

		auto prot = Protein(HEAD);
		auto profile = ProteinProfile(prot);

		profile.Read_FromFile(db_Profiles, 1 + 2 + 4 + 8);
		if (profile.RemainingProfiles() != 0)
		{
			profile.Read_FromFile(sFDir, 32);
			profile.CalculateProfiles( profile.RemainingProfiles() );
			profile.Write_ToFile(db_Profiles, 1 + 2 + 4 + 8 + 32);
		}

		std::cerr << HEAD << " " << prot.md5 << "\n";
		retProfiles.push_back(profile);
	}

	return retProfiles;
}


void DeepAI::GenerateDataset(std::vector<ProteinProfile>& vecProfiles,
                             std::string fCSV,
                             int maxPerProt)
{
	std::ofstream outFile(fCSV);
	auto vecDB = utility::CATH_ListFiles(db_CATH);

	std::map<int, std::vector<Protein>> mpProtEqLen;
	std::map<std::string, int> mpUsedCount;

	for (auto& fPath : vecDB)
	{
		auto prot = Protein(fPath);
		mpProtEqLen[ prot.length() ].push_back(prot);
		mpUsedCount[ prot.md5 ] = 0;
	}
	for (auto& profile : vecProfiles)
	{
		int len = profile._refProtein.length();
		auto vecProtEq = mpProtEqLen[ len ];

		if (vecProtEq.size() < 3)
			continue;
		int nToGen = std::min((int)vecProtEq.size(), maxPerProt) - 1;
		auto bound = std::partition (vecProtEq.begin(),
		                             vecProtEq.end(),
		[&](Protein & p) { return mpUsedCount[ p.md5 ] < 10; }
		                            );

		decltype(vecProtEq) vecSample;
		std::sample(vecProtEq.begin(), bound,
		            std::back_inserter(vecSample),
		            maxPerProt + 1,
		            std::mt19937{std::random_device{}()});
		for (int i = 0; i < nToGen; ++i)
		{
			if (vecSample[i].md5 == profile._refProtein.md5) {
				nToGen++;
				continue;
			}
			GenerateScores(profile, vecSample[i], outFile);
			mpUsedCount[ vecSample[i].md5 ]++;
		}
		GenerateScores(profile, profile._refProtein, outFile);
		mpUsedCount[ profile._refProtein.md5 ]++;
	}
}


void DeepAI::GenerateScores(ProteinProfile& profile, Protein& prot,
                            std::ofstream& outFile, std::string sep)
{
	static DFIRE2 dDFIRE_Inst(db_DFIRE2);
	outFile << utility::FileBasename(profile._refProtein.fPath) << '-'
	        << utility::FileBasename(prot.fPath) << '-'
	        << prot.length() << sep;

	outFile << Evaluator::Solvent_Score(prot, profile._aSolvent_Profile) << sep
	        << Evaluator::Secondary_Struct(prot, profile._aSec_Profile) << sep
	        << Evaluator::AlignmentScore(prot, profile._aAlgn_Profile) << sep
	        << dDFIRE_Inst.Calc_CFE(prot) / prot.length() << sep;

	auto retPot = Evaluator::PotScore(prot,
	                                  profile._aPot_Bar,
	                                  profile._aPot_Stdev,
	                                  profile._dPotS_Param);
	for (auto sc : retPot)
		outFile << sc << sep;

	auto retFreq = Evaluator::FrequencyScore(prot, profile._aAA_Freq_Mean, profile._aAA_Freq_Stdev);
	for (auto sc : retFreq)
		outFile << sc << sep;

	outFile << PairAlgnScore(profile._refProtein.fPath, prot.fPath) << "\n";
}


double DeepAI::PairAlgnScore(std::string PDB_1, std::string PDB_2)
{
	std::string stdout, line;
	auto oBuffer = subprocess::check_output({ex_TMALIGN.c_str(),
	                                        PDB_1.c_str(),
	                                        PDB_2.c_str(), "-a"
	                                        });
	std::stringstream ret(oBuffer.buf.data());

	double score;
	getline(ret, line); getline(ret, line);
	ret >> line >> score;
	return score;
}