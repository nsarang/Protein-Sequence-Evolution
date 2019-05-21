#include "deep.h"


void DeepAI::CheckFamilies(std::string sFDir)
{
	auto vecDB = utility::CATH_ListFiles(db_CATH);
	std::set<std::string> setSearch;
	for (auto& prot : vecDB)
		setSearch.insert(prot);

	for (auto& eEntry : std::filesystem::directory_iterator(sFDir))
	{
		auto fPath = eEntry.path().string();
		auto baseName = eEntry.path().filename().string();

		if (baseName.rfind("Family", 0) != 0)
			continue;

		std::ifstream ifsFamily(fPath);
		std::string HEAD, NAME;
		ifsFamily >> HEAD >> NAME;

		// utility::Progress_Indicator(NAME, 0, 1);
		while (HEAD != sFileStartSym)
			ifsFamily >> HEAD;

		while (std::getline(ifsFamily, HEAD)) {
			if (HEAD != "") {
				if (!utility::FileExists(HEAD)) {
					std::cerr << HEAD << "\tNOT FOUND!\n";
					exit(1);
				}
				else if (setSearch.count(HEAD) == 0) {
					std::cout << "Problem in " << fPath << ":" << NAME << "\t" << "File: " << HEAD << "\n";
					// std::cout << "Problem in " << fPath << "\t" << "File: " << HEAD << "\n";

				}
			}
		}
		// utility::Progress_Indicator(NAME, 1, 1);
	}
}


void DeepAI::GenerateDataset(std::string sFamDir,
                             std::string fCSV,
                             int nFamCutoff)
{
	std::ofstream trFile(fCSV); trFile.close(); // Truncate file
	auto vecDB = utility::CATH_ListFiles(db_CATH);
	ThreadPool pool(10); // Creates 10 threads.

	std::map<int, std::vector<Protein>> mpProtEqLen;
	int countNow = 0;

	// prepare maps
	for (auto& fPath : vecDB)
	{
		utility::Progress_Indicator("Loading protein DB", countNow++, vecDB.size());
		auto prot = Protein(fPath, 0 , Pot_S_Constant, false);
		mpProtEqLen[ prot.length() ].push_back(prot);
	}
	utility::Progress_Indicator("Loading protein DB", 1, 1);


	int totProfs = utility::CountFilesInDir(sFamDir);
	countNow = 0;

	for (auto& eEntry : std::filesystem::directory_iterator(sFamDir))
	{
		utility::Progress_Indicator("Generating dataset", countNow++, totProfs);
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
			profile.Read_FromFile(sFamDir, 32);
			profile.CalculateProfiles( profile.RemainingProfiles() );
			profile.Write_ToFile(db_Profiles, 1 + 2 + 4 + 8 + 32);
		}

		if (profile.FamilySize() < nFamCutoff)
			continue;

		std::ofstream outFile(fCSV, std::ios_base::app);
		for (auto protein : mpProtEqLen[ prot.length() ])
			pool.AddJob( [this, protein, &outFile, &profile] { GenerateData(profile, protein, outFile); } );

		pool.WaitAll();
	}
	utility::Progress_Indicator("Generating dataset", 1, 1);
}


void DeepAI::GenerateData(ProteinProfile& profile, Protein prot,
                          std::ofstream& outFile, std::string sep)
{
	static DFIRE2 dDFIRE_Inst(db_DFIRE2);
	static std::mutex _mtx_write;
	std::stringstream outSS;

	outSS << utility::FileBasename(profile._refProtein.fPath) << sep
	      << utility::FileBasename(prot.fPath) << sep
	      << prot.length() << sep;

	outSS << Evaluator::Solvent_Score(prot, profile._aSolvent_Profile) << sep
	      << Evaluator::Secondary_Struct(prot, profile._aSec_Profile) << sep
	      << Evaluator::AlignmentScore(prot, profile._aAlgn_Profile) << sep
	      << dDFIRE_Inst.Calc_CFE(prot) / prot.length() << sep;

	auto retPot = Evaluator::PotScore(prot,
	                                  profile._aPot_Bar,
	                                  profile._aPot_Stdev,
	                                  profile._dPotS_Param);
	for (auto sc : retPot)
		outSS << sc << sep;

	auto retFreq = Evaluator::FrequencyScore(prot, profile._aAA_Freq_Mean, profile._aAA_Freq_Stdev);
	for (auto sc : retFreq)
		outSS << sc << sep;

	outSS << PairAlgnScore(profile._refProtein.fPath, prot.fPath) << "\n";
	std::lock_guard<std::mutex> lock(_mtx_write);
	outFile << outSS.rdbuf();
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