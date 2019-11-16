#include "dataset.h"

void Dataset::CheckFamilies(std::string sFDir)
{
	auto vecDB = utility::CATH_ListFiles(db_CATH);
	std::set<std::string> setSearch;
	for (auto &prot : vecDB)
		setSearch.insert(prot);

	for (auto &eEntry : std::filesystem::directory_iterator(sFDir))
	{
		auto fPath = eEntry.path().string();
		auto baseName = eEntry.path().filename().string();

		if (baseName.rfind("Family", 0) != 0)
			continue;

		std::ifstream ifsFamily(fPath);
		std::string HEAD, NAME;
		ifsFamily >> HEAD >> NAME;

		utility::Progress_Indicator(NAME, 0, 1);
		while (HEAD != sFileStartSym)
			ifsFamily >> HEAD;

		while (std::getline(ifsFamily, HEAD))
		{
			if (HEAD != "")
			{
				if (!utility::FileExists(HEAD))
				{
					std::cerr << HEAD << "\tNOT FOUND!\n";
					exit(1);
				}
				else if (setSearch.count(HEAD) == 0)
				{
					std::cout << "Problem in " << fPath << ":" << NAME << "\t"
							  << "File: " << HEAD << "\n";
				}
			}
		}
		utility::Progress_Indicator(NAME, 1, 1);
	}
}

void Dataset::GenerateDataset(std::string sFamDir,
							  std::string fCSV,
							  int nFamCutoff)
{
	std::ofstream trFile(fCSV);
	trFile.close(); // Truncate file
	auto vecDB = utility::CATH_ListFiles(db_CATH);
	ThreadPool pool(10); // Creates 10 threads.

	std::map<int, std::vector<Protein>> mpProtEqLen;
	int countNow = 0;

	// prepare maps
	for (auto &fPath : vecDB)
	{
		utility::Progress_Indicator("Loading protein DB", countNow++, vecDB.size());
		auto prot = Protein(fPath, 0, Pot_S_Constant, false);
		mpProtEqLen[prot.length()].push_back(prot);
	}
	utility::Progress_Indicator("Loading protein DB", 1, 1);

	int totProfs = utility::CountFilesInDir(sFamDir);
	countNow = 0;

	for (auto &eEntry : std::filesystem::directory_iterator(sFamDir))
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
			profile.CalculateProfiles(profile.RemainingProfiles());
			profile.Write_ToFile(db_Profiles, 1 + 2 + 4 + 8 + 32);
		}

		if (profile.FamilySize() < nFamCutoff)
			continue;

		std::ofstream outFile(fCSV, std::ios_base::app);
		for (auto protein : mpProtEqLen[prot.length()])
			pool.AddJob([this, protein, &outFile, &profile] { GenerateData(profile, protein, outFile); });

		pool.WaitAll();
	}
	utility::Progress_Indicator("Generating dataset", 1, 1);
}


/**
 * @brief 				Engineer the features to be given as input to the ML model.
 * 
 * @param profile 		Profile of the target protein structure
 * @param candidProt 	A protein that it's sequence is used as a candidate. It's 
 * 					    structure is also used to generate the groundtruth (aka
 * 						TM-Score value)
 * @param outFile 		Buffer to write the values to
 * @param sep 			The seperator used to split the values
 */
void Dataset::GenerateData(ProteinProfile &profile, Protein candidProt,
						   std::ostream &outFile, std::string sep)
{
	assert(profile._refProtein.length() == candidProt.length());

	auto newProt = profile._refProtein; // make a copy of the target structure
	newProt.sequence = candidProt.sequence; // replace with the new sequence
	newProt.Calculate_Pot(Pot_S_Constant, true);

	static DFIRE2 dDFIRE_Inst(db_DFIRE2);
	static std::mutex _mtx_write;
	std::stringstream outSS;

	outSS << utility::FileBasename(newProt.fPath) << sep
		  << utility::FileBasename(candidProt.fPath) << sep
		  << newProt.length() << sep;

	outSS << Evaluator::Solvent_Score(newProt, profile._aSolvent_Profile) << sep
		  << Evaluator::Secondary_Struct(newProt, profile._aSec_Profile) << sep
		  << Evaluator::AlignmentScore(newProt, profile._aAlgn_Profile) << sep
		  << dDFIRE_Inst.Calc_CFE(newProt) / newProt.length() << sep;

	auto retPot = Evaluator::PotScore(newProt,
									  profile._aPot_Bar,
									  profile._aPot_Stdev,
									  profile._dPotS_Param);
	for (auto sc : retPot)
		outSS << sc << sep;

	auto retFreq = Evaluator::FrequencyScore(newProt, profile._aAA_Freq_Mean, profile._aAA_Freq_Stdev);
	for (auto sc : retFreq)
		outSS << sc << sep;

	outSS << PairAlgnScore(newProt.fPath, candidProt.fPath) << "\n";
	std::lock_guard<std::mutex> lock(_mtx_write);
	outFile << outSS.rdbuf();
}

double Dataset::PairAlgnScore(std::string PDB_1, std::string PDB_2)
{
	std::string stdout, line;
	auto oBuffer = subprocess::check_output({ex_TMALIGN.c_str(),
											 PDB_1.c_str(),
											 PDB_2.c_str(), "-a"});
	std::stringstream ret(oBuffer.buf.data());

	double score;
	getline(ret, line);
	getline(ret, line);
	ret >> line >> score;
	return score;
}