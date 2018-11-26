
#include <iostream>
#include <tuple>
#include <vector>
#include <fstream>
#include <dirent.h>
#include <cassert>
#include <sstream>
#include <array>
#include <mutex>
#include <thread>


#define RESET   "\033[0m"
#define BOLDRED     "\033[1m\033[31m"
#define BOLDGREEN   "\033[1m\033[32m"
#define BUFFERSIZE 1024
std::string ex_STRIDE;



// TYPE DEFINITIONS
struct AminoAcid {
	AminoAcid(std::string iname, double x = 0, double y = 0, double z = 0,
	          int count = 0)
		: name( iname ), neighbour_count( count ), cords( std::tie(x, y, z)) {}

	std::string name;
	int neighbour_count;
	std::tuple<double, double, double> cords;
};




void ProgressIndicator(std::string text, long long current, long long total) {
	std::cout << (current != 0 ?  "\r\033[F" : "" ) << text << ": ";
	int percent = current * 100 / total;
	if (percent == 100)
		std::cout << BOLDGREEN << "OK " << RESET << std::endl;
	else
		std::cout << BOLDRED << percent << "%" << RESET << std::endl;
}


int SystemCall_Err(std::string command, std::string& stdout) {
	std::array<char, BUFFERSIZE> buffer;

	FILE* pipe = popen(command.c_str(), "r");
	if (!pipe) throw std::runtime_error("popen() failed!\ncommand: " + command);
	try {
		while (!feof(pipe))
			if (fgets(buffer.data(), BUFFERSIZE, pipe) != NULL)
				stdout += buffer.data();
	} catch (...) {
		pclose(pipe);
		throw;
	}
	int wstat = pclose(pipe);
	return WEXITSTATUS(wstat);
}


bool FileExists(std::string fName) {
	std::ifstream inFile(fName);
	return inFile.good();
}


bool IsStandardAA(std::string abrv)
{
	return !(abrv == "ASX" || abrv == "GLX" || abrv == "SEC" || abrv == "PYL" || abrv == "UNK");
}


std::tuple<std::string, std::string> SecondaryStructure(std::string fPath) {
	std::string stdout, line, prev_line,
	    seq_ret, sec_str;
	while (SystemCall_Err(ex_STRIDE + " -o " + fPath + " 2>&1", stdout) != 0);
	std::stringstream ret(stdout);

	while (getline(ret, line)) {
		if (line.substr(0, 3) == "STR")
		{
			for (int i = 10; i < 60; ++i)
			{
				if (isspace(prev_line[i]) || prev_line[i] == 'X')
					continue;
				if (isspace(line[i]))
					line[i] = 'C';
				if (line[i] == 'b')
					line[i] = 'B';

				seq_ret += prev_line[i];
				sec_str += line[i];
			}
		}
		if (line.substr(0, 3) == "LOC")
			break;

		prev_line = line;
	}
	return tie(seq_ret, sec_str);
}


double DataNoise(std::string sequence) {
	int count = 0;
	for (auto c : sequence) {
		if (c == 'X')
			count++;
	}
	return (double)count / sequence.length();
}


void DB_ListFiles(std::string path, std::vector<std::string>& db_protein_list) {
	db_protein_list.clear();
	DIR *hDir;
	dirent *hFile;
	assert(hDir = opendir(path.c_str()));
	path += (path.back() != '/' ? '/' : 0);

	while ((hFile = readdir(hDir))) {
		std::string fName = hFile->d_name;
		if (fName.size() != 7 || !isdigit(fName[0])) // Characteristics of CATH database
			continue;

		db_protein_list.push_back(path + fName);
	}
}


void pThread(int s, int t, std::vector<std::string>& ProteinList, std::ofstream& outFile) {
	static int count = 0;
	static std::mutex db_mtx_print;
	std::string seq, sec_str;

	for (int i = s; i < t; ++i)
	{
		tie(seq, sec_str) = SecondaryStructure(ProteinList[i]);
		std::lock_guard<std::mutex> lock(db_mtx_print);
		if (DataNoise(seq) < .1 )	
			outFile << seq << " " << sec_str << "\n";
		ProgressIndicator("Processing", ++count, ProteinList.size());
	}
}


void Parse_Database(std::string sDB_Path, std::string outPath, int thread_num = 12) {
	std::ofstream outFile(outPath, std::ios::trunc);
	std::vector<std::string> ProteinList;
	std::vector<std::thread> vecThreads;

	DB_ListFiles(sDB_Path, ProteinList);
	int size = ProteinList.size();

	ProgressIndicator("Processing", 0, 1);
	for (int i = 0; i < thread_num; ++i)
		vecThreads.push_back(std::thread(pThread,
		                                 i * size / thread_num, (i + 1) * (size) / thread_num,
		                                 std::ref(ProteinList),
		                                 std::ref(outFile)
		                                ));

	for (int i = 0; i < thread_num; ++i)
		vecThreads[i].join();
}


int main(int argc, const char * argv[]) {
	if (argc != 4) {
		std::cerr << "Usage: parser database_path stride_path output_name\n";
		return 0;
	}
	ex_STRIDE = std::string(argv[2]);
	Parse_Database(std::string(argv[1]), std::string(argv[3]));
	std::cerr << "Done.\n";
}