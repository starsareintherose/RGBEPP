// #include <sys/stat.h>

#include <filesystem>
#include <fstream>
#include <iostream>
// #include <sstream>
#include <string>
// #include <vector>

namespace fs = std::filesystem;

std::string removeExtension(const std::string& filename) {
	size_t lastdot = filename.find_last_of(".");
	if (lastdot != std::string::npos) {
		return filename.substr(0, lastdot);
	}
	return filename;
}

void splitFasta(const std::string& input_fasta) {
	std::ifstream infile(input_fasta);
	if (!infile) {
		std::cerr << "Error: Unable to open input file " << input_fasta
			  << std::endl;
		return;
	}

	std::string line;
	std::string dir_name;
	std::ofstream outfile;
	bool in_sequence = false;

	while (std::getline(infile, line)) {
		if (line.empty()) continue;

		if (line[0] == '>') {
			// New sequence header
			if (in_sequence) {  // if found new sequence, close
				outfile.close();  // previous output file
			}
			dir_name = line.substr(1);  // Remove '>'
			// directory
			fs::create_directories(dir_name);
			fs::path output_file =
			    fs::path(dir_name) /
			    input_fasta;  // suitable to many os
			outfile.open(output_file);
			outfile << ">" << removeExtension(input_fasta)
				<< std::endl;
			// will enter sequence
			in_sequence = true;
		} else if (in_sequence) {
			// Inside sequence content
			outfile << line << std::endl;
		}
	}

	if (in_sequence) {
		outfile.close();
	}

	if (infile.eof()) {
		std::cout << "Sequences have been split into individual files."
			  << std::endl;
	} else {
		std::cerr << "Error occurred while reading file." << std::endl;
	}

	infile.close();
}

int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <input_fasta>"
			  << std::endl;
		return 1;
	}

	splitFasta(argv[1]);

	return 0;
}

