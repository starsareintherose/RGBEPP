/*
 * This script is for deleting stop codon
 * in both AA and original_NT files
 * Copyright
 * Guoyi Zhang, UNSW & Australian Museum
 */

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

class FastaProcessor {
      public:
	FastaProcessor(const std::string& inputFile,
		       const std::string& outputFile)
	    : input_file(inputFile), output_file(outputFile) {}

	void processFile() {
		// bakup input file
		std::string backup_file = input_file + ".bak";
		fs::copy_file(input_file, backup_file,
			      fs::copy_options::overwrite_existing);

		std::ifstream infile(input_file);
		std::ofstream outfile(input_file + ".tmp");  // temp file
		std::string line;

		while (std::getline(infile, line)) {
			if (line[0] == '>') {
				if (!sequence.empty()) {
					positions_map[taxa] =
					    processSequence(sequence);
					outfile << taxa << "\n"
						<< sequence << "\n";
				}
				taxa = line;
				sequence.clear();
			} else {
				sequence = line;
			}
		}

		// process final taxa
		if (!sequence.empty()) {
			positions_map[taxa] = processSequence(sequence);
			outfile << taxa << "\n" << sequence << "\n";
		}

		infile.close();
		outfile.close();

		// overwrite with temp file
		fs::rename(input_file + ".tmp", input_file);

		// process same name file
		if (fs::exists(output_file)) {
			// backup outputfile
			std::string output_backup_file = output_file + ".bak";
			fs::copy_file(output_file, output_backup_file,
				      fs::copy_options::overwrite_existing);

			processOutputFile();
			fs::rename(output_file + ".tmp", output_file);
		}
	}

      private:
	std::string input_file;
	std::string output_file;
	std::string taxa;
	std::string sequence;
	std::unordered_map<std::string, std::vector<size_t>> positions_map;

	std::vector<size_t> processSequence(std::string& sequence) {
		std::vector<size_t> positions;
		size_t non_dash_index = 0;

		for (size_t i = 0; i < sequence.size(); ++i) {
			if (sequence[i] != '-') {
				non_dash_index++;
				if (sequence[i] == '!' || sequence[i] == '*') {
					positions.push_back(
					    non_dash_index);  // pos without -
					sequence[i] = '-';    // replace with -
				}
			}
		}
		return positions;
	}

	void processOutputFile() {
		std::ifstream infile(output_file);
		std::ofstream outfile(output_file + ".tmp");
		std::string line;
		std::string current_taxa;
		std::string current_sequence;

		while (std::getline(infile, line)) {
			if (line[0] == '>') {
				if (!current_sequence.empty() &&
				    positions_map.find(current_taxa) !=
					positions_map.end()) {
					replacePositions(
					    current_sequence,
					    positions_map[current_taxa]);
					outfile << current_taxa << "\n"
						<< current_sequence << "\n";
				}
				current_taxa = line;
				current_sequence.clear();
			} else {
				current_sequence = line;
			}
		}

		// process final taxa
		if (!current_sequence.empty() &&
		    positions_map.find(current_taxa) != positions_map.end()) {
			replacePositions(current_sequence,
					 positions_map[current_taxa]);
			outfile << current_taxa << "\n"
				<< current_sequence << "\n";
		}

		infile.close();
		outfile.close();
	}

	void replacePositions(std::string& sequence,
			      const std::vector<size_t>& positions) {
		std::string new_sequence = sequence;
		size_t sequence_length = sequence.size();

		for (size_t pos : positions) {
			size_t start = (pos * 3) - 3;  //  3n-2
			size_t end = pos * 3;	       //  3n

			if (start < sequence_length) {
				for (size_t i = start;
				     i < end && i < sequence_length; ++i) {
					new_sequence[i] = '-';
				}
			}
		}

		sequence = new_sequence;
	}
};

int main(int argc, char* argv[]) {
	if (argc != 3) {
		std::cerr << "Usage: " << argv[0]
			  << " <input_folder> <output_folder>" << std::endl;
		return 1;
	}

	std::string input_folder = argv[1];
	std::string output_folder = argv[2];

	// get all fasta
	for (const auto& entry : fs::directory_iterator(input_folder)) {
		if (entry.path().extension() == ".fasta" ||
		    entry.path().extension() == ".fas") {
			std::string input_file = entry.path().string();
			std::string output_file =
			    output_folder + "/" +
			    entry.path().filename().string();

			FastaProcessor processor(input_file, output_file);
			processor.processFile();
		}
	}

	return 0;
}

