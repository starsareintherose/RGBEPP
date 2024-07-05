#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace std;

// Function to generate reverse complement of a DNA sequence
string revcomp(const string &seq) {
	string revseq;
	for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
		switch (*it) {
			case 'A':
				revseq += 'T';
				break;
			case 'T':
				revseq += 'A';
				break;
			case 'C':
				revseq += 'G';
				break;
			case 'G':
				revseq += 'C';
				break;
			case 'R':
				revseq += 'Y';
				break;
			case 'Y':
				revseq += 'R';
				break;
			case 'S':
				revseq += 'S';
				break;
			case 'W':
				revseq += 'W';
				break;
			case 'K':
				revseq += 'M';
				break;
			case 'M':
				revseq += 'K';
				break;
			case 'B':
				revseq += 'V';
				break;
			case 'D':
				revseq += 'H';
				break;
			case 'H':
				revseq += 'D';
				break;
			case 'V':
				revseq += 'B';
				break;
			case 'N':
				revseq += 'N';
				break;
			default:
				break;
		}
	}
	return revseq;
}

void readInputFile(const string &filename,
		   map<string, pair<double, string>> &max_map) {
	ifstream infile(filename);
	if (!infile) {
		cerr << "Error opening input file: " << filename << endl;
		return;
	}

	string line;
	while (getline(infile, line)) {
		istringstream iss(line);
		string fields[20];
		int i = 0;
		while (iss >> fields[i] && i < 20) {
			i++;
		}
		string key = fields[1];
		double value = stod(fields[11]);
		// Check if the key already exists in the map
		if (max_map.find(key) != max_map.end()) {
			// If the new value is greater, update the map
			if (value > max_map[key].first) {
				max_map[key] = make_pair(value, line);
			}
		} else {
			// If the key does not exist, insert the new key-value
			// pair
			max_map[key] = make_pair(value, line);
		}
	}
	infile.close();
}

void processMap(const map<string, pair<double, string>> &max_map,
		vector<pair<string, string>> &result) {
	for (const auto &entry : max_map) {
		istringstream iss(entry.second.second);
		string fields[20];
		int i = 0;
		while (iss >> fields[i] && i < 20) {
			i++;
		}
		if (stoi(fields[6]) > stoi(fields[7])) {
			fields[17] = revcomp(fields[17]);
		}
		result.push_back(make_pair(">" + fields[1], fields[17]));
	}
	sort(result.begin(), result.end());
}

void writeOutputFile(const string &filename,
		     const vector<pair<string, string>> &result) {
	ofstream outfile(filename);
	if (!outfile) {
		cerr << "Error opening output file: " << filename << endl;
		return;
	}

	for (const auto &entry : result) {
		outfile << entry.first << "\n" << entry.second << "\n";
	}
	outfile.close();
}

int main(int argc, char *argv[]) {
	if (argc != 3) {
		cerr << "Usage: " << argv[0] << " <input_file> <output_file>"
		     << endl;
		return 1;
	}

	string in_name = argv[1];
	string ot_name = argv[2];

	map<string, pair<double, string>> max_map;
	readInputFile(in_name, max_map);

	vector<pair<string, string>> result;
	processMap(max_map, result);

	writeOutputFile(ot_name, result);

	return 0;
}

