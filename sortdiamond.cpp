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

int main(int argc, char *argv[]) {
	if (argc != 3) {
		cerr << "Usage: " << argv[0] << " <input_file> <output_file>"
		     << endl;
		return 1;
	}

	string in_name = argv[1];
	string ot_name = argv[2];

	ifstream infile(in_name);
	if (!infile) {
		cerr << "Error opening input file: " << in_name << endl;
		return 1;
	}

	ofstream outfile(ot_name);
	if (!outfile) {
		cerr << "Error opening output file: " << ot_name << endl;
		return 1;
	}

	map<string, pair<double, string>>
	    max_map;  // Key: sseqid, Value: pair<score, line>
	string line;

	while (getline(infile, line)) {
		istringstream iss(line);
		string fields[20];  // Adjust size if needed
		int i = 0;

		while (iss >> fields[i] && i < 20) {
			i++;
		}

		string key = fields[1];
		double value = stod(fields[11]);

		if (max_map.find(key) != max_map.end()) {
			if (value > max_map[key].first) {
				max_map[key] = make_pair(value, line);
			}
		} else {
			max_map[key] = make_pair(value, line);
		}
	}

	vector<pair<string, string>> result;

	for (const auto &entry : max_map) {
		istringstream iss(entry.second.second);
		string fields[20];  // Adjust size if needed
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

	for (const auto &entry : result) {
		outfile << entry.first << "\n" << entry.second << endl;
	}

	infile.close();
	outfile.close();

	return 0;
}

