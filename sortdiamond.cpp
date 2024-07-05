#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace std;

const int n_intarr = 5;

string revcomp(const string &seq);
int maxInts(int intnums[n_intarr]);

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
		   map<string, pair<double, string>> &max_map,
		   const int *intnums, const int intmax) {
	ifstream infile(filename);
	if (!infile) {
		cerr << "Error opening input file: " << filename << endl;
		return;
	}

	string line;
	while (getline(infile, line)) {
		istringstream iss(line);
		string fields[intmax + 1];
		int i = 0;
		while (iss >> fields[i] && i < (intmax + 1)) {
			i++;
		}
		// subject seq id
		string key = fields[intnums[0]];
		// bit score
		double value = stod(fields[intnums[3]]);
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
		vector<pair<string, string>> &result, const int *intnums,
		const int intmax) {
	for (const auto &entry : max_map) {
		istringstream iss(entry.second.second);
		string fields[intmax + 1];
		int i = 0;
		while (iss >> fields[i] && i < (intmax + 1)) {
			i++;
		}
		// check if qstart is larger than qend
		if (stoi(fields[intnums[1]]) > stoi(fields[intnums[2]])) {
			fields[intnums[4]] = revcomp(fields[intnums[4]]);
		}
		result.push_back(
		    make_pair(">" + fields[intnums[0]], fields[intnums[4]]));
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

void splitInts(const std::string &str, int intnums[n_intarr]) {
	istringstream iss(str);
	int numi = 0;
	string tmpstr;
	while (std::getline(iss, tmpstr, ',') && numi < n_intarr) {
		intnums[numi] = std::stoi(tmpstr);
		cout << intnums[numi] << endl;
		numi++;
	}
}

int maxInts(int intnums[n_intarr]) {
	int intmax = intnums[0];
	for (int i = 1; i < n_intarr; i++) {
		if (intnums[i] > intmax) {
			intmax = intnums[i];
		}
	}
	cout << intmax << endl;
	return intmax;
}

int main(int argc, char *argv[]) {
	int intnums[n_intarr] = {1, 6, 7, 11, 17};
	int intmax = 17;

	if (argc == 4) {
		splitInts(argv[3], intnums);
		intmax = maxInts(intnums);

	} else if (argc != 3) {
		cerr << "Usage: " << argv[0]
		     << " <input_file> <output_file> "
			"<sseq,qstart,qend,bitscore,qseq>"
		     << endl;
		return 1;
	}

	if (argc == 3 || argc == 4) {
		string in_name = argv[1];
		string ot_name = argv[2];

		map<string, pair<double, string>> max_map;
		readInputFile(in_name, max_map, intnums, intmax);

		vector<pair<string, string>> result;
		processMap(max_map, result, intnums, intmax);

		writeOutputFile(ot_name, result);
	}

	return 0;
}

