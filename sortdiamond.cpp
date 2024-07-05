#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace std;

const int n_intarr = 5;
bool use_bitscore = true;

// Struct to store the data fields
struct SeqData {
	string sseqid;
	int qstart;
	int qend;
	double bit_or_e;
	string qseq;
};

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

void readInputFile(const string &filename, vector<SeqData> &data_vector,
		   const int *intnums, const int intmax) {
	ifstream infile(filename);
	if (!infile) {
		cerr << "Error opening input file: " << filename << endl;
		return;
	}
	string line;
	while (getline(infile, line)) {
		istringstream iss(line);
		vector<string> fields(intmax + 1);
		// string fields[intmax + 1];
		int i = 0;
		while (iss >> fields[i] && i < (intmax + 1)) {
			i++;
		}
		// Extract fields
		SeqData data;
		data.sseqid = fields[intnums[0]];
		data.qstart = stoi(fields[intnums[1]]);
		data.qend = stoi(fields[intnums[2]]);
		data.bit_or_e = stod(fields[intnums[3]]);
		data.qseq = fields[intnums[4]];

		// Store the data
		data_vector.push_back(data);
	}
	infile.close();
}

void processCompare(const vector<SeqData> &data_vector,
		    map<string, SeqData> &best_map) {
	for (const auto &data : data_vector) {
		double bitE = data.bit_or_e;
		const string &sseqid = data.sseqid;

		if (use_bitscore) {
			// bit score
			// Check if the seqid already exists in the map
			if (best_map.find(sseqid) != best_map.end()) {
				// If the new bitscore is greater, update the
				// map
				if (bitE > best_map[sseqid].bit_or_e) {
					best_map[sseqid] = data;
				}
			} else {
				// If the seqid does not exist, insert the new
				// seqid-bitE pair
				best_map[sseqid] = data;
			}
		} else {
			// evalue
			// Check if the seqid already exists in the map
			if (best_map.find(sseqid) != best_map.end()) {
				// If the new evalue is greater, update the
				// map
				if (bitE < best_map[sseqid].bit_or_e) {
					best_map[sseqid] = data;
				}
			} else {
				// If the seqid does not exist, insert the new
				// seqid-bitE pair
				best_map[sseqid] = data;
			}
		}
	}
}

void processRevert(const map<string, SeqData> &best_map,
		   vector<pair<string, string>> &result) {
	for (const auto &entry : best_map) {
		const SeqData &data = entry.second;

		// Check if qstart is larger than qend
		string qseq = data.qseq;
		if (data.qstart > data.qend) {
			qseq = revcomp(qseq);
		}
		result.push_back(make_pair(">" + data.sseqid, qseq));
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
	return intmax;
}

int main(int argc, char *argv[]) {
	int intnums[n_intarr] = {1, 6, 7, 11, 17};
	int intmax = 17;

	if (argc == 4 || argc == 5) {
		splitInts(argv[3], intnums);
		intmax = maxInts(intnums);
		if (argc == 5) {
			string tmpstri = argv[4];
			if (tmpstri == "bitscore") {
				use_bitscore = true;
			} else if (tmpstri == "evalue") {
				use_bitscore = false;
			} else {
				cout << "Unknown argument: " << argv[4] << endl;
			}
		}
	} else if (argc != 3) {
		cerr << "Usage: " << argv[0]
		     << " <input_file> <output_file> "
			"<sseq,qstart,qend,bitscore/evalue,qseq> "
			"<bitscore(default)/evalue>\nthe column number starts "
			"at 0"
		     << endl;
		return 1;
	}

	if (argc <= 5 && argc >= 3) {
		string in_name = argv[1];
		string ot_name = argv[2];

		// get useful information from reading input file
		vector<SeqData> data_vector;
		readInputFile(in_name, data_vector, intnums, intmax);

		// calculate best bitscore or evalue
		map<string, SeqData> best_map;
		processCompare(data_vector, best_map);

		// write the result
		vector<pair<string, string>> result;
		processRevert(best_map, result);

		writeOutputFile(ot_name, result);
	}

	return 0;
}

