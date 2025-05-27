import std.stdio;
import std.file;
import std.algorithm;
import std.conv;
import std.string;
import std.array;
import std.typecons;

enum n_intarr = 5;
bool use_bitscore = true;

// Struct to store the data fields
struct SeqData {
    string sseqid;
    int qstart;
    int qend;
    double bit_or_e;
    string qseq;
}

// Function to generate reverse complement of a DNA sequence
string revcomp(const string seq) {
    string revseq;
    foreach_reverse (base; seq) {
        switch (base) {
            case 'A': revseq ~= 'T'; break;
            case 'T': revseq ~= 'A'; break;
            case 'C': revseq ~= 'G'; break;
            case 'G': revseq ~= 'C'; break;
            case 'R': revseq ~= 'Y'; break;
            case 'Y': revseq ~= 'R'; break;
            case 'S': revseq ~= 'S'; break;
            case 'W': revseq ~= 'W'; break;
            case 'K': revseq ~= 'M'; break;
            case 'M': revseq ~= 'K'; break;
            case 'B': revseq ~= 'V'; break;
            case 'D': revseq ~= 'H'; break;
            case 'H': revseq ~= 'D'; break;
            case 'V': revseq ~= 'B'; break;
            case 'N': revseq ~= 'N'; break;
            default: break; // 添加这一行
        }
    }
    return revseq;
}

// read input file and parse the data
void readInputFile(string filename, ref SeqData[] data_vector, int[] intnums, int intmax) {
    foreach(line; File(filename, "r").byLine) {
        auto fields = line.strip.split;
        if (fields.length <= intmax)
            continue;
        SeqData data;
        data.sseqid = fields[cast(size_t)intnums[0]].idup;
        data.qstart = fields[cast(size_t)intnums[1]].to!int;
        data.qend = fields[cast(size_t)intnums[2]].to!int;
        data.bit_or_e = fields[cast(size_t)intnums[3]].to!double;
        data.qseq = fields[cast(size_t)intnums[4]].idup;
        data_vector ~= data;
    }
}

// choose best bitscore/evalue
void processCompare(const SeqData[] data_vector, ref SeqData[string] best_map) {
    foreach (data; data_vector) {
        double bitE = data.bit_or_e;
        auto sseqid = data.sseqid;
        if (use_bitscore) {
            if (sseqid in best_map) {
                if (bitE > best_map[sseqid].bit_or_e)
                    best_map[sseqid] = data;
            } else {
                best_map[sseqid] = data;
            }
        } else {
            if (sseqid in best_map) {
                if (bitE < best_map[sseqid].bit_or_e)
                    best_map[sseqid] = data;
            } else {
                best_map[sseqid] = data;
            }
        }
    }
}

// check qst and qend, reverse complement if necessary
void processRevert(const SeqData[string] best_map, ref Tuple!(string, string)[] result) {
    foreach (entry; best_map.byKeyValue) {
        auto data = entry.value;
        string qseq = data.qseq;
        if (data.qstart > data.qend)
            qseq = revcomp(qseq);
        result ~= tuple(">" ~ data.sseqid, qseq);
    }
    result.sort!((a, b) => a[0] < b[0]);
}

// 写输出文件
void writeOutputFile(string filename, const Tuple!(string, string)[] result) {
    auto file = File(filename, "w");
    foreach (entry; result) {
        file.writeln(entry[0]);
        file.writeln(entry[1]);
    }
}

// use comma to split the input string into integers
void splitInts(string str, ref int[n_intarr] intnums) {
    auto parts = str.split(",");
    foreach (i, val; parts) if (i < n_intarr) intnums[i] = val.to!int;
}

// calculate the maximum column number from the array
int maxInts(int[n_intarr] intnums) {
    int intmax = intnums[0];
    foreach (i; 1 .. n_intarr)
        if (intnums[i] > intmax)
            intmax = intnums[i];
    return intmax;
}

void main(string[] args) {
    int[n_intarr] intnums = [1, 6, 7, 11, 17];
    int intmax = 17;

    if (args.length == 4 || args.length == 5) {
        splitInts(args[3], intnums);
        intmax = maxInts(intnums);
        if (args.length == 5) {
            string tmpstri = args[4];
            if (tmpstri == "bitscore")
                use_bitscore = true;
            else if (tmpstri == "evalue")
                use_bitscore = false;
            else
                writeln("Unknown argument: ", args[4]);
        }
    } else if (args.length != 3) {
        writeln("sortDiamond\nAuthor: Guoyi Zhang\nLicense:GPL-2.0-only\nUsage: ",
            args[0], " <input_file> <output_file> <sseq,qstart,qend,bitscore/evalue,qseq> <bitscore(default)/evalue>\nthe column number starts at 0");
        return;
    }

    if (args.length <= 5 && args.length >= 3) {
        string in_name = args[1];
        string ot_name = args[2];

        SeqData[] data_vector;
        readInputFile(in_name, data_vector, intnums[], intmax);

        SeqData[string] best_map;
        processCompare(data_vector, best_map);

        Tuple!(string, string)[] result;
        processRevert(best_map, result);

        writeOutputFile(ot_name, result);
    }
}
