import std.file;
import std.stdio;
import std.string;
import std.path;

string removeExtension(const string filename) {
    size_t lastdot = filename.lastIndexOf(".");
    if (lastdot != -1) {
        return filename[0 .. lastdot];
    }
    return filename;
}

void splitFasta(const string input_fasta) {
    File infile;
    try {
        infile = File(input_fasta, "r");
    } catch (FileException e) {
        stderr.writeln("Error: Unable to open input file ", input_fasta);
        return;
    }

    string line;
    string dir_name;
    File outfile;
    bool in_sequence = false;

    foreach (lineContent; infile.byLine()) {
        if (lineContent.empty) continue;

        if (lineContent[0] == '>') {
            // New sequence header
            if (in_sequence) {  // if found new sequence, close
                outfile.close();  // previous output file
            }
            dir_name = cast(string)lineContent[1 .. $];  // Remove '>'
            // directory
            mkdirRecurse(dir_name);
            auto output_file = buildPath(dir_name, input_fasta);  // suitable to many os
            outfile = File(output_file, "w");
            outfile.writeln(">", removeExtension(input_fasta));
            // will enter sequence
            in_sequence = true;
        } else if (in_sequence) {
            // Inside sequence content
            outfile.writeln(lineContent);
        }
    }

    if (in_sequence) {
        outfile.close();
    }

    if (infile.eof) {
        writeln("Sequences have been split into individual files.");
    } else {
        stderr.writeln("Error occurred while reading file.");
    }

    infile.close();
}

void main(string[] args) {
    if (args.length != 2) {
        stderr.writeln("splitFasta\nAuthor: Guoyi Zhang\nLicense:GPL-2.0-only\nUsage: ",
                       args[0], " <input_fasta>");
        return;
    }

    splitFasta(args[1]);
}
