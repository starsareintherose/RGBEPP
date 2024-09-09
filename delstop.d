import std.file;
import std.stdio;
import std.string;
import std.algorithm;
import std.array;
import std.conv;

// contains taxa and pos
struct TaxaInfo {
    string name;
    int[] positions;
}

// include TaxaInfo and foundSpecialChar
struct FastaData {
    TaxaInfo[] taxaList;  // every taxa correspeonds to TaxaInfo
    bool foundSpecialChar;
}

void backup(string filePath) {
    string backupPath = filePath ~ ".bak";
    copy(filePath, backupPath);
}

string replaceSpecialChars(string sequence, const(int)[] positions) {
    char[] mutableSequence = sequence.dup;  // change to editable
    foreach (pos; positions) {
        mutableSequence[pos - 1] = '-';  // replace `!` and `*` with `-`
    }
    sequence = mutableSequence.idup; // change to not editable
    return sequence;
}

int[] findSpecialPositions(const string sequence) {
    int[] positions;
    foreach (i, char c; sequence) {
        if (c == '*' || c == '!') {
            positions ~= cast(int) (i + 1); // record `int` pos
        }
    }
    return positions;
}

FastaData processFastaAA(string filePath) {
    FastaData fastaData;  // create struct
    string[] lines = File(filePath).byLine().map!(line=>line.to!string).array;  // read file
    string currentTaxa;
    string currentSequence;
    string tmpContent;

    for (int i = 0; i < lines.length; i++) {
        if (lines[i].startsWith(">")) {
            if (!currentTaxa.empty) {
                // keep last sequence
                TaxaInfo taxaInfo;
                taxaInfo.name = currentTaxa;
                taxaInfo.positions = findSpecialPositions(currentSequence);
                fastaData.taxaList ~= taxaInfo;  // add TaxaInfo to taxaList 
                
                if (taxaInfo.positions.length > 0) {
                    fastaData.foundSpecialChar = true;  // check any special char
                }

                currentSequence = replaceSpecialChars(currentSequence, taxaInfo.positions);  // replace
		tmpContent ~= ">" ~ currentTaxa ~ "\n" ~ currentSequence ~ "\n";
                currentSequence = "";  // clean prepare next
            }
            currentTaxa = lines[i][1..$];  // get taxa
        } else {
            currentSequence ~= lines[i];  // get seq
        }
    }

    // the final one
    if (!currentTaxa.empty) {
        TaxaInfo taxaInfo;
        taxaInfo.name = currentTaxa;
        taxaInfo.positions = findSpecialPositions(currentSequence);
        fastaData.taxaList ~= taxaInfo;

        if (taxaInfo.positions.length > 0) {
            fastaData.foundSpecialChar = true;  
        }

        currentSequence = replaceSpecialChars(currentSequence, taxaInfo.positions);  
    }
    if(fastaData.foundSpecialChar){
        backup(filePath);
        // write back
        std.file.write(filePath, tmpContent);
    }
    return fastaData; 
}

void processFastaNT(string filePath, const FastaData fastaData) {
    string[] lines = File(filePath).byLine().map!(line=>line.to!string).array; // get fasta_nt 

    for (int i = 0; i < lines.length; i++) {
        if (lines[i].startsWith(">")) {
            string currentTaxa = lines[i][1..$];
            auto taxaIndex = fastaData.taxaList.countUntil!(t => t.name == currentTaxa);  // use countUntil find taxaName

            if (taxaIndex != fastaData.taxaList.length) {  // if find taxa
                int[] positions = fastaData.taxaList[taxaIndex].positions.dup;  

                int lineIndex = i + 1;
                string sequence;
                while (lineIndex < lines.length && !lines[lineIndex].startsWith(">")) {
                    sequence ~= lines[lineIndex];  // get full
                    lineIndex++;
                }

                // replace
                char[] mutableSequence = sequence.dup;
		foreach (pos; positions) {
    		    int startPos = 3 * pos - 2;  //  3n-2 
    		    int endPos = 3 * pos;        //  3n 
    		    if (startPos - 1 < mutableSequence.length) {
        		for (int j = startPos - 1; j < endPos && j < mutableSequence.length; j++) {
            			mutableSequence[j] = '-';  //  3n-2 to 3n  '-'
        		}
    		    }
		}
                sequence = mutableSequence.idup;

                // write to fasta_nt
                int sequencePos = 0;
                lineIndex = i + 1;
                while (lineIndex < lines.length && !lines[lineIndex].startsWith(">")) {
                    int len = cast(int) lines[lineIndex].length;  // convert to int
                    lines[lineIndex] = sequence[sequencePos..sequencePos + len];
                    sequencePos += len;
                    lineIndex++;
                }
            }
        }
    }
    // write back
    std.file.write(filePath, lines.join("\n"));  
}

void processFasta(string fasta_aa, string fasta_nt){
    // get pos from fasta_aa & modify
    FastaData fastaData = processFastaAA(fasta_aa);

    if (fastaData.foundSpecialChar) {
        // if found special, proces fasta_nt
        // backup fasta_aa
        backup(fasta_nt);  // modify fasta_nt
        processFastaNT(fasta_nt, fastaData);  // modify fasta_nt
    } 
}

void main(string[] args) {
    if (args.length < 3) {
        writeln("Usage: program <fasta_aa> <fasta_nt>");
        return;
    }

    string fasta_aa = args[1];
    string fasta_nt = args[2];

    processFasta(fasta_aa, fasta_nt);
}

