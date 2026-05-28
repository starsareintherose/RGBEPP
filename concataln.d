import std.stdio;
import std.string;
import std.array;
import std.conv;
import std.algorithm;
import std.exception;
import std.file;
import std.path;

struct GenePartition {
    string file;
    size_t start1; // 1-based inclusive
    size_t end1;   // 1-based inclusive
    size_t len;
}

struct GeneData {
    string file;
    size_t alnLen;
    string[string] seqBySpecies; // speciesKey -> sequence (no whitespace)
}

struct Supermatrix {
    string[] species;              // ordered list of species keys (include leading '>')
    size_t[string] speciesIndex;   // speciesKey -> index
    string[][] segByGeneThenSpecies; // [gene][species] => segment or ""
    size_t[] geneLens;             // [gene] => alnLen
    GenePartition[] partitions;    // ranges in concatenated matrix
    size_t totalLen;
}

bool isHeaderLine(string line)
{
    string s = line.stripLeft;
    return s.length > 0 && s[0] == '>';
}

string speciesKeyFromHeader(string headerLine)
{
    auto s = headerLine.strip;
    if (s.length == 0 || s[0] != '>')
        throw new Exception("Invalid FASTA header (does not start with '>'): " ~ headerLine);

    auto rest = s[1 .. $].stripLeft;
    enforce(rest.length > 0, "Empty FASTA id in header: " ~ headerLine);

    auto toks = rest.split(); // whitespace-delimited tokens
    enforce(toks.length > 0, "Empty FASTA id in header: " ~ headerLine);

    return ">" ~ toks[0]; // store with '>' for easy output
}

string stripSeqWhitespace(string line)
{
    char[] buf;
    buf.reserve(line.length);
    foreach (ch; line) {
        if (ch != ' ' && ch != '\t' && ch != '\r' && ch != '\n')
            buf ~= ch;
    }
    return cast(string)buf.idup;
}

GeneData readFastaFile(string file)
{
    enforce(exists(file), "File not found: " ~ file);

    string text = cast(string)readText(file);
    string[] lines = text.splitLines();

    string[string] seqsForFile;

    string currentSpecies = "";
    bool ignoreCurrent = false;
    char[] seqBuf;

    void flushRecord()
    {
        if (currentSpecies.length == 0) return;
        if (!ignoreCurrent) {
            seqsForFile[currentSpecies] = cast(string)seqBuf.idup;
        }
        seqBuf.length = 0;
    }

    foreach (line; lines) {
        if (isHeaderLine(line)) {
            flushRecord();

            currentSpecies = speciesKeyFromHeader(line);
            if (currentSpecies in seqsForFile) {
                ignoreCurrent = true;
                stderr.writeln("Found an extra copy of ", currentSpecies, " in file ", file, " ... ignoring this copy");
            } else {
                ignoreCurrent = false;
            }
        } else {
            if (!ignoreCurrent) {
                auto cleaned = stripSeqWhitespace(line);
                // append
                foreach (ch; cleaned) seqBuf ~= ch;
            }
        }
    }
    flushRecord();

    enforce(seqsForFile.length > 0, "No sequences parsed from file: " ~ file);

    // Determine alignment length and validate consistency
    size_t alnLen = size_t.max;
    foreach (sp, seq; seqsForFile) {
        if (alnLen == size_t.max) alnLen = seq.length;
        else enforce(seq.length == alnLen,
            "ERROR: not all sequences the same length in " ~ file ~
            " (" ~ to!string(seq.length) ~ " != " ~ to!string(alnLen) ~ ")");
    }
    enforce(alnLen != size_t.max, "Could not determine alignment length for file: " ~ file);

    return GeneData(file, alnLen, seqsForFile);
}

void addSpeciesIfMissing(ref Supermatrix sm, string sp)
{
    if (sp in sm.speciesIndex) return;
    sm.speciesIndex[sp] = sm.species.length;
    sm.species ~= sp;

    // IMPORTANT: when a new species appears, extend all existing gene rows
    foreach (ref row; sm.segByGeneThenSpecies) {
        row.length = sm.species.length;
    }
}

void addGene(ref Supermatrix sm, GeneData gene)
{
    // Ensure all species are registered
    foreach (sp; gene.seqBySpecies.keys) {
        addSpeciesIfMissing(sm, sp);
    }

    // Build a row (gene segments per species)
    string[] row = new string[](sm.species.length);
    foreach (sp, seq; gene.seqBySpecies) {
        row[sm.speciesIndex[sp]] = seq;
    }

    sm.segByGeneThenSpecies ~= row;
    sm.geneLens ~= gene.alnLen;

    size_t start1 = sm.totalLen + 1;
    size_t end1 = sm.totalLen + gene.alnLen;
    string gn = baseName(gene.file);
    string gnNoExt = stripExtension(gn); 
    sm.partitions ~= GenePartition(gnNoExt, start1, end1, gene.alnLen);
    sm.totalLen += gene.alnLen;
}

Supermatrix buildSupermatrix(string[] files)
{
    Supermatrix sm;
    sm.totalLen = 0;

    foreach (file; files) {
        GeneData gene = readFastaFile(file);
        addGene(sm, gene);
    }
    return sm;
}

void writePartitions(string path, const(GenePartition)[] parts)
{
    File fout = File(path, "w");
    foreach (p; parts) {
        fout.writefln("%s\t=\t%d-%d;", p.file, p.start1, p.end1);
    }
}

void writeSupermatrixFasta(string path, const Supermatrix sm)
{
    File fout = File(path, "w");

    foreach (spIdx, sp; sm.species) {
        fout.writeln(sp);

        foreach (geneIdx; 0 .. sm.segByGeneThenSpecies.length) {
            string seg = sm.segByGeneThenSpecies[geneIdx][spIdx];
            if (seg.length != 0) {
                enforce(seg.length == sm.geneLens[geneIdx],
                    "Internal error: segment length mismatch for " ~ sp);
                fout.write(seg);
            } else {
                fout.write(replicate("?", sm.geneLens[geneIdx]));
            }
        }

        fout.writeln();
    }
}

void main(string[] args)
{
    if (args.length < 3) {
	stderr.writeln("ConcatAln: Concatenate multiple sequence alignments");
        stderr.writeln("Usage: " ~ args[0] ~ " <prefix> file1 file2 ...");
        stderr.writeln("Example: " ~ args[0]  ~ " run1 gene*.fas");
        return;
    }

    string prefix = args[1];
    string[] files = args[2 .. $];

    Supermatrix sm = buildSupermatrix(files);

    string partPath = prefix ~ ".partitions.txt";
    string fasPath = prefix ~ ".fasta";

    writePartitions(partPath, sm.partitions);
    writeSupermatrixFasta(fasPath, sm);

    stdout.writeln("num_species = ", sm.species.length, ", num_files = ", files.length);
    stdout.writeln("expected concatenated sequence length ", sm.totalLen);
    stdout.writeln("Partitions written to \"", partPath, "\"");
    stdout.writeln("Concatenated alignment written to \"", fasPath, "\"");
}
