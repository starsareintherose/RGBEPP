import std.stdio;
import std.file;
import std.path;
import std.algorithm;
import std.array;
import std.string;
import std.getopt;

/*
  Minimal FASTA parser:
  - record header line starts with '>'
  - sequence lines until next '>' or EOF
  - we remove a record if header contains any taxa substring
*/

struct FastaRecord {
    string header; // without '>'
    string seq;    // includes newlines as stored (normalized to '\n')
}

FastaRecord[] readFasta(string path) {
    string text = cast(string) std.file.read(path);
    text = text.replace("\r\n", "\n").replace("\r", "\n");

    FastaRecord[] recs;
    string curHeader;
    auto seqBuf = appender!string();
    bool hasSeqLine = false;

    foreach (rawLine; text.splitter('\n')) {
        string line = rawLine.stripLeft();

        if (line.length && line[0] == '>') {
            if (curHeader.length) {
                recs ~= FastaRecord(curHeader, seqBuf.data);
                seqBuf = appender!string();
                hasSeqLine = false;
            }
            curHeader = line[1 .. $];
        } else {
            if (curHeader.length) {
                // skip empty lines in sequence, which can be common in some FASTA files
                if (line.length == 0) continue;

                // only add newline if we already have sequence lines, to avoid leading newline
                if (hasSeqLine) seqBuf.put("\n");
                seqBuf.put(line);
                hasSeqLine = true;
            }
        }
    }

    if (curHeader.length) {
        recs ~= FastaRecord(curHeader, seqBuf.data);
    }
    return recs;
}

void writeFasta(string path, const(FastaRecord)[] recs) {
    auto buf = appender!string();

    foreach (i, r; recs) {
        buf.put(">");
        buf.put(r.header);
        buf.put("\n");
        buf.put(r.seq);

        // finall record doesn't need extra newline
        if (i + 1 < recs.length) buf.put("\n");
    }

    std.file.write(path, buf.data);
}

bool shouldRemove(string header, string[] taxa) {
    foreach (t; taxa) {
        if (t.length == 0) continue;
        if (header.canFind(t)) return true;
    }
    return false;
}

int main(string[] args) {
    string fastaPath;
    string[] taxa;
    if (args.length < 3) {
        stderr.writeln("Usaege: " ~ args[0] ~ " <xxx.fasta> <taxa1> [taxa2 ...]");
        return 2;
    } else {
        fastaPath = args[1];
        taxa = args[2 .. $];
    }

    if (!exists(fastaPath)) {
        stderr.writeln("Error: can't find ", fastaPath);
        return 0;
    }

    string bakPath = fastaPath ~ ".bak";

    // 1) create backup if not exists
    if (!exists(bakPath)) {
        copy(fastaPath, bakPath);
        writeln("Backup file created: ", bakPath);
    } else {
	writeln("Backup file already exists: ", bakPath);
    } 

    // 2) read from backup and filter
    auto recs = readFasta(bakPath);

    FastaRecord[] kept;
    size_t removed = 0;

    foreach (r; recs) {
        if (shouldRemove(r.header, taxa)) {
            removed++;
        } else {
            kept ~= r;
        }
    }

    // 3) write back
    if (removed > 0) {
        writeFasta(fastaPath, kept);
        writeln("[", baseName(fastaPath), "] deleted ", removed, " sequences (taxa=", taxa, ")");
    } else {
        writeln("[", baseName(fastaPath), "] mismatch (taxa=", taxa, ")");
    }

    return 0;
}
