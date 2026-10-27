#!/usr/bin/env rdmd

import std.stdio;
import std.file;
import std.process;
import std.algorithm;
import std.conv;
import std.array;
import std.path;
import std.parallelism;
import std.regex;
import std.string;
import std.csv;
import core.atomic : atomicLoad, atomicOp, atomicStore;

shared size_t stageCommandCount;
shared size_t stageCommandFailures;
string activeCommandStage;

string checkpointDir() {
    return buildPath(std.file.getcwd(), ".rgbepp", "checkpoints");
}

void beginCommandBatch(string stage) {
    activeCommandStage = stage;
    atomicStore(stageCommandCount, 0);
    atomicStore(stageCommandFailures, 0);
    auto dir = checkpointDir();
    createDir(buildPath(std.file.getcwd(), ".rgbepp"));
    createDir(dir);
    foreach (entry; dirEntries(dir, stage ~ ".cmdfail.*", SpanMode.shallow)) {
        std.file.remove(entry.name);
    }
    auto done = buildPath(dir, stage ~ ".done");
    auto warning = buildPath(dir, stage ~ ".warning");
    if (done.exists) std.file.remove(done);
    if (warning.exists) std.file.remove(warning);
    std.file.write(buildPath(dir, stage ~ ".running"), stage ~ "\n");
}

void recordCommandResult(bool succeeded, string description) {
    atomicOp!"+="(stageCommandCount, cast(size_t) 1);
    if (!succeeded) {
        auto failureNumber = atomicOp!"+="(stageCommandFailures, cast(size_t) 1);
        auto failureFile = buildPath(checkpointDir(), activeCommandStage ~ ".cmdfail." ~ failureNumber.to!string);
        std.file.write(failureFile, description ~ "\n");
        stderr.writeln("Warning: command failed during ", activeCommandStage, ": ", description);
    }
}

string commandFailureReport(string checkpointPath, string stage, size_t failures, size_t attempts) {
    string report = to!string(failures) ~ " of " ~ to!string(attempts) ~ " command(s) failed\n";
    auto failureFiles = dirEntries(checkpointPath, stage ~ ".cmdfail.*", SpanMode.shallow)
        .map!(entry => entry.name)
        .array;
    sort(failureFiles);
    foreach (failureFile; failureFiles) {
        report ~= "Failed command: " ~ readText(failureFile).strip ~ "\n";
    }
    return report;
}

void clearCommandFailureFiles(string checkpointPath, string stage) {
    foreach (entry; dirEntries(checkpointPath, stage ~ ".cmdfail.*", SpanMode.shallow)) {
        std.file.remove(entry.name);
    }
}

void finishCommandBatch() {
    auto dir = checkpointDir();
    auto running = buildPath(dir, activeCommandStage ~ ".running");
    auto failed = buildPath(dir, activeCommandStage ~ ".failed");
    auto done = buildPath(dir, activeCommandStage ~ ".done");
    auto warning = buildPath(dir, activeCommandStage ~ ".warning");
    auto attempts = atomicLoad(stageCommandCount);
    auto failures = atomicLoad(stageCommandFailures);
    if (running.exists) std.file.remove(running);
    if (failures != 0) {
        if (failures == attempts) {
            std.file.write(failed, commandFailureReport(dir, activeCommandStage, failures, attempts));
            clearCommandFailureFiles(dir, activeCommandStage);
            throw new Exception(activeCommandStage ~ " failed: all " ~ to!string(attempts) ~ " command(s) failed.");
        }
        if (failed.exists) std.file.remove(failed);
        std.file.write(warning, commandFailureReport(dir, activeCommandStage, failures, attempts));
        std.file.write(done, "Completed with partial command failures\n");
        clearCommandFailureFiles(dir, activeCommandStage);
        stderr.writeln("Warning: ", activeCommandStage, " completed with ", failures,
            " of ", attempts, " command(s) failing; checkpoint marked complete.");
    } else {
        if (failed.exists) std.file.remove(failed);
        std.file.write(done, "Completed\n");
    }
    activeCommandStage = "";
}

void show_help(string pkgver) {
    writeln("\t\t\t\t\t\033[0;47;31mR\033[0m\033[0;47;92mG\033[0m\033[0;47;94mB\033[0m\033[0;47m \033[0m\033[0;47;33mE\033[0m\033[0;47;94mP\033[0m\033[0;47;33mP\033[0m
\t\t\tReference Genome based Exon Phylogeny Pipeline
	    Version: ", pkgver, "
	    License: GPL-2.0-only
	    Author: Guoyi Zhang
	    General options:
	    --codon\t\tOnly use the codon region (optional)
	    -c\t--config\tconfig file for software path (optional; takes precedence over tool path options)
	    -f\t--functions\tfunctions type (optional): all clean assembly map postmap
	      \t           \tvarcall consen codon align ortholog trim concat
	    -g\t--genes\t\tgene file path (optional, if -r is specified)
	    -h\t--help\t\tshow this information
	    -l\t--list\t\tlist file path
	    --log\t\tlog file path (optional, default rgbepp.log)
	    -m\t--memory\tmemory settings (optional, default 16 GB)
	    -r\t--reference\treference genome path
	    --resume\t\tresume using stage checkpoints (optional)
	    -t\t--threads\tthreads setting (optional, default 8 threads)
	    Tool paths (optional; command names are resolved through PATH):
	    --bcftools\t\tBcftools path
	    --bowtie2\t\tBowtie2 path
	    --concataln\t\tConcataln path
	    --delstop\t\tDelstop path
	    --deltaxa\t\tDeltaxa path
	    --diamond\t\tDiamond path
	    --exonerate\tExonerate path
	    --fastp\t\tFastp path
	    --macse\t\tMacse jarfile path
	    --samtools\tSamtools path
	    --sortdiamond\tSortDiamond path
	    --spades\t\tSpades python path
	    --trimal\t\tTrimal path
	    Directories (optional; relative paths use the run directory):
	    --raw_dir\t\tInput FASTQ directory (default: 00_raw)
	    --fastp_dir\tOutput directory (default: 01_fastp)
	    --spades_dir\tOutput directory (default: 02_spades)
	    --bowtie2_dir\tOutput directory (default: 03_bowtie2)
	    --bam_dir\t\tOutput directory (default: 04_bam)
	    --vcf_dir\t\tOutput directory (default: 05_vcf)
	    --consen_dir\tOutput directory (default: 06_consen)
	    --macse_dir\tOutput directory (default: 07_macse)
	    --trimal_dir\tOutput directory (default: 08_trimal)
	    When a directory is set both on the command line and in config, the config value wins with a warning.
	    for example: ./RGBEPP -f all -l list -t 8 -r reference.fasta \n");
}

bool testJava() {
    bool pass = true;
    auto result = execute(["java", "-version"]);
    if (result.status != 0) {
	pass = false;
	writeln("Error: Java is not found");
    } 
    return pass;
}

bool testFiles(string[] filePaths) {
    bool pass = true;
    foreach(filePath; filePaths){
	if (filePath.empty) {
	    writeln("Error: a required tool path was not provided.");
	    pass = false;
	    continue;
	}
	// Explicit paths are checked directly; command names are resolved through PATH.
	bool found = filePath.canFind('/') || filePath.canFind('\\')
	    ? exists(filePath)
	    : environment.get("PATH", "").split(pathSeparator)
		.any!(directory => exists(buildPath(directory, filePath)));
	if (!found) {
	    writeln("Error: required tool was not found: " ~ filePath);
	    pass = false;
        } 
    }
    return pass;
} 

bool testInputFiles(string[] filePaths) {
    bool pass = true;
    foreach (filePath; filePaths) {
        if (filePath.empty || !filePath.isFile) {
            writeln("Error: required input file was not found: ", filePath);
            pass = false;
        }
    }
    return pass;
}

bool testString(string input) {
    return input.length != 0;
}

bool testStringArray(string[] input) {
    return !input.empty;
}

bool filesExist(string[] paths) {
    if (paths.empty) return false;
    foreach (path; paths) {
        if (!path.isFile) return false;
    }
    return true;
}

bool skipCompletedStage(string stage, string[] expectedOutputs, bool resume) {
    auto dir = checkpointDir();
    bool hasIncompleteCheckpoint = buildPath(dir, stage ~ ".failed").exists
        || buildPath(dir, stage ~ ".running").exists;
    bool hasCheckpoint = buildPath(dir, stage ~ ".done").exists;
    if (resume && !hasIncompleteCheckpoint && (hasCheckpoint || filesExist(expectedOutputs))) {
        writeln(stage, "::Skipped (completed checkpoint/outputs found; --resume)");
        return true;
    }
    return false;
}

void createDir(string path) {
    if (!exists(path)) {
        mkdir(path);
    }
}

void moveDir (string oldPath, string newPath) {
    try {
        rename(oldPath, newPath);
        writeln("Directory renamed successfully.");
    } catch (Exception e) {
        writeln("Error renaming directory: ", e.msg);
    }
}

void executeCommand(string[] cmd) {
    auto process = spawnProcess(cmd);

    bool succeeded = wait(process) == 0;
    recordCommandResult(succeeded, cmd.join(" "));
}

void executeCommandToFile(string[] cmd, string outputFile) {
    // Create a pipe for the command's output
    auto pipe = pipe();

    // Spawn the process
    auto pid = spawnProcess(cmd, stdin, pipe.writeEnd);

    // Close the write end of the pipe to signal EOF
    pipe.writeEnd.close();

    // Read the output from the pipe
    auto output = cast(string) pipe.readEnd.byChunk(4096).joiner.array;

    // Wait for the process to finish
    bool succeeded = wait(pid) == 0;
    recordCommandResult(succeeded, cmd.join(" "));

    // Write the output to the specified file
    std.file.write(outputFile, cast(ubyte[])output);
}

void initializeLog(string logPath) {
    auto logFile = File(logPath, "a");
    stdout.flush();
    stderr.flush();
    stdout = logFile;
    stderr = logFile;
    writeln("\n=== RGBEPP run started ===");
}


void executeCommandPipe(string[][] cmds) {

    Pid[] pids;
    bool failed;
    scope(exit) {
	    foreach (pid; pids) {
		    if (wait(pid) != 0) failed = true;
	    }
        recordCommandResult(!failed, cmds.map!(cmd => cmd.join(" ")).join(" | "));
    }

    // pipe init
    auto temp_pipe = pipe();
    // process first
    pids ~= spawnProcess(cmds[0], stdin, temp_pipe.writeEnd);

    // process cmd2 ~ cmdN-1
    for (int i = 1; i < cmds.length - 1; i++) {
        auto new_pipe = pipe(); // create next pipe
        pids ~= spawnProcess(cmds[i], temp_pipe.readEnd, new_pipe.writeEnd);
        temp_pipe = new_pipe; // update the pipe
    }

    // process final, output to stdout
    pids ~= spawnProcess(cmds[$-1], temp_pipe.readEnd, stdout);
}

string[] readArrFromFile(string filename) {
    string[] arr;

    try {
        arr = filename.readText().splitter.array;
    } catch (FileException ex) {
        writeln("Error reading file: ", ex.msg);
    } catch (Exception ex) {
        writeln("Exception: ", ex.msg);
    }

    return arr;
}

string getBaseName(string ARG_R){
    string ARG_R_extension = extension(ARG_R); // get extension
    string baseNameRef = baseName(ARG_R, ARG_R_extension); //rm dir and extension
    return baseNameRef;
}

string[] getRef(string ARG_R, string DirMap){
    string baseNameRef = getBaseName(ARG_R);
    string ARG_R_index = buildPath(DirMap, "index",  baseNameRef); // bt2_index_base
    string ARG_R_refer = ARG_R_index ~ ".fasta"; //reference_in fasta file
    string[] Refs = [ARG_R_index, ARG_R_refer];
    return Refs;
}

string[] getARG_G(string ARG_R){
    string[] ARG_G;
    // if ARG_G is empty
    if (ARG_G.length == 0) {
        auto file = File(ARG_R, "r");
        ARG_G = file.byLine
                .filter!(line => line.startsWith(">")) // flitering
                .map!(line => line[1..$].idup) // convert to word
                .array;
    }
    return ARG_G;
}

string getValueFromConfig(string file, string key) {
    string content = readText(file);
    string value;
    auto regex = regex("^" ~ key ~ r"\s*=\s*(.+)");

    foreach (line; content.splitter("\n")) {
        if (auto match = matchFirst(line, regex)) {
            value = match.captures[1].strip;;
            break;
        }
    }

    return value;
}

string getConfiguredValue(string file, string key, string currentValue, string[string] commandLineOptions) {
    string configValue = getValueFromConfig(file, key);
    if (!configValue.empty) {
        if (auto commandLineValue = key in commandLineOptions) {
            if (*commandLineValue != configValue) {
                writeln("Warning: --", key, " (", *commandLineValue,
                    ") conflicts with config value (", configValue,
                    "); using the config value.");
            }
        }
        return configValue;
    }
    return currentValue;
}


string[] getValueArray(string value) {
    return value.length
        ? value.splitter(" ").filter!(a => a.length).array
        : [];
}

void processQcTrim(string[] ARG_L, int ARG_T, string DirRaw, string DirQcTrim, string PathFastp, string[] fastpParas) {
    // Prepare directory
    createDir(DirQcTrim);
    writeln("QcTrimming::Start");
    foreach (string file; ARG_L) {
        string baseName = getBaseName(file);
        string inputFileR1 = buildPath(DirRaw, baseName ~ "_R1.fastq.gz");
        string inputFileR2 = buildPath(DirRaw, baseName ~ "_R2.fastq.gz");
        string outputFileR1 = buildPath(DirQcTrim, baseName ~ "_R1.fastq.gz");
        string outputFileR2 = buildPath(DirQcTrim, baseName ~ "_R2.fastq.gz");
        string jsonFile = buildPath(DirQcTrim, baseName ~ ".json");
        string htmlFile = buildPath(DirQcTrim, baseName ~ ".html");

        // Perform quality control and trimming using external program `fastp`
        string[] cmdQcTrim = [PathFastp, "-i", inputFileR1, "-I", inputFileR2,
                        "-o", outputFileR1, "-O", outputFileR2,
                        "-j", jsonFile, "-h", htmlFile,
                        "-w", ARG_T.to!string] ~ fastpParas;
        executeCommand(cmdQcTrim);
    }
    writeln("QcTrimming::End");
}

void processAssembly(string[] ARG_L, int ARG_M, int ARG_T, string DirQcTrim, string DirAssembly, string PathSpades, string[] spadesParas){
    writeln("Assembly::Start");
    createDir(DirAssembly);
    foreach (string file; ARG_L) {
       string baseName = getBaseName(file);
       string DirAss = buildPath(DirAssembly, baseName);
       createDir(DirAss);
       string inputFileR1 = buildPath(DirQcTrim, baseName ~ "_R1.fastq.gz");
       string inputFileR2 = buildPath(DirQcTrim, baseName ~ "_R2.fastq.gz");
       string[] cmdAssembly = [PathSpades, "--pe1-1", inputFileR1, "--pe1-2", inputFileR2, "-t", ARG_T.to!string, "-m", ARG_M.to!string] ~ spadesParas ~ ["-o", DirAss];
    	executeCommand(cmdAssembly);
    }
    writeln("Assembly::End");
}

void processAssemMv(string[] ARG_L,string DirAssembly){
    // Prepare
    string DirAssemblySca = buildPath(DirAssembly, "scaffolds");
    string DirAssemblyCont = buildPath(DirAssembly, "contigs");
    writeln("Assembly_Move::Start");
    createDir(DirAssemblySca);
    createDir(DirAssemblyCont);
    foreach (string file; ARG_L ){
        string baseName = getBaseName(file);
	string DirAssemblyInd = buildPath(DirAssembly, baseName);
	string inputSca = buildPath(DirAssemblyInd, "scaffolds.fasta");
	string inputCont = buildPath(DirAssemblyInd, "contigs.fasta");
	string outputSca = buildPath(DirAssemblySca,  baseName ~ ".fasta");
	string outputCont = buildPath(DirAssemblyCont, baseName ~ ".fasta");
	if (!exists(inputSca)) {
            writeln("File not found: ", inputSca);
            continue;
        } else {
	    copy(inputSca, outputSca);
	}
        if (!exists(inputCont)) {
            writeln("File not found: ", inputCont);
            continue;
        } else {
	    copy(inputCont, outputCont);
	}
    }
    writeln("Assembly_Move::End");
}

void processMappingDenovo(string[] ARG_L, string ARG_R, int ARG_T, string DirQcTrim, string DirAssembly, string DirMap, string PathBowtie2, string PathDiamond, string PathSamtools, string PathSortDiamond, string[] diamondMakedbParas, string[] diamondBlastxParas, string[] bowtie2BuildParas,  string[] bowtie2Paras, string[] samtoolsViewParas){
    // Prepare directory
    writeln("Mapping::Start");
    createDir(DirMap);
    createDir(buildPath(DirMap, "index"));
    string DirAssemblySca = buildPath(DirAssembly, "scaffolds");
    string DirAssemblyFas = buildPath(DirAssembly, "fasta");
    createDir(DirAssemblyFas);
    
    string ARG_R_Base = getBaseName(ARG_R);
    string ARG_R_Ref = buildPath(DirAssemblyFas, ARG_R_Base ~ ".fasta");
    copy(ARG_R, ARG_R_Ref);
    string [] cmdDmMakeDB = [ PathDiamond, "makedb", "--db", "Reference", "--in", ARG_R_Ref] ~ diamondMakedbParas;
    executeCommand(cmdDmMakeDB);
    string ReferDmnd = buildPath(DirAssemblyFas, "Reference.dmnd");
    string PathBowtie2_build = PathBowtie2 ~ "-build";

    foreach (string file; ARG_L) {
	string baseName = getBaseName(file);
	string inputM8 = buildPath(DirAssemblySca, baseName ~ ".m8");
	string inputFasta = buildPath(DirAssemblySca,  baseName ~ ".fasta");
	string outputSort = buildPath(DirAssemblyFas, baseName ~ ".fasta");
	string outputIndex = buildPath(DirAssemblyFas, baseName);
        string inputFileR1 = buildPath(DirQcTrim, baseName ~ "_R1.fastq.gz");
        string inputFileR2 = buildPath(DirQcTrim, baseName ~ "_R2.fastq.gz");
	string outputBam = buildPath(DirMap, baseName ~ ".bam");

	string[] cmdDiamond = [PathDiamond, "blastx", "-d", "Reference.dmnd", "-q", inputFasta, "-o", inputM8] ~ diamondBlastxParas ~ ["--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "gaps", "ppos", "qframe", "qseq"];
	string[] cmdSortDiamond = [PathSortDiamond, inputM8, outputSort];
    	string[] cmdBuildDB = [PathBowtie2_build, "--threads", ARG_T.to!string, outputSort, outputIndex];
        string[] cmdMap = [PathBowtie2, "-x", outputIndex, "-1", inputFileR1, "-2", inputFileR2, "-p", ARG_T.to!string];
        string[] cmdSam2Bam = [PathSamtools, "view", "-bS", "-@", ARG_T.to!string, "-o", outputBam] ~ samtoolsViewParas;
	executeCommand(cmdDiamond);
	executeCommand(cmdSortDiamond);
    	executeCommand(cmdBuildDB);
	executeCommandPipe([cmdMap, cmdSam2Bam]);
    }

    writeln("Mapping::End");
}

void processPostMap(string[] ARG_L, int ARG_T, string DirMap, string DirBam, string PathSamtools, string[] samtoolsFixmateParas, string[] samtoolsSortParas, string[] samtoolsMarkdupParas, string[] samtoolsIndexParas) {

    createDir(DirBam);
    writeln("PostMapping::Start");

    foreach (string file; ARG_L) {
        string baseName = getBaseName(file);
        string inputBam = buildPath(DirMap, baseName ~ ".bam");
        string outputBam = buildPath(DirBam, baseName ~ ".bam");

        // Convert SAM to BAM, sort and remove duplicates using Samtools
	string[] cmdFixmate = [PathSamtools, "fixmate", "-@", ARG_T.to!string, "-m", inputBam, ] ~ samtoolsFixmateParas ~ ["-"];
        string[] cmdSort = [PathSamtools, "sort", "-@", ARG_T.to!string] ~ samtoolsSortParas ~ ["-"];
        string[] cmdMarkdup = [PathSamtools, "markdup", "-@", ARG_T.to!string] ~ samtoolsMarkdupParas ~ ["-", outputBam];
	executeCommandPipe([cmdFixmate, cmdSort, cmdMarkdup]);

        string [] cmdIndexBam = [PathSamtools, "index", "-@", ARG_T.to!string] ~ samtoolsIndexParas ~ [outputBam];
        executeCommand(cmdIndexBam);
    }

    writeln("PostMapping::End");
}

void processVarCallDenovo(string[] ARG_L, int ARG_T, string DirAssembly, string DirMap, string DirBam, string DirVcf, string PathBcftools, string[] bcftoolsMpileupParas, string[] bcftoolsCallParas, string[] bcftoolsNormParas, string[] bcftoolsFilterParas) {
    writeln("VarCalling::Start");

    string DirAssemblyFas = buildPath(DirAssembly, "fasta");
    createDir(DirVcf);

    foreach (string file; parallel(ARG_L, 1)) {
        string baseName = getBaseName(file);
        string inputBam = buildPath(DirBam, baseName ~ ".bam");
        string outputVcf = buildPath(DirVcf, baseName ~ ".vcf.gz");
	string referFasta = buildPath(DirAssemblyFas, baseName ~ ".fasta");
        // Variant calling using bcftools
        string[] cmdPileup = [PathBcftools, "mpileup", "-Oz", "--threads", ARG_T.to!string, "-f", referFasta] ~ bcftoolsMpileupParas ~ [inputBam];
	string[] cmdVarCall = [PathBcftools, "call"] ~ bcftoolsCallParas ~ ["-Oz", "--threads", ARG_T.to!string];
	string[] cmdNorm = [PathBcftools, "norm", "--threads", ARG_T.to!string, "-f", referFasta] ~ bcftoolsNormParas ~ ["-Oz"];
	string[] cmdFilter = [PathBcftools, "filter", "--threads", ARG_T.to!string, "-Oz"] ~ bcftoolsFilterParas ~ ["-o", outputVcf];
        executeCommandPipe([cmdPileup, cmdVarCall, cmdNorm, cmdFilter]); 
    }

    writeln("VarCalling::End");

}


void processConDenovo(string[] ARG_G, string[] ARG_L, int ARG_T, string DirAssembly,  string DirVcf, string DirConsensus, string PathBcftools, string[] bcftoolsIndexParas, string[] bcftoolsConParas) {
    createDir(DirConsensus);

    string DirConTaxa = buildPath(DirConsensus, "taxa");
    string DirAssemblyFas = buildPath(DirAssembly, "fasta");
    createDir(DirConTaxa);

    writeln("Consensus::Start");
    // Extract fasta from vcf file
    foreach (string file; ARG_L) {
        string baseName = getBaseName(file);
        string inputVcf = buildPath(DirVcf, baseName ~ ".vcf.gz");
        string outputFasta = buildPath(DirConTaxa, baseName ~ ".fasta");
	string referFasta = buildPath(DirAssemblyFas, baseName ~ ".fasta");
	// index vcf.gz
	string[] cmdIndexVcf = [PathBcftools, "index"] ~ bcftoolsIndexParas ~ [inputVcf];
	executeCommand(cmdIndexVcf);
	
        // Generate consensus sequences using bcftools
        string[] cmdCon = [PathBcftools, "consensus", "-f", referFasta] ~ bcftoolsConParas ~ [inputVcf, "-o", outputFasta];
	executeCommand(cmdCon);
    }
    // Recombine the sequences based on genes
    writeln("Consensus::End");
}

void processCombFasta(string[] ARG_G, string[] ARG_L, string DirConsensus, string strTaxa, string strGene) {

    string DirConTaxa = buildPath(DirConsensus, strTaxa);
    string DirConGene = buildPath(DirConsensus, strGene);    
    createDir(DirConGene);

    // create a dictory
    string[string] geneSequences;

    writeln("ConvertFasta::Start");
    // read first
    foreach (file; ARG_L) {
        string inputFile = buildPath(DirConTaxa, file ~ ".fasta");
        if (!exists(inputFile)) {
            writeln("File not found: ", inputFile);
            continue;
        }
        string content = cast(string) readText(inputFile);
        bool inSequence = false;
        string currentGene;

        foreach (line; content.splitter("\n")) {
            if (line.empty) continue;

            if (line[0] == '>') {
                string header = line[1 .. $];
                if (ARG_G.canFind(header)) {
                    currentGene = header;
                    geneSequences[currentGene] ~= ">" ~ file ~ "\n";
                    inSequence = true;
                } else {
                    inSequence = false;
                }
            } else if (inSequence) {
                geneSequences[currentGene] ~= line ~ "\n";
            }
        }
    }
    // write different files
    foreach (gene; ARG_G) {
        string outputFile = buildPath(DirConGene, gene ~ ".fasta");
        File output = File(outputFile, "w");
        if (gene in geneSequences) {
            output.write(geneSequences[gene]);
        }
    }
    writeln("ConvertFasta::End");
}

void splitFasta(string inputFasta, string DirOut) {
    File infile;
    try {
        infile = File(inputFasta, "r");
    } catch (FileException e) {
        stderr.writeln("Error: Unable to open input file ", inputFasta);
        return;
    }

    string line;
    string seqName;
    File outfile;
    bool in_sequence = false;

    foreach (lineContent; infile.byLine()) {
        if (lineContent.empty) continue;

        if (lineContent[0] == '>') {
            // New sequence header
            if (in_sequence) {  // if found new sequence, close
                outfile.close();  // previous output file
            }
            seqName = cast(string)lineContent[1 .. $];  // Remove '>'
            string outputFile = buildPath(DirOut, seqName ~ ".fasta");  // suitable to many os
            outfile = File(outputFile, "w");
            outfile.writeln(">", getBaseName(inputFasta));
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

void processCodon(string[] ARG_G, string ARG_R, string DirConsensus, string PathExonerate){

    string DirConGene = buildPath(DirConsensus , "gene"); 

    string ARG_R_Base = getBaseName(ARG_R);
    string ARG_R_Ref = buildPath(DirConsensus, ARG_R_Base ~ ".fasta");
    copy(ARG_R, ARG_R_Ref);
    splitFasta(ARG_R_Ref, DirConsensus);

    if (!exists(DirConGene ~ "_bak")) {
        moveDir(DirConGene, DirConGene ~ "_bak");
    }
    if (!exists(DirConGene)) {
        createDir(DirConGene);
    }
    
    writeln("GetCodon::Start");

    foreach (gene; ARG_G) {
	string inputFile = buildPath(DirConGene ~ "_bak", gene ~ ".fasta");
        string outputFile = buildPath(DirConGene, gene ~ ".fasta");
	string referFile = buildPath(DirConsensus, gene ~ ".fasta");

        if (!exists(inputFile)) {
            writeln("File not found: ", inputFile);
            continue;
        } else {
            string[] cmdExonerate = [PathExonerate, inputFile, referFile, "--showalignment", "no", "--showvulgar", "no", "--showtargetgff", "no", "--ryo", ">%qi\n%qcs\n", "--verbose", "0"];
	    executeCommandToFile(cmdExonerate, outputFile);
	}
	std.file.remove(referFile);
    }
    
    rmdirRecurse(DirConGene ~ "_bak");

    writeln("GetCodon::End");
}

void processAlign(string[] ARG_G, string DirConsensus, string DirAlign, string PathMacse, string[] macseParas){

    string DirConGene = buildPath(DirConsensus, "gene");
    string DirAlignAA = buildPath(DirAlign, "AA");
    string DirAlignNT = buildPath(DirAlign, "NT");

    writeln("Align::Start");
    createDir(DirAlign);
    createDir(DirAlignAA);
    createDir(DirAlignNT);
    foreach (gene; parallel(ARG_G, 1)) {
    	string inputFasta = buildPath(DirConGene, gene ~ ".fasta");
    	string outAA = buildPath(DirAlignAA, gene ~ ".fasta");
    	string outNT = buildPath(DirAlignNT, gene ~ ".fasta");
	if (!exists(inputFasta)) {
            writeln("File not found: ", inputFasta);
            continue;
        } else{
    	    string[] cmdAlign = ["java", "-jar", PathMacse, "-prog", "alignSequences", "-seq" , inputFasta, "-out_AA", outAA, "-out_NT", outNT ] ~ macseParas;
    	    executeCommand(cmdAlign);
	}
    }
    writeln("Align::End");

}

void processReplaceQue(string DirAlign, string[] ARG_L, string strTaxa)
{
    // replace questio with minus of fasta files in specific folder
    foreach (taxa; ARG_L)
    {
	string fnameBase = taxa ~ ".fasta";
	string fname = buildPath(DirAlign, strTaxa, fnameBase);
        if (!fname.isFile) {
	    writeln("File not found: ", fname);
	    continue;
	} else {
            string text = readText(fname);
            string textReplaced = text.replace("!", "-");

            // Once change, write it
            if (textReplaced != text)
                std.file.write(fname, textReplaced);
	}
    }
}

void processCompareTSV(string[] ARG_L, string DirAlign, string strTaxa, string strParalog)
{
    // compare the tsv files and find the paralogous genes, then write to a csv file
    File outfile;
    outfile = File(strParalog, "w");

    auto writeLine = (string s) {
        outfile.writeln(s);
    };
    foreach (taxa; ARG_L)
    {
	string fnameBase = taxa ~ ".tsv";
        string fname = buildPath(DirAlign, strTaxa, fnameBase);
        string displayName;
        File f = File(fname, "r");
        foreach (line; f.byLine)
        {
            string s = line.idup.chomp;
            if (s.length == 0) continue;

            string[] cols = s.split('\t');
            if (cols.length >= 2)
            {
                string col1 = cols[0];
                string col2 = cols[1];
                if (col1 != col2)
                    writeLine(taxa ~ "," ~ col1 ~ "," ~ col2);
            }
        }
    }
    outfile.close();
}

void processMismatch(string strParalog, string PathDeltaxa, string DirAlign)
{
    // Delete the paralogs

    // store information
    alias TaxaSet = bool[string];
    TaxaSet[string] geneToTaxaSet;

    string text = readText(strParalog);

    foreach (rowRange; csvReader!string(text))
    {
        string[] row = rowRange.array;
        if (row.length < 2) continue;

        string taxa = row[0].strip;
        string gene = row[1].strip;

        if (taxa.length == 0 || gene.length == 0) continue;

        geneToTaxaSet[gene][taxa] = true; // deduplication
    }

    // run for everygene
    immutable string[2] modes = ["NT", "AA"];
    
    foreach (gene, taxaSet; geneToTaxaSet)
    {
        string[] taxaList = taxaSet.byKey.array;
        sort(taxaList);
    
        foreach (mode; modes)
        {
            string fastaPath = buildPath(DirAlign, mode, gene ~ ".fasta");
    
            // PathDeltaxa fastaPath taxa1 taxa2 ...
            auto cmd = appender!(string[])();
            cmd.put(PathDeltaxa);
            cmd.put(fastaPath);
            foreach (t; taxaList) cmd.put(t);
    
            executeCommand(cmd.data);
        }
    }
}

void processOrtholog(string[] ARG_L, string ARG_R, string DirAlign, string PathDeltaxa, string PathDiamond,  string strParalog){
    string DirAlignTaxa = buildPath(DirAlign, "taxa"); 

    writeln("Preparing reference::Start");
    processReplaceQue(DirAlign, ARG_L, "taxa");
    writeln("Preparing reference::End");


    writeln("Ortholog::Start");
    foreach (taxa; parallel(ARG_L, 1)) {
    	string inputFasta = buildPath(DirAlignTaxa, taxa ~ ".fasta");
	string outputTSV = buildPath(DirAlignTaxa, taxa ~ ".tsv");
	string ARG_R_Ref = buildPath(DirAlignTaxa, taxa);

	if (!exists(inputFasta)) {
            writeln("File not found: ", inputFasta);
            continue;
        } else{
    	    string[] cmdDmMakeDB2 = [PathDiamond, "makedb", "--db", ARG_R_Ref, "--in", inputFasta] ;
	    string[] cmdDmBlastp = [PathDiamond, "blastp", "-q", ARG_R, "--max-target-seqs", "1", "-d", ARG_R_Ref, "-o", outputTSV, "-f", "6", "sseqid", "qseqid"];
    	    executeCommand(cmdDmMakeDB2);
	    executeCommand(cmdDmBlastp);
	}
    }
    writeln("Ortholog::End");

    writeln("Getting paralog::Start");
    processCompareTSV( ARG_L, DirAlign, "taxa", strParalog);
    writeln("Getting paralog::End");

    writeln("Removing paralog::End");
    processMismatch(strParalog, PathDeltaxa, DirAlign);
    writeln("Removing paralog::End");
}

void processTrimming(string[] ARG_G, string DirAlign, string DirTrim, string PathDelstop, string PathTrimal, string[] trimalParas){
    writeln("Trimming::Start");

    string DirAA = buildPath(DirAlign, "AA");
    string DirNT = buildPath(DirAlign, "NT");
    string DirAA_out = buildPath(DirAlign, "AA_out");
    string DirNT_out = buildPath(DirAlign, "NT_out");

    createDir(DirAA_out);
    createDir(DirNT_out);
    
    // copy file firstly
    foreach (gene; parallel(ARG_G,1)){
   	string inputFastaAA = buildPath(DirAA, gene ~ ".fasta");
	string outputFastaAA = buildPath(DirAA_out, gene ~ ".fasta");
   	string inputFastaNT = buildPath(DirNT, gene ~ ".fasta");
	string outputFastaNT = buildPath(DirNT_out, gene ~ ".fasta");
  
       	if (!exists(inputFastaNT)) {
            writeln("File not found: ", inputFastaNT);
            continue;
        } else{	
            copy(inputFastaNT, outputFastaNT);
	}
	if (!exists(inputFastaAA)) {
            writeln("File not found: ", inputFastaAA);
            continue;
        } else{
            copy(inputFastaAA, outputFastaAA);
	}  
        // del stop codon
        string[] cmdDelStop = [PathDelstop, outputFastaAA, outputFastaNT, "--delete"];
        executeCommand(cmdDelStop);	
    }

    string DirTrimNT = buildPath(DirTrim, "NT");
    createDir(DirTrim);
    createDir(DirTrimNT);
    foreach (gene; parallel(ARG_G,1)){
        string inputFastaAA = buildPath(DirAA_out, gene ~ ".fasta");
	string inputBackTransNT = buildPath(DirNT_out, gene ~ ".fasta");
	string outputFastaNT = buildPath(DirTrimNT, gene ~ ".fasta");
	if (exists(inputFastaAA) && exists(inputBackTransNT)) {
            string[] cmdTrim = [PathTrimal, "-in", inputFastaAA, "-backtrans", inputBackTransNT, "-out", outputFastaNT] ~ trimalParas;
            executeCommand(cmdTrim);
        } else {
            writeln("Skipping gene: ", gene, " as files are missing.");
        }
    }
    writeln("Trimming::End");

}

void processConcat(string[] ARG_G, string DirTrim, string PathConcataln){
    writeln("Concatenate::Start");
    string PathOutput = buildPath(DirTrim, "concat.fasta");
    string[] cmdConcat = [PathConcataln, PathOutput];

     foreach (a; ARG_G)
     {
         auto fname = buildPath(DirTrim, "NT", a ~ ".fasta");
	 if (exists(fname)) {
             cmdConcat ~= fname;
	 } else {
            writeln("Skipping gene: ", a, " as files are missing.");
	 }
     }
    executeCommand(cmdConcat);
    writeln("Concatenate::End");
}


void main(string[] args) {
    string pkgver = "0.0.5";

    string DirHome = std.file.getcwd();
    string DirRaw = buildPath(DirHome, "00_raw");
    string DirQcTrim = buildPath(DirHome, "01_fastp");
    string DirAssembly = buildPath(DirHome, "02_spades");
    string DirMap = buildPath(DirHome, "03_bowtie2");
    string DirBam = buildPath(DirHome, "04_bam");
    string DirVcf = buildPath(DirHome, "05_vcf");
    string DirConsensus = buildPath(DirHome, "06_consen");
    string DirAlign = buildPath(DirHome, "07_macse"); 
    string DirTrim = buildPath(DirHome, "08_trimal"); 

    string PathFastp = "fastp";
    string PathSpades = "spades.py";
    string PathDiamond = "diamond";
    string PathSortDiamond = "sortdiamond";
    string PathBowtie2 = "bowtie2";
    string PathSamtools = "samtools";
    string PathBcftools = "bcftools";
    string PathExonerate = "exonerate";
    string PathMacse = "/usr/share/java/macse.jar";
    string PathDelstop = "delstop";
    string PathDeltaxa = "deltaxa";
    string PathTrimal = "trimal";
    string PathConcataln = "concataln";

    string[] fastpParas = [];
    string[] spadesParas = ["--careful", "--phred-offset", "33"];
    string[] diamondMakedbParas = [];
    string[] diamondBlastxParas = ["--ultra-sensitive"];
    string[] bowtie2BuildParas = [];
    string[] bowtie2Paras = [];
    string[] samtoolsViewParas = [];
    string[] samtoolsFixmateParas = [];
    string[] samtoolsSortParas = [];
    string[] samtoolsMarkdupParas = [];
    string[] samtoolsIndexParas = [];
    string[] bcftoolsMpileupParas = [];
    string[] bcftoolsCallParas = ["-mv"];
    string[] bcftoolsNormParas = ["--check-ref", "s"];
    string[] bcftoolsFilterParas = [];
    string[] bcftoolsIndexParas = [];
    string[] bcftoolsConParas = [];
    string[] exonerateParas = [];
    string[] macseParas = [];
    string[] trimalParas = ["-gt", "0.7"];

    int ARG_T = 8;
    int ARG_M = 16;
    string[] ARG_G;
    string[] ARG_L;
    string[] ARG_LIST_PATHS;
    string[] ARG_GENE_PATHS;
    string ARG_C;
    string ARG_F;
    string ARG_R;
    string ARG_LOG = "rgbepp.log";
    string[string] commandLineOptions;
    string strTaxa = "taxa";
    string strGene = "gene";
    string strParalog = "paralog.csv"; 
    bool enableCodon = false;
    bool enableResume = false;

    if (args.length > 1){
        foreach (int i; 0 .. cast(int)args.length) {
            switch (args[i]) {
                case "-c", "--config":
		    i++;
                    ARG_C = args[i];
                    break;
                case "-f", "--functions":
		    i++;
                    ARG_F = args[i];
                    break;
                case "-g", "--gene", "--genes":
		    i++;
		    ARG_GENE_PATHS ~= args[i];
                    ARG_G ~= readArrFromFile(args[i]);
		    break;
                case "-h", "--help":
                    show_help(pkgver);
                    return;
                case "-l", "--list":
		    i++;
		    ARG_LIST_PATHS ~= args[i];
                    ARG_L ~= readArrFromFile(args[i]);
                    break;
                case "--log":
		    i++;
                    ARG_LOG = args[i];
                    break;
                case "--resume":
                    enableResume = true;
                    break;
                case "-m", "--memory":
		    i++;
                    ARG_M = args[i].to!int;
                    break;
                case "-r", "--reference":
		    i++;
                    ARG_R = args[i];
                    break;
                case "-t", "--threads":
		    i++;
                    ARG_T = args[i].to!int;
                    break;
                case "--codon":
                    enableCodon = true;
                    break;
                case "--fastp":
		    i++;
                    PathFastp = args[i];
                    commandLineOptions["fastp"] = PathFastp;
                    break;
                case "--spades":
		    i++;
                    PathSpades = args[i];
                    commandLineOptions["spades"] = PathSpades;
                    break;
                case "--diamond":
		    i++;
                    PathDiamond = args[i];
                    commandLineOptions["diamond"] = PathDiamond;
                    break;
		case "--sortdiamond":
		    i++;
                    PathSortDiamond = args[i];
                    commandLineOptions["sortdiamond"] = PathSortDiamond;
                    break;
                case "--bowtie2":
		    i++;
                    PathBowtie2 = args[i];
                    commandLineOptions["bowtie2"] = PathBowtie2;
                    break;
                case "--samtools":
		    i++;
                    PathSamtools = args[i];
                    commandLineOptions["samtools"] = PathSamtools;
                    break;
                case "--bcftools":
		    i++;
                    PathBcftools = args[i];
                    commandLineOptions["bcftools"] = PathBcftools;
                    break;
                case "--exonerate":
		    i++;
                    PathExonerate = args[i];
                    commandLineOptions["exonerate"] = PathExonerate;
                    break;
                case "--macse":
		    i++;
                    PathMacse = args[i];
                    commandLineOptions["macse"] = PathMacse;
                    break;
                case "--delstop":
		    i++;
                    PathDelstop = args[i];
                    commandLineOptions["delstop"] = PathDelstop;
                    break;
                case "--deltaxa":
		    i++;
                    PathDeltaxa = args[i];
                    commandLineOptions["deltaxa"] = PathDeltaxa;
                    break;
                case "--trimal":
		    i++;
                    PathTrimal = args[i];
                    commandLineOptions["trimal"] = PathTrimal;
                    break;
                case "--concataln":
		    i++;
                    PathConcataln = args[i];
                    commandLineOptions["concataln"] = PathConcataln;
                    break;
                case "--raw_dir", "--fastp_dir", "--spades_dir", "--bowtie2_dir",
                    "--bam_dir", "--vcf_dir", "--consen_dir", "--macse_dir", "--trimal_dir":
                    i++;
                    string key = args[i - 1][2 .. $];
                    string value = args[i];
                    commandLineOptions[key] = value;
                    switch (key) {
                        case "raw_dir": DirRaw = value; break;
                        case "fastp_dir": DirQcTrim = value; break;
                        case "spades_dir": DirAssembly = value; break;
                        case "bowtie2_dir": DirMap = value; break;
                        case "bam_dir": DirBam = value; break;
                        case "vcf_dir": DirVcf = value; break;
                        case "consen_dir": DirConsensus = value; break;
                        case "macse_dir": DirAlign = value; break;
                        case "trimal_dir": DirTrim = value; break;
                        default: break;
                    }
                    break;
                default:
                    break;
            }
        }
    } else {
        show_help(pkgver);
        return;
    }

    initializeLog(ARG_LOG);
    writeln("Command: ", args.join(" "));

    // get gene from ARG_R reference fasta
    if (ARG_R.length != 0 && ARG_G.length == 0 ){
    	ARG_G = getARG_G(ARG_R);
    }
  
    if (ARG_F.length == 0 ){
    	ARG_F = "all";
    }

    // get pathXXX form config file
    if (ARG_C.length != 0){
    
        PathFastp = getConfiguredValue(ARG_C, "fastp", PathFastp, commandLineOptions);
	PathSpades = getConfiguredValue(ARG_C, "spades", PathSpades, commandLineOptions);
	PathDiamond = getConfiguredValue(ARG_C, "diamond", PathDiamond, commandLineOptions);
	PathSortDiamond = getConfiguredValue(ARG_C, "sortdiamond", PathSortDiamond, commandLineOptions);
        PathBowtie2 = getConfiguredValue(ARG_C, "bowtie2", PathBowtie2, commandLineOptions);
        PathSamtools = getConfiguredValue(ARG_C, "samtools", PathSamtools, commandLineOptions);
        PathBcftools = getConfiguredValue(ARG_C, "bcftools", PathBcftools, commandLineOptions);
	PathExonerate = getConfiguredValue(ARG_C, "exonerate", PathExonerate, commandLineOptions);
        PathMacse = getConfiguredValue(ARG_C, "macse", PathMacse, commandLineOptions);
	PathDelstop = getConfiguredValue(ARG_C, "delstop", PathDelstop, commandLineOptions);
	PathDeltaxa = getConfiguredValue(ARG_C, "deltaxa", PathDeltaxa, commandLineOptions);
        PathTrimal = getConfiguredValue(ARG_C, "trimal", PathTrimal, commandLineOptions);
        PathConcataln = getConfiguredValue(ARG_C, "concataln", PathConcataln, commandLineOptions);
        
        DirRaw = getConfiguredValue(ARG_C, "raw_dir", DirRaw, commandLineOptions);
        DirQcTrim = getConfiguredValue(ARG_C, "fastp_dir", DirQcTrim, commandLineOptions);
        DirAssembly = getConfiguredValue(ARG_C, "spades_dir", DirAssembly, commandLineOptions);
        DirMap = getConfiguredValue(ARG_C, "bowtie2_dir", DirMap, commandLineOptions);
        DirBam = getConfiguredValue(ARG_C, "bam_dir", DirBam, commandLineOptions);
        DirVcf = getConfiguredValue(ARG_C, "vcf_dir", DirVcf, commandLineOptions);
        DirConsensus = getConfiguredValue(ARG_C, "consen_dir", DirConsensus, commandLineOptions);
        DirAlign = getConfiguredValue(ARG_C, "macse_dir", DirAlign, commandLineOptions);
        DirTrim = getConfiguredValue(ARG_C, "trimal_dir", DirTrim, commandLineOptions);

	fastpParas = getValueArray(getValueFromConfig(ARG_C, "fastp_paras"));
	spadesParas = getValueArray(getValueFromConfig(ARG_C, "spades_paras"));
	diamondMakedbParas = getValueArray(getValueFromConfig(ARG_C, "diamond_makedb_paras"));
	diamondBlastxParas = getValueArray(getValueFromConfig(ARG_C, "diamond_blastx_paras"));
	bowtie2BuildParas = getValueArray(getValueFromConfig(ARG_C, "bowtie2_build_paras"));
	bowtie2Paras = getValueArray(getValueFromConfig(ARG_C, "bowtie2_paras"));
	samtoolsViewParas = getValueArray(getValueFromConfig(ARG_C, "samtools_view_paras"));
	samtoolsFixmateParas = getValueArray(getValueFromConfig(ARG_C, "samtools_fixmate_paras"));
	samtoolsSortParas = getValueArray(getValueFromConfig(ARG_C, "samtools_sort_paras"));
	samtoolsMarkdupParas = getValueArray(getValueFromConfig(ARG_C, "samtools_markdup_paras"));
	samtoolsIndexParas = getValueArray(getValueFromConfig(ARG_C, "samtools_index_paras"));
	bcftoolsMpileupParas = getValueArray(getValueFromConfig(ARG_C, "bcftools_mpileup_paras"));
	bcftoolsCallParas = getValueArray(getValueFromConfig(ARG_C, "bcftools_call_paras"));
	bcftoolsNormParas = getValueArray(getValueFromConfig(ARG_C, "bcftools_norm_paras"));
	bcftoolsFilterParas = getValueArray(getValueFromConfig(ARG_C, "bcftools_filter_paras"));
	bcftoolsIndexParas = getValueArray(getValueFromConfig(ARG_C, "bcftools_index_paras"));
	bcftoolsConParas = getValueArray(getValueFromConfig(ARG_C, "bcftools_consensus_paras"));
	exonerateParas = getValueArray(getValueFromConfig(ARG_C, "exonerate_paras"));
	macseParas = getValueArray(getValueFromConfig(ARG_C, "macse_paras"));
	trimalParas = getValueArray(getValueFromConfig(ARG_C, "trimal_paras"));
    }

    // Validate the entire pipeline before starting expensive stages. This catches
    // missing tools (for example sortdiamond or the MACSE JAR) before SPAdes runs.
    if (ARG_F == "all") {
        bool valid = true;
        string[] requiredTools = [PathFastp, PathSpades, PathDiamond, PathSortDiamond,
            PathBowtie2, PathBowtie2 ~ "-build", PathSamtools, PathBcftools, "java",
            PathMacse, PathDelstop, PathDeltaxa, PathTrimal, PathConcataln];
        valid = testFiles(requiredTools) && valid;
        valid = testJava() && valid;
        valid = testInputFiles(ARG_LIST_PATHS ~ [ARG_R] ~ ARG_GENE_PATHS) && valid;
        if (ARG_L.empty) {
            writeln("Error: the sample list is empty.");
            valid = false;
        }
        foreach (sample; ARG_L) {
            string sampleName = getBaseName(sample);
            valid = testInputFiles([
                buildPath(DirRaw, sampleName ~ "_R1.fastq.gz"),
                buildPath(DirRaw, sampleName ~ "_R2.fastq.gz")
            ]) && valid;
        }
        if (enableCodon) valid = testFiles([PathExonerate]) && valid;
        if (!valid) {
            throw new Exception("Preflight validation failed; no pipeline stages were started.");
        }
    }

    writeln("RGBEPP::Start");
    // Perform steps based on provided function argument
    string[] cleanOutputs;
    string[] assemblyOutputs;
    string[] mapOutputs;
    string[] postmapOutputs;
    string[] varcallOutputs;
    string[] consensusOutputs;
    string[] codonOutputs;
    string[] alignOutputs;
    string[] orthologOutputs;
    string[] trimOutputs;
    foreach (sample; ARG_L) {
        string sampleName = getBaseName(sample);
        cleanOutputs ~= buildPath(DirQcTrim, sampleName ~ "_R1.fastq.gz");
        cleanOutputs ~= buildPath(DirQcTrim, sampleName ~ "_R2.fastq.gz");
        assemblyOutputs ~= buildPath(DirAssembly, "scaffolds", sampleName ~ ".fasta");
        assemblyOutputs ~= buildPath(DirAssembly, "contigs", sampleName ~ ".fasta");
        mapOutputs ~= buildPath(DirAssembly, "fasta", sampleName ~ ".fasta");
        mapOutputs ~= buildPath(DirMap, sampleName ~ ".bam");
        postmapOutputs ~= buildPath(DirBam, sampleName ~ ".bam");
        postmapOutputs ~= buildPath(DirBam, sampleName ~ ".bam.bai");
        varcallOutputs ~= buildPath(DirVcf, sampleName ~ ".vcf.gz");
        consensusOutputs ~= buildPath(DirConsensus, "taxa", sampleName ~ ".fasta");
        orthologOutputs ~= buildPath(DirAlign, "taxa", sampleName ~ ".tsv");
    }
    foreach (gene; ARG_G) {
        consensusOutputs ~= buildPath(DirConsensus, "gene", gene ~ ".fasta");
        codonOutputs ~= buildPath(DirConsensus, "gene", gene ~ ".fasta");
        alignOutputs ~= buildPath(DirAlign, "AA", gene ~ ".fasta");
        alignOutputs ~= buildPath(DirAlign, "NT", gene ~ ".fasta");
        trimOutputs ~= buildPath(DirTrim, "NT", gene ~ ".fasta");
    }
    orthologOutputs ~= strParalog;
    string[] concatOutputs = [buildPath(DirTrim, "concat.fasta.fasta")];

    if ((ARG_F == "all" || ARG_F == "clean") && !skipCompletedStage("QcTrimming", cleanOutputs, enableResume)) {
	if(testFiles([PathFastp]) && testStringArray(ARG_L)){
          beginCommandBatch("QcTrimming");
          processQcTrim(ARG_L, ARG_T, DirRaw, DirQcTrim, PathFastp, fastpParas); //ARG_L
          finishCommandBatch();
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if ((ARG_F == "all" || ARG_F == "assembly") && !skipCompletedStage("Assembly", assemblyOutputs, enableResume)) {
        if(testFiles([PathSpades]) && testStringArray(ARG_L)){
          beginCommandBatch("Assembly");
	  processAssembly(ARG_L, ARG_M, ARG_T, DirQcTrim, DirAssembly, PathSpades, spadesParas); //ARG_L
	  processAssemMv(ARG_L, DirAssembly);
          finishCommandBatch();
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if ((ARG_F == "all" || ARG_F == "map") && !skipCompletedStage("Mapping", mapOutputs, enableResume)) {
	if(testFiles([PathBowtie2, PathDiamond, PathSamtools, PathSortDiamond]) && testStringArray(ARG_L) && testString(ARG_R)  ){
          beginCommandBatch("Mapping");
	  processMappingDenovo(ARG_L, ARG_R, ARG_T, DirQcTrim, DirAssembly, DirMap, PathBowtie2, PathDiamond, PathSamtools, PathSortDiamond, diamondMakedbParas, diamondBlastxParas, bowtie2BuildParas, bowtie2Paras, samtoolsViewParas); //ARG_L, ARG_R
          finishCommandBatch();
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if ((ARG_F == "all" || ARG_F == "postmap") && !skipCompletedStage("PostMapping", postmapOutputs, enableResume)) {
	if(testFiles([PathSamtools]) && testStringArray(ARG_L) ){
          beginCommandBatch("PostMapping");
          processPostMap(ARG_L, ARG_T, DirMap, DirBam, PathSamtools, samtoolsFixmateParas, samtoolsSortParas, samtoolsMarkdupParas, samtoolsIndexParas); //ARG_L
          finishCommandBatch();
        } else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if ((ARG_F == "all" || ARG_F == "varcall") && !skipCompletedStage("VarCalling", varcallOutputs, enableResume)) {
	if(testFiles([PathBcftools]) && testStringArray(ARG_L) ){
            beginCommandBatch("VarCalling");
	    processVarCallDenovo(ARG_L, ARG_T, DirAssembly, DirMap, DirBam, DirVcf, PathBcftools, bcftoolsMpileupParas, bcftoolsCallParas, bcftoolsNormParas, bcftoolsFilterParas); //ARG_L
            finishCommandBatch();
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if ((ARG_F == "all" || ARG_F == "consen") && !skipCompletedStage("Consensus", consensusOutputs, enableResume)) {
	if(testFiles([PathBcftools]) && testStringArray(ARG_L) && testStringArray(ARG_G) ){
          beginCommandBatch("Consensus");
	  processConDenovo(ARG_G, ARG_L, ARG_T, DirAssembly, DirVcf, DirConsensus, PathBcftools, bcftoolsIndexParas, bcftoolsConParas); //ARG_G ARG_L 
	  processCombFasta(ARG_G, ARG_L, DirConsensus, strTaxa, strGene); //ARG_G ARG_L
          finishCommandBatch();
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if (((ARG_F == "all" && enableCodon) || ARG_F == "codon") && !skipCompletedStage("Codon", codonOutputs, enableResume)) {
	if(testFiles([PathExonerate]) && testStringArray(ARG_G) && testString(ARG_R)){
          beginCommandBatch("Codon");
	  processCodon(ARG_G, ARG_R, DirConsensus, PathExonerate); //ARG_G
          finishCommandBatch();
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if ((ARG_F == "all" || ARG_F == "align") && !skipCompletedStage("Align", alignOutputs, enableResume)) {
	if(testFiles([PathMacse]) && testJava && testStringArray(ARG_G)){
          beginCommandBatch("Align");
	  processAlign(ARG_G, DirConsensus, DirAlign, PathMacse, macseParas); //ARG_G
          finishCommandBatch();
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if ((ARG_F == "all" || ARG_F == "ortholog") && !skipCompletedStage("Ortholog", orthologOutputs, enableResume)) {
	if(testFiles([PathDeltaxa]) && testFiles([PathDiamond]) && testStringArray(ARG_G) && testStringArray(ARG_L)){
          beginCommandBatch("Ortholog");
	  processCombFasta(ARG_L, ARG_G, DirAlign, "AA", strTaxa);
	  processOrtholog(ARG_L, ARG_R, DirAlign, PathDeltaxa, PathDiamond, strParalog); //ARG_L ARG_R
          finishCommandBatch();
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if ((ARG_F == "all" || ARG_F == "trim") && !skipCompletedStage("Trimming", trimOutputs, enableResume)) {
	if(testFiles([PathTrimal]) && testStringArray(ARG_G) ){
          beginCommandBatch("Trimming");
	  processTrimming(ARG_G, DirAlign, DirTrim, PathDelstop, PathTrimal, trimalParas); //ARG_G
          finishCommandBatch();
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if ((ARG_F == "all" || ARG_F == "concat") && !skipCompletedStage("Concatenate", concatOutputs, enableResume)) {
        if(testFiles([PathConcataln]) && testStringArray(ARG_G)){
          beginCommandBatch("Concatenate");
	  processConcat(ARG_G, DirTrim, PathConcataln);
          finishCommandBatch();
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    writeln("RGBEPP::End");
}
