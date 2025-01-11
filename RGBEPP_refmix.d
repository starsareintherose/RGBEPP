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

void show_help(string pkgver) {
    writeln("\t\t\t\t\t\033[0;47;31mR\033[0m\033[0;47;92mG\033[0m\033[0;47;94mB\033[0m\033[0;47m \033[0m\033[0;47;33mE\033[0m\033[0;47;94mP\033[0m\033[0;47;33mP\033[0m
\t\t\tReference Genome based Exon Phylogeny Pipeline
	    Version: ", pkgver, "
	    License: GPL-2.0-only
	    Author: Guoyi Zhang
	    -c\t--config\tconfig file for software path (optional)
	    -f\t--functions\tfunctions type (optional): all clean assembly 
	      \t           \t map postmap varcall consen codon align trim
	    -g\t--genes\t\tgene file path (optional, if -r is specified)
	    -h\t--help\t\tshow this information
	    -l\t--list\t\tlist file path
	    -m\t--memory\tmemory settings (optional, default 16 GB)
	    -r\t--reference\treference genome path
	    -t\t--threads\tthreads setting (optional, default 8 threads)
	    --codon\t\tOnly use the codon region (optional)
	    --fastp\t\tFastp path (optional)
	    --spades\t\tSpades python path (optional)
	    --diamond\t\tDiamond python path (optional)
	    --sortdiamond\tSortDiamond python path (optional)
	    --bowtie2\t\tBowtie2 path (optional)
	    --samtools\t\tSamtools path (optional)
	    --bcftools\t\tBcftools path (optional)
	    --exonerate\t\tExonerate path (optional)
	    --macse\t\tMacse jarfile path (optional)
	    --delstop\t\tDelstop path (optional)
	    --trimal\t\tTrimal path (optional)
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
        if (!exists(filePath) && filePath != "") {
            writeln("Error: " ~ filePath ~ " does not exists.");
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

    if (wait(process) != 0) {
        writeln("Error executing command: ", cmd.join(" "));
    }
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
    wait(pid);

    // Write the output to the specified file
    std.file.write(outputFile, cast(ubyte[])output);
}


void executeCommandPipe(string[][] cmds) {

    Pid[] pids;
    scope(exit) {
	    foreach (pid; pids) {
		    wait(pid);
	    }
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
    auto regex = regex(key ~ r"\s*=\s*(.+)");

    foreach (line; content.splitter("\n")) {
        if (auto match = matchFirst(line, regex)) {
            value = match.captures[1];
            break;
        }
    }

    return value;
}

void processQcTrim(string[] ARG_L, int ARG_T, string DirRaw, string DirQcTrim, string PathFastp) {
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
                        "-w", ARG_T.to!string];
        executeCommand(cmdQcTrim);
    }
    writeln("QcTrimming::End");
}

void processReduce(string[] ARG_L, int ARG_T, string DirQcTrim, string PathFastUniq, string PathPigz, string PathUnpgiz, string PathMkfifo){
    writeln("Reducing::Start");
    foreach (string file; parallel(ARG_L, 1)) {
        string baseName = getBaseName(file);
	string fastq1 = buildPath(DirQcTrim, baseName ~ "_R1.fastq.gz");
	string fastq2 = buildPath(DirQcTrim, baseName ~ "_R2.fastq.gz");
	string fastuniqTxt = buildPath(DirQcTrim, baseName ~ ".txt");
	if(exists(fastq1) && exists(fastq2)){

	string[] cmdDmMakeUnzip1 = [PathUnpgiz, baseName ~ "_R1.fastq.gz"];
        string[] cmdDmMakeUnzip2 = [PathUnpgiz, baseName ~ "_R2.fastq.gz"];
        string[] cmdReduce = [PathFastUniq, "-i", baseName, "-o", baseName  ~ "_R2.fastq" , "-p", baseName ~ "_R1.fastq" ];
	} else {
            writeln("Please check: ", fastq1, " and ", fastq2);
            continue;	
	}
    }
    writeln("Reducing::End");
}

void processAssembly(string[] ARG_L, int ARG_M, int ARG_T, string DirQcTrim, string DirAssembly, string PathSpades){
    writeln("Assembly::Start");
    createDir(DirAssembly);
    foreach (string file; ARG_L) {
       string baseName = getBaseName(file);
       string DirAss = buildPath(DirAssembly, baseName);
       createDir(DirAss);
       string inputFileR1 = buildPath(DirQcTrim, baseName ~ "_R1.fastq.gz");
       string inputFileR2 = buildPath(DirQcTrim, baseName ~ "_R2.fastq.gz");
       string[] cmdAssembly = [PathSpades, "--pe1-1", inputFileR1, "--pe1-2", inputFileR2, "-t", ARG_T.to!string, "-m", ARG_M.to!string, "--careful", "--phred-offset", "33", "-o", DirAss];
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

void processMapping(string[] ARG_L, string ARG_R, int ARG_T, string DirQcTrim, string DirMap, string PathBowtie2, string PathSamtools) {
    writeln("Mapping::Start");

    // Prepare directory
    createDir(DirMap);

    createDir(DirMap ~ "/index");
    string PathBowtie2_build = PathBowtie2 ~ "-build";
    string[] Refs = getRef(ARG_R, DirMap);
    string ARG_R_index = Refs[0]; // bt2_index_base
    string ARG_R_refer = Refs[1]; //reference_in fasta file

    copy(ARG_R, ARG_R_refer);

    string[] cmdBuildDB = [PathBowtie2_build, "--threads", ARG_T.to!string, ARG_R_refer,  ARG_R_index];
    executeCommand(cmdBuildDB);

    foreach (string file; ARG_L) {
        string baseName = baseName(file, ".fastq.gz");
        string outputBam = DirMap ~ "/" ~ baseName ~ ".bam";
        string inputFileR1 = DirQcTrim ~ "/" ~ baseName ~ "_R1.fastq.gz";
        string inputFileR2 = DirQcTrim ~ "/" ~ baseName ~ "_R2.fastq.gz";

        // Perform mapping using Bowtie2 and converted to Bam using samtools
        string[] cmdMap = [PathBowtie2, "-x", ARG_R_index, "-1", inputFileR1, "-2", inputFileR2,
                        "-p", ARG_T.to!string];
        string[] cmdSam2Bam = [PathSamtools, "view", "-bS", "-@", ARG_T.to!string, "-o", outputBam];

	executeCommandPipe([cmdMap, cmdSam2Bam]);

    }
    writeln("Mapping::End");
}

void processMappingDenovo(string[] ARG_L, string ARG_R, int ARG_T, string DirQcTrim, string DirAssembly, string DirMap, string PathBowtie2, string PathDiamond, string PathSamtools, string PathSortDiamond){
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
    string [] cmdDmMakeDB = [ PathDiamond, "makedb", "--db", "Reference", "--in", ARG_R_Ref];
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

	string[] cmdDiamond = [PathDiamond, "blastx", "-d", "Reference.dmnd", "-q", inputFasta, "-o", inputM8, "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "gaps", "ppos", "qframe", "qseq"];
	string[] cmdSortDiamond = [PathSortDiamond, inputM8, outputSort];
    	string[] cmdBuildDB = [PathBowtie2_build, "--threads", ARG_T.to!string, outputSort, outputIndex];
        string[] cmdMap = [PathBowtie2, "-x", outputIndex, "-1", inputFileR1, "-2", inputFileR2, "-p", ARG_T.to!string];
        string[] cmdSam2Bam = [PathSamtools, "view", "-bS", "-@", ARG_T.to!string, "-o", outputBam];
	executeCommand(cmdDiamond);
	executeCommand(cmdSortDiamond);
    	executeCommand(cmdBuildDB);
	executeCommandPipe([cmdMap, cmdSam2Bam]);
    }

    writeln("Mapping::End");
}

void processPostMap(string[] ARG_L, int ARG_T, string DirMap, string DirBam, string PathSamtools) {

    createDir(DirBam);
    writeln("PostMapping::Start");

    foreach (string file; ARG_L) {
        string baseName = getBaseName(file);
        string inputBam = buildPath(DirMap, baseName ~ ".bam");
        string outputBam = buildPath(DirBam, baseName ~ ".bam");

        // Convert SAM to BAM, sort and remove duplicates using Samtools
	string[] cmdFixmate = [PathSamtools, "fixmate", "-@", ARG_T.to!string, "-m", inputBam, "-"];
        string[] cmdSort = [PathSamtools, "sort", "-@", ARG_T.to!string, "-"];
        string[] cmdMarkdup = [PathSamtools, "markdup", "-@", ARG_T.to!string, "-", outputBam];
	executeCommandPipe([cmdFixmate, cmdSort, cmdMarkdup]);

        string [] cmdIndexBam = [PathSamtools, "index", "-@", ARG_T.to!string, outputBam];
        executeCommand(cmdIndexBam);
    }

    writeln("PostMapping::End");
}

void processVarCall(string[] ARG_L, string ARG_R, int ARG_T, string DirMap, string DirBam, string DirVcf, string PathBcftools) {
    writeln("VarCalling::Start");

    string[] Refs = getRef(ARG_R, DirMap);
    string ARG_R_refer = Refs[1]; //reference_in fasta file

    createDir(DirVcf);

    foreach (string file; parallel(ARG_L, 1)) {
        string baseName = getBaseName(file);
        string inputBam = DirBam ~ "/" ~ baseName ~ ".bam";
        string outputVcf = DirVcf ~ "/" ~ baseName ~ ".vcf.gz";

        // Variant calling using bcftools
        string[] cmdPileup = [PathBcftools, "mpileup", "-Oz", "--threads", ARG_T.to!string, "-f", ARG_R_refer, inputBam];
	string[] cmdVarCall = [PathBcftools, "call", "-mv", "-Oz", "--threads", ARG_T.to!string];
	string[] cmdNorm = [PathBcftools, "norm", "--threads", ARG_T.to!string, "-f", ARG_R_refer, "--check-ref", "s", "-Oz"];
	string[] cmdFilter = [PathBcftools, "filter", "--threads", ARG_T.to!string, "--IndelGap", "5", "-Oz", "-o", outputVcf];
        executeCommandPipe([cmdPileup, cmdVarCall, cmdNorm, cmdFilter]); 
    }

    writeln("VarCalling::End");

}

void processVarCallDenovo(string[] ARG_L, int ARG_T, string DirAssembly, string DirMap, string DirBam, string DirVcf, string PathBcftools) {
    writeln("VarCalling::Start");

    string DirAssemblyFas = buildPath(DirAssembly, "fasta");
    createDir(DirVcf);

    foreach (string file; parallel(ARG_L, 1)) {
        string baseName = getBaseName(file);
        string inputBam = buildPath(DirBam, baseName ~ ".bam");
        string outputVcf = buildPath(DirVcf, baseName ~ ".vcf.gz");
	string referFasta = buildPath(DirAssemblyFas, baseName ~ ".fasta");
        // Variant calling using bcftools
        string[] cmdPileup = [PathBcftools, "mpileup", "-Oz", "--threads", ARG_T.to!string, "-f", referFasta, inputBam];
	string[] cmdVarCall = [PathBcftools, "call", "-mv", "-Oz", "--threads", ARG_T.to!string];
	string[] cmdNorm = [PathBcftools, "norm", "--threads", ARG_T.to!string, "-f", referFasta, "-Oz"];
	string[] cmdFilter = [PathBcftools, "filter", "--threads", ARG_T.to!string, "--IndelGap", "5", "-Oz", "-o", outputVcf];
        executeCommandPipe([cmdPileup, cmdVarCall, cmdNorm, cmdFilter]); 
    }

    writeln("VarCalling::End");

}

void processCon(string[] ARG_G, string[] ARG_L, string ARG_R, int ARG_T, string DirMap,  string DirVcf, string DirConsensus, string PathBcftools) {
    createDir(DirConsensus);

    string DirConTaxa = DirConsensus ~ "/" ~ "taxa";

    createDir(DirConTaxa);

    string[] Refs = getRef(ARG_R, DirMap);
    string ARG_R_refer = Refs[1]; //reference_in fasta file

    writeln("Consensus::Start");
    // Extract fasta from vcf file
    foreach (string file; ARG_L) {
        string baseName = getBaseName(file);
        string inputVcf = DirVcf ~ "/" ~ baseName ~ ".vcf.gz";
        string outputFasta = DirConTaxa ~ "/" ~ baseName ~ ".fasta";

	// index vcf.gz
	string[] cmdIndexVcf = [PathBcftools, "index", inputVcf];
	executeCommand(cmdIndexVcf);
	
        // Generate consensus sequences using bcftools
        string[] cmdCon = [PathBcftools, "consensus", "-f", ARG_R, inputVcf, "-o", outputFasta];
	executeCommand(cmdCon);
    }
    // Recombine the sequences based on genes
    writeln("Consensus::End");
}

void processConDenovo(string[] ARG_G, string[] ARG_L, int ARG_T, string DirAssembly,  string DirVcf, string DirConsensus, string PathBcftools) {
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
	string[] cmdIndexVcf = [PathBcftools, "index", inputVcf];
	executeCommand(cmdIndexVcf);
	
        // Generate consensus sequences using bcftools
        string[] cmdCon = [PathBcftools, "consensus", "-f", referFasta, inputVcf, "-o", outputFasta];
	executeCommand(cmdCon);
    }
    // Recombine the sequences based on genes
    writeln("Consensus::End");
}


void processCombFasta(string[] ARG_G, string[] ARG_L, string DirConsensus) {

    string DirConTaxa = buildPath(DirConsensus, "taxa");
    string DirConGene = buildPath(DirConsensus, "gene");    
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

void processAlign(string[] ARG_G, string DirConsensus, string DirAlign, string PathMacse){

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
    	    string[] cmdAlign = ["java", "-jar", PathMacse, "-prog", "alignSequences", "-seq" , inputFasta, "-out_AA", outAA, "-out_NT", outNT ];
    	    executeCommand(cmdAlign);
	}
    }
    writeln("Align::End");

}

void processTrimming(string[] ARG_G, string DirAlign, string DirTrim, string PathDelstop, string PathTrimal){
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
            string[] cmdTrim = [PathTrimal, "-in", inputFastaAA, "-backtrans", inputBackTransNT, "-out", outputFastaNT, "-automated1"];
            executeCommand(cmdTrim);
        } else {
            writeln("Skipping gene: ", gene, " as files are missing.");
        }
    }
    writeln("Trimming::End");

}

void main(string[] args) {
    string pkgver = "0.0.3";

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

    string PathFastp = "/usr/bin/fastp";
    string PathSpades = "/usr/bin/spades.py";
    string PathDiamond = "/usr/bin/diamond";
    string PathSortDiamond = "/usr/bin/sortdiamond";
    string PathBowtie2 = "/usr/bin/bowtie2";
    string PathSamtools = "/usr/bin/samtools";
    string PathBcftools = "/usr/bin/bcftools";
    string PathExonerate = "/usr/bin/exonerate";
    string PathMacse = "/usr/share/java/macse.jar";
    string PathDelstop = "/usr/bin/delstop";
    string PathTrimal = "/usr/bin/trimal";

    int ARG_T = 8;
    int ARG_M = 16;
    string[] ARG_G;
    string[] ARG_L;
    string ARG_C;
    string ARG_F;
    string ARG_R;
    bool enableCodon = false;

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
		case "-g", "--gene":
		    i++;
                    ARG_G ~= readArrFromFile(args[i]);
		    break;
                case "-h", "--help":
                    show_help(pkgver);
                    return;
                case "-l", "--list":
		    i++;
                    ARG_L ~= readArrFromFile(args[i]);
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
                    break;
                case "--spades":
		    i++;
                    PathSpades = args[i];
                    break;
                case "--diamond":
		    i++;
                    PathDiamond = args[i];
                    break;
		case "--sortdiamond":
		    i++;
                    PathSortDiamond = args[i];
                    break;
                case "--bowtie2":
		    i++;
                    PathBowtie2 = args[i];
                    break;
                case "--samtools":
		    i++;
                    PathSamtools = args[i];
                    break;
                case "--bcftools":
		    i++;
                    PathBcftools = args[i];
                    break;
                case "--exonerate":
		    i++;
                    PathExonerate = args[i];
                    break;
                case "--macse":
		    i++;
                    PathMacse = args[i];
                    break;
                case "--delstop":
		    i++;
                    PathDelstop = args[i];
                    break;
                case "--trimal":
		    i++;
                    PathTrimal = args[i];
                    break;
                default:
                    break;
            }
        }
    } else {
        show_help(pkgver);
        return;
    }

    // get gene from ARG_R reference fasta
    if (ARG_R.length != 0 && ARG_G.length == 0 ){
    	ARG_G = getARG_G(ARG_R);
    }
  
    if (ARG_F.length == 0 ){
    	ARG_F = "all";
    }

    // get pathXXX form config file
    if (ARG_C.length != 0){
    
        PathFastp = getValueFromConfig(ARG_C, "fastp");
	PathSpades = getValueFromConfig(ARG_C, "spades");
	PathDiamond = getValueFromConfig(ARG_C, "diamond");
	PathSortDiamond = getValueFromConfig(ARG_C, "sortdiamond");
        PathBowtie2 = getValueFromConfig(ARG_C, "bowtie2");
        PathSamtools = getValueFromConfig(ARG_C, "samtools");
        PathBcftools = getValueFromConfig(ARG_C, "bcftools");
	PathExonerate = getValueFromConfig(ARG_C, "exonerate");
        PathMacse = getValueFromConfig(ARG_C, "macse");
	PathDelstop = getValueFromConfig(ARG_C, "delstop");
        PathTrimal = getValueFromConfig(ARG_C, "trimal");
        
	DirRaw = getValueFromConfig(ARG_C, "raw_dir");
        DirQcTrim = getValueFromConfig(ARG_C, "fastp_dir");
        DirAssembly = getValueFromConfig(ARG_C, "spades_dir");
        DirMap = getValueFromConfig(ARG_C, "bowtie2_dir");
        DirBam = getValueFromConfig(ARG_C, "bam_dir");
        DirVcf = getValueFromConfig(ARG_C, "vcf_dir");
        DirConsensus = getValueFromConfig(ARG_C, "consen_dir");
        DirAlign = getValueFromConfig(ARG_C, "macse_dir"); 
        DirTrim = getValueFromConfig(ARG_C, "trimal_dir"); 

    }

    writeln("RGBEPP::Start");
    // Perform steps based on provided function argument
    if (ARG_F == "all" || ARG_F == "clean") {
	if(testFiles([PathFastp]) && testStringArray(ARG_L)){
          processQcTrim(ARG_L, ARG_T, DirRaw, DirQcTrim, PathFastp); //ARG_L
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if (ARG_F == "all" || ARG_F == "assembly") {
        if(testFiles([PathSpades]) && testStringArray(ARG_L)){
	  processAssembly(ARG_L, ARG_M, ARG_T, DirQcTrim, DirAssembly, PathSpades); //ARG_L
	  processAssemMv(ARG_L, DirAssembly);
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if (ARG_F == "all" || ARG_F == "map") {
	if(testFiles([PathBowtie2, PathDiamond, PathSamtools, PathSortDiamond]) && testStringArray(ARG_L) && testString(ARG_R)  ){
	  processMapping(ARG_L, ARG_R, ARG_T, DirQcTrim, DirMap, PathBowtie2, PathSamtools); 
	  //processMappingDenovo(ARG_L, ARG_R, ARG_T, DirQcTrim, DirAssembly, DirMap, PathBowtie2, PathDiamond, PathSamtools, PathSortDiamond); //ARG_L, ARG_R
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if (ARG_F == "all" || ARG_F == "postmap") {
	if(testFiles([PathSamtools]) && testStringArray(ARG_L) ){
          processPostMap(ARG_L, ARG_T, DirMap, DirBam, PathSamtools); //ARG_L
        } else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if (ARG_F == "all" || ARG_F == "varcall") {
	if(testFiles([PathBcftools]) && testStringArray(ARG_L) ){
	    processVarCall(ARG_L, ARG_R, ARG_T, DirMap, DirBam, DirVcf, PathBcftools);
	    //processVarCallDenovo(ARG_L, ARG_T, DirAssembly, DirMap, DirBam, DirVcf, PathBcftools); //ARG_L
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if (ARG_F == "all" || ARG_F == "consen") {
	if(testFiles([PathBcftools]) && testStringArray(ARG_L) && testStringArray(ARG_G) ){
	  processCon(ARG_G, ARG_L, ARG_R, ARG_T, DirMap, DirVcf, DirConsensus, PathBcftools);
	  //processConDenovo(ARG_G, ARG_L, ARG_T, DirAssembly, DirVcf, DirConsensus, PathBcftools); //ARG_G ARG_L 
	  processCombFasta(ARG_G, ARG_L, DirConsensus); //ARG_G ARG_L
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if (ARG_F == "all" && enableCodon || ARG_F == "codon") {
	if(testFiles([PathExonerate]) && testStringArray(ARG_G) && testString(ARG_R)){
	  processCodon(ARG_G, ARG_R, DirConsensus, PathExonerate); //ARG_G
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if (ARG_F == "all" || ARG_F == "align") {
	if(testFiles([PathMacse]) && testJava && testStringArray(ARG_G)){
	  processAlign(ARG_G, DirConsensus, DirAlign, PathMacse); //ARG_G
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    if (ARG_F == "all" || ARG_F == "trim") {
	if(testFiles([PathTrimal]) && testStringArray(ARG_G) ){
	  processTrimming(ARG_G, DirAlign, DirTrim, PathDelstop, PathTrimal); //ARG_G
	} else {
	  throw new Exception("please confirm paramenters are correct"); 
	}
    }

    writeln("RGBEPP::End");
}

