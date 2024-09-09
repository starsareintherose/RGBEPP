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
	    -g\t--genes\t\tgene file path (optional, if -r is specified)
	    -f\t--functions\tfunctions type (optional): all clean map 
	      \t           \tpostmap varcall consen align
	    -h\t--help\t\tshow this information
	    -l\t--list\t\tlist file path
	    -m\t--memory\tmemory settings (optional, default 16 GB)
	    -r\t--reference\treference genome path
	    -t\t--threads\tthreads setting (optional, default 8 threads)
	    --fastp\t\tFastp path (optional)
	    --bowtie2\t\tBowtie2 path (optional)
	    --samtools\t\tSamtools path (optional)
	    --bcftools\t\tBcftools path (optional)
	    --macse\t\tMacse jarfile path (optional)
	    --trimal\t\tTrimal path (optional)
	    --spades\t\tSpades python path (optional)
	    for example: ./RGBEPP -f all -l list -t 8 -r reference.fasta \n");
}

void createDir(string path) {
    if (!exists(path)) {
        mkdir(path);
    }
}

void executeCommand(string[] cmd) {
    auto process = spawnProcess(cmd);

    if (wait(process) != 0) {
        writeln("Error executing command: ", cmd.join(" "));
    }
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
    string ARG_R_index = DirMap ~ "/index/" ~ baseNameRef; // bt2_index_base
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
        string inputFileR1 = DirRaw ~ "/" ~ baseName ~ "_R1.fastq.gz";
        string inputFileR2 = DirRaw ~ "/" ~ baseName ~ "_R2.fastq.gz";
        string outputFileR1 = DirQcTrim ~ "/" ~ baseName ~ "_R1.fastq.gz";
        string outputFileR2 = DirQcTrim ~ "/" ~ baseName ~ "_R2.fastq.gz";
        string jsonFile = DirQcTrim ~ "/" ~ baseName ~ ".json";
        string htmlFile = DirQcTrim ~ "/" ~ baseName ~ ".html";

        // Perform quality control and trimming using external program `fastp`
        string[] cmdQcTrim = [PathFastp, "-i", inputFileR1, "-I", inputFileR2,
                        "-o", outputFileR1, "-O", outputFileR2,
                        "-j", jsonFile, "-h", htmlFile,
                        "-w", ARG_T.to!string];
        executeCommand(cmdQcTrim);
    }
    writeln("QcTrimming::End");
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

void processPostMap(string[] ARG_L, int ARG_T, string DirMap, string DirBam, string PathSamtools) {

    createDir(DirBam);
    writeln("PostMapping::Start");

    foreach (string file; ARG_L) {
        string baseName = getBaseName(file);
        string inputBam = DirMap ~ "/" ~ baseName ~ ".bam";
        string outputBam = DirBam ~ "/" ~ baseName ~ ".bam";

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

    foreach (string file; ARG_L) {
        string baseName = getBaseName(file);
        string inputBam = DirBam ~ "/" ~ baseName ~ ".bam";
        string outputVcf = DirVcf ~ "/" ~ baseName ~ ".vcf.gz";

        // Variant calling using bcftools
        string[] cmdPileup = [PathBcftools, "mpileup", "-Oz", "--threads", ARG_T.to!string, "-f", ARG_R_refer, inputBam];
	string[] cmdVarCall = [PathBcftools, "call", "-mv", "-Oz", "--threads", ARG_T.to!string];
	string[] cmdNorm = [PathBcftools, "norm", "--threads", ARG_T.to!string, "-f", ARG_R_refer, "-Oz"];
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

void processCombFasta(string[] ARG_G, string[] ARG_L, string DirConsensus) {

    string DirConTaxa = DirConsensus ~ "/" ~ "taxa";
    string DirConGene = DirConsensus ~ "/" ~ "gene";    
    createDir(DirConGene);

    // create a dictory
    string[string] geneSequences;

    writeln("ConvertFasta::Start");
    // read first
    foreach (file; ARG_L) {
        string inputFile = DirConTaxa ~ "/" ~ file ~ ".fasta";
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
        string outputFile = DirConGene ~ "/" ~ gene ~ ".fasta";
        File output = File(outputFile, "w");
        if (gene in geneSequences) {
            output.write(geneSequences[gene]);
        }
    }
    writeln("ConvertFasta::End");
}

void processAlign(string[] ARG_G, string DirConsensus, string DirAlign, string PathMacse){

    string DirConGene = DirConsensus ~ "/" ~ "gene";
    string DirAlignAA = DirAlign ~ "/" ~ "AA";
    string DirAlignNT = DirAlign ~ "/" ~ "NT";

    writeln("Align::Start");
    createDir(DirAlign);
    createDir(DirAlignAA);
    createDir(DirAlignNT);
    foreach (gene; parallel(ARG_G, 1)) {
    	string inputFasta = DirConGene ~ "/" ~ gene ~ ".fasta";
    	string outAA = DirAlignAA ~ "/" ~ gene ~ ".fasta";
    	string outNT = DirAlignNT ~ "/" ~ gene ~ ".fasta";
    	string[] cmdAlign = ["java", "-jar", PathMacse, "-prog", "alignSequences", "-seq" , inputFasta, "-out_AA", outAA, "-out_NT", outNT ];
    	executeCommand(cmdAlign);
    }
    writeln("Align::End");

}

void processTrimming(string[] ARG_G, string DirAlign, string DirTrim, string PathDelstop, string PathTrimal){
    createDir(DirAlign ~ "/" ~ "AA_out");
    createDir(DirAlign ~ "/" ~ "NT_out");

    string DirAA = DirAlign ~ "/" ~ "AA";
    string DirNT = DirAlign ~ "/" ~ "NT";

    // copy file firstly
    foreach (gene; ARG_G){
   	string inputFastaAA = DirAA ~ "/" ~ gene ~ ".fasta";
	string outputFastaAA = DirAA ~ "_out" ~ "/" ~ gene ~ ".fasta";
   	string inputFastaNT = DirNT ~ "/" ~ gene ~ ".fasta";
	string outputFastaNT = DirNT ~ "_out" ~ "/" ~ gene ~ ".fasta";
   
        copy(inputFastaNT, outputFastaNT);
        copy(inputFastaAA, outputFastaAA);  
        // del stop codon
        string[] cmdDelStop = [PathDelstop, outputFastaAA, outputFastaNT, "--delete"];
        executeCommand(cmdDelStop);	
    }

    createDir(DirTrim);
    foreach (gene; ARG_G){
        string inputFastaAA = DirAA ~ "_out" ~ "/" ~ gene ~ ".fasta";
	string inputBackTransNT = DirNT ~ "_out" ~ "/" ~ gene ~ ".fasta";
	string outputFastaNT = DirTrim ~ "/" ~ gene ~ ".fasta";
	
   	string[] cmdTrim = [PathTrimal, "-in", inputFastaAA, "-backtrans", inputBackTransNT, "-out", outputFastaNT, "-gappyout"];
        executeCommand(cmdTrim);	
    }
}

void processAssembly(string[] ARG_L, int ARG_M, int ARG_T, string DirQcTrim, string DirAssembly, string PathSpades){
    writeln("Assembly::Start");
    createDir(DirAssembly);
    foreach (string file; ARG_L) {
       string baseName = getBaseName(file);
       string DirAss = DirAssembly ~ "/" ~ baseName;
       createDir(DirAss);
       string inputFileR1 = DirQcTrim ~ "/" ~ baseName ~ "_R1.fastq.gz";
       string inputFileR2 = DirQcTrim ~ "/" ~ baseName ~ "_R2.fastq.gz";
       string[] cmdAssembly = [PathSpades, "--pe1-1", inputFileR1, "--pe1-2", inputFileR2, "-t", ARG_T.to!string, "-m", ARG_M.to!string, "--careful", "--phred-offset", "33", "-o", DirAss];
    	executeCommand(cmdAssembly);
    }
    writeln("Assembly::End");
}

void main(string[] args) {
    string pkgver = "0.0.3";

    string DirHome = std.file.getcwd();
    string DirRaw = DirHome ~ "/00_raw";
    string DirQcTrim = DirHome ~ "/01_fastp";
    string DirMap = DirHome ~ "/02_bowtie2";
    string DirAssembly = DirHome ~ "/02_spades";
    string DirBam = DirHome ~ "/03_bam";
    string DirVcf = DirHome ~ "/04_vcf";
    string DirConsensus = DirHome ~ "/05_consen";
    string DirAlign = DirHome ~ "/06_macse"; 
    string DirTrim = DirHome ~ "/07_trimal"; 

    string PathFastp = "/usr/bin/fastp";
    string PathBowtie2 = "/usr/bin/bowtie2";
    string PathSamtools = "/usr/bin/samtools";
    string PathBcftools = "/usr/bin/bcftools";
    string PathMacse = "/usr/share/java/macse.jar";
    string PathDelstop = "/usr/bin/delstop";
    string PathTrimal = "/usr/bin/trimal";
    string PathSpades = "/usr/bin/spades.py";

    int ARG_T = 8;
    int ARG_M = 16;
    string[] ARG_G;
    string[] ARG_L;
    string ARG_C;
    string ARG_F;
    string ARG_R;
   
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
                case "-r", "--reference":
		    i++;
                    ARG_R = args[i];
                    break;
                case "-t", "--threads":
		    i++;
                    ARG_T = args[i].to!int;
                    break;
                case "--fastp":
		    i++;
                    PathFastp = args[i];
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
                case "--macse":
		    i++;
                    PathMacse = args[i];
                    break;
                case "--trimal":
		    i++;
                    PathTrimal = args[i];
                    break;
                case "--spades":
		    i++;
                    PathSpades = args[i];
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
    if (ARG_R.length != 0 ){
    	ARG_G = getARG_G(ARG_R);
    }
   
    // get pathXXX form config file
    if (ARG_C != ""){
    
        PathFastp = getValueFromConfig(ARG_C, "fastp");
        PathBowtie2 = getValueFromConfig(ARG_C, "bowtie2");
        PathSamtools = getValueFromConfig(ARG_C, "samtools");
        PathBcftools = getValueFromConfig(ARG_C, "bcftools");
        PathMacse = getValueFromConfig(ARG_C, "macse");
        PathTrimal = getValueFromConfig(ARG_C, "trimal");

    }

    writeln("RGBEPP::Start");
    // Perform steps based on provided function argument
    if (ARG_F == "all" || ARG_F == "clean") {
        processQcTrim(ARG_L, ARG_T, DirRaw, DirQcTrim, PathFastp);
    }

    if (ARG_F == "assembly") {
	processAlign(ARG_G, DirConsensus, DirAlign, PathMacse);
    }

    if (ARG_F == "all" || ARG_F == "map") {
        processMapping(ARG_L, ARG_R, ARG_T, DirQcTrim, DirMap, PathBowtie2, PathSamtools);
    }

    if (ARG_F == "all" || ARG_F == "postmap") {
        processPostMap(ARG_L, ARG_T, DirMap, DirBam, PathSamtools);
    }

    if (ARG_F == "all" || ARG_F == "varcall") {
        processVarCall(ARG_L, ARG_R, ARG_T, DirMap, DirBam, DirVcf, PathBcftools);
    }

    if (ARG_F == "all" || ARG_F == "consen") {
        processCon(ARG_G, ARG_L, ARG_R, ARG_T, DirMap, DirVcf, DirConsensus, PathBcftools);
	processCombFasta(ARG_G, ARG_L, DirConsensus);
    }

    if (ARG_F == "all" || ARG_F == "align") {
	processAlign(ARG_G, DirConsensus, DirAlign, PathMacse);
    }

    if (ARG_F == "all" || ARG_F == "trim") {
	processTrimming(ARG_G, DirAlign, DirTrim, PathDelstop, PathTrimal);
    }

    writeln("RGBEPP::End");
}

