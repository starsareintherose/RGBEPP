#!/usr/bin/env rdmd

import std.stdio;
import std.file;
import std.process;
import std.algorithm;
import std.conv;
import std.array;
import std.path;

void show_help(string pkgver) {
    writeln("\t\t\t\t\t\033[0;47;31mR\033[0m\033[0;47;92mG\033[0m\033[0;47;94mB\033[0m\033[0;47m \033[0m\033[0;47;33mE\033[0m\033[0;47;94mP\033[0m\033[0;47;33mP\033[0m
\t\t\tReference Genome based Exon Phylogeny Pipeline
	    Version: ", pkgver, "
	    License: GPL-2.0-only
	    Author: Guoyi Zhang
	    -c\t--contigs\tcontigs type: scaffolds or contigs
	    -g\t--genes\t\tgene file path
	    -f\t--functions\tfunctions type (optional): all clean map postmap varcall consen
	    -h\t--help\t\tshow this information
	    -l\t--list\t\tlist file path
	    -m\t--memory\tmemory settings (optional, default 16 GB)
	    -r\t--reference\treference genome path
	    -t\t--threads\tthreads setting (optional, default 8 threads)
	    --macse\t\tMacse jarfile path
	    for example: ./your_program -c scaffolds -f all -l list -g genes \\ 
	    -r reference.fasta \n");
}

void createDir(string path) {
    if (!exists(path)) {
        mkdir(path);
    }
}

void executeCommand(string[] cmd, string workingDir = "") {
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

void performQualityControl(string[] ARG_L, string DirRaw, string DirQcTrim, int ARG_T) {
    // Prepare directory
    createDir(DirQcTrim);

    foreach (string file; ARG_L) {
        string baseName = getBaseName(file);
        string inputFileR1 = DirRaw ~ "/" ~ baseName ~ "_R1.fastq.gz";
        string inputFileR2 = DirRaw ~ "/" ~ baseName ~ "_R2.fastq.gz";
        string outputFileR1 = DirQcTrim ~ "/" ~ baseName ~ "_R1.fastq.gz";
        string outputFileR2 = DirQcTrim ~ "/" ~ baseName ~ "_R2.fastq.gz";
        string jsonFile = DirQcTrim ~ "/" ~ baseName ~ ".json";
        string htmlFile = DirQcTrim ~ "/" ~ baseName ~ ".html";

        // Perform quality control and trimming using external program `fastp`
        string[] cmd = ["fastp", "-i", inputFileR1, "-I", inputFileR2,
                        "-o", outputFileR1, "-O", outputFileR2,
                        "-j", jsonFile, "-h", htmlFile,
                        "-w", ARG_T.to!string];
        executeCommand(cmd);
    }

}

void performMapping(string ARG_R, string[] ARG_L, string DirQcTrim, string DirMap, int ARG_T) {
    writeln("Mapping::Start");

    // Prepare directory
    createDir(DirMap);

    createDir(DirMap ~ "/index");

    string[] Refs = getRef(ARG_R, DirMap);
    string ARG_R_index = Refs[0]; // bt2_index_base
    string ARG_R_refer = Refs[1]; //reference_in fasta file

    copy(ARG_R, ARG_R_refer);

    string[] cmdBuildDB = ["bowtie2-build", "--threads", ARG_T.to!string, ARG_R_refer,  ARG_R_index];
    executeCommand(cmdBuildDB);

    foreach (string file; ARG_L) {
        string baseName = baseName(file, ".fastq.gz");
        string outputBam = DirMap ~ "/" ~ baseName ~ ".bam";
        string inputFileR1 = DirQcTrim ~ "/" ~ baseName ~ "_R1.fastq.gz";
        string inputFileR2 = DirQcTrim ~ "/" ~ baseName ~ "_R2.fastq.gz";

        // Perform mapping using Bowtie2 and converted to Bam using samtools
        string[] cmdMap = ["bowtie2", "-x", ARG_R_index, "-1", inputFileR1, "-2", inputFileR2,
                        "-p", ARG_T.to!string];
        string[] cmdSam2Bam = ["samtools", "view", "-bS", "-@", ARG_T.to!string, "-o", outputBam];

	executeCommandPipe([cmdMap, cmdSam2Bam]);

    }
    writeln("Mapping::End");
}

void processPostMap(string[] ARG_L, int ARG_T, string DirMap, string DirBam) {

    createDir(DirBam);
    writeln("PostMapping::Start");

    foreach (string file; ARG_L) {
        string baseName = getBaseName(file);
        string inputBam = DirMap ~ "/" ~ baseName ~ ".bam";
        string outputBam = DirBam ~ "/" ~ baseName ~ ".bam";

        // Convert SAM to BAM, sort and remove duplicates using Samtools
	string[] cmdFixmate = ["samtools", "fixmate", "-@", ARG_T.to!string, "-m", inputBam, "-"];
        string[] cmdSort = ["samtools", "sort", "-@", ARG_T.to!string, "-"];
        string[] cmdMarkdup = ["samtools", "markdup", "-@", ARG_T.to!string, "-", outputBam];
	executeCommandPipe([cmdFixmate, cmdSort, cmdMarkdup]);

        string [] cmdIndexBam = ["samtools", "index", "-@", ARG_T.to!string, outputBam];
        executeCommand(cmdIndexBam);
    }

    writeln("PostMapping::End");
}


void processVarCall(string[] ARG_L, string ARG_R, int ARG_T, string DirBam, string DirVcf, string DirMap) {
    writeln("VarCalling::Start");

    string[] Refs = getRef(ARG_R, DirMap);
    string ARG_R_refer = Refs[1]; //reference_in fasta file

    createDir(DirVcf);

    foreach (string file; ARG_L) {
        string baseName = getBaseName(file);
        string inputBam = DirBam ~ "/" ~ baseName ~ ".bam";
        string outputVcf = DirVcf ~ "/" ~ baseName ~ ".vcf.gz";

        // Variant calling using bcftools
        string[] cmdPileup = ["bcftools", "mpileup", "-Oz", "--threads", ARG_T.to!string, "-f", ARG_R_refer, inputBam];
	string[] cmdVarCall = ["bcftools", "call", "-mv", "-Oz", "--threads", ARG_T.to!string];
	string[] cmdNorm = ["bcftools", "norm", "--threads", ARG_T.to!string, "-f", ARG_R_refer, "-Oz"];
	string[] cmdFilter = ["bcftools", "filter", "--threads", ARG_T.to!string, "--IndelGap", "5", "-Oz", "-o", outputVcf];
        executeCommandPipe([cmdPileup, cmdVarCall, cmdNorm, cmdFilter]); 
    }

    writeln("VarCalling::End");

}

void processCon(string[] ARG_L, string ARG_R, int ARG_T, string[] ARG_G, string DirVcf, string DirConsensus, string DirMap) {
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
	string[] cmdIndexVcf = ["bcftools", "index", inputVcf];
	executeCommand(cmdIndexVcf);
	
        // Generate consensus sequences using bcftools
        string[] cmdCon = ["bcftools", "consensus", "-f", ARG_R, inputVcf, "-o", outputFasta];
	executeCommand(cmdCon);
    }
    // Recombine the sequences based on genes
    writeln("Consensus::End");
}

void processCombFasta(string[] ARG_L, string[] ARG_G, string DirConsensus) {

    string DirConTaxa = DirConsensus ~ "/" ~ "taxa";
    string DirConGene = DirConsensus ~ "/" ~ "gene";    
    
    // create a dictory
    string[string] geneSequences;

    // read first
    foreach (file; ARG_L) {
        string inputFile = DirConTaxa ~ "/" ~ file ~ ".fas";
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
        File output = File(outputFile, "a");
        if (gene in geneSequences) {
            output.write(geneSequences[gene]);
        }
    }

}

void main(string[] args) {
    string pkgver = "0.0.3";

    string DirHome = std.file.getcwd();
    string DirRaw = DirHome ~ "/00_raw";
    string DirQcTrim = DirHome ~ "/01_fastp";
    string DirMap = DirHome ~ "/02_bowtie2";
    string DirBam = DirHome ~ "/03_bam";
    string DirVcf = DirHome ~ "/04_vcf";
    string DirConsensus = DirHome ~ "/05_consen";
    string DirAlign = DirHome ~ "/06_macse";

    string PathMacse = "/usr/share/java/macse.jar";

    int ARG_T = 8;
    string[] ARG_G;
    string[] ARG_L;
    string ARG_F;
    string ARG_R;
   
    if (args.length > 1){
        foreach (int i; 0 .. cast(int)args.length) {
            switch (args[i]) {
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
                case "--macse":
		    i++;
                    PathMacse = args[i];
                    break;
                default:
                    break;
            }
        }
    } else {
        show_help(pkgver);
        return;
    }

    writeln("RGBEPP::Start");
    // Perform steps based on provided function argument
    if (ARG_F == "all" || ARG_F == "clean") {
        performQualityControl(ARG_L, DirRaw, DirQcTrim, ARG_T);
    }

    if (ARG_F == "all" || ARG_F == "map") {
        performMapping(ARG_R, ARG_L, DirQcTrim, DirMap, ARG_T);
    }

    if (ARG_F == "all" || ARG_F == "postmap") {
        processPostMap(ARG_L, ARG_T, DirMap, DirBam);
    }

    if (ARG_F == "all" || ARG_F == "varcall") {
        processVarCall(ARG_L, ARG_R, ARG_T, DirBam, DirVcf, DirMap);
    }

    if (ARG_F == "all" || ARG_F == "consen") {
        processCon(ARG_L, ARG_R, ARG_T, ARG_G, DirVcf, DirConsensus,DirMap);
	processCombFasta(ARG_L, ARG_G, DirConsensus);
    }

    if (ARG_F == "all" || ARG_F == "align") {

    }

    writeln("RGBEPP::End");
}

