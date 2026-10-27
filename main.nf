nextflow.enable.dsl=2

params.help = false
params.rgbepp = params.rgbepp ?: 'RGBEPP'
params.list = params.list ?: 'list'
params.reference = params.reference ?: null
params.genes = params.genes ?: null
params.config = params.config ?: null
params.outdir = params.outdir ?: 'results'
params.threads = params.threads ?: 8
params.memory = params.memory ?: 16
params.codon = params.codon ?: false
params.function = params.function ?: 'all'
params.log = params.log ?: 'rgbepp.log'
params.publish_mode = params.publish_mode ?: 'copy'
params.fastp = params.fastp ?: 'fastp'
params.spades = params.spades ?: 'spades.py'
params.diamond = params.diamond ?: 'diamond'
params.sortdiamond = params.sortdiamond ?: 'sortdiamond'
params.bowtie2 = params.bowtie2 ?: 'bowtie2'
params.samtools = params.samtools ?: 'samtools'
params.bcftools = params.bcftools ?: 'bcftools'
params.exonerate = params.exonerate ?: 'exonerate'
params.macse = params.macse ?: '/usr/share/java/macse.jar'
params.delstop = params.delstop ?: 'delstop'
params.deltaxa = params.deltaxa ?: 'deltaxa'
params.trimal = params.trimal ?: 'trimal'
params.concataln = params.concataln ?: 'concataln'
params.raw_dir = params.raw_dir ?: '00_raw'
params.fastp_dir = params.fastp_dir ?: '01_fastp'
params.spades_dir = params.spades_dir ?: '02_spades'
params.bowtie2_dir = params.bowtie2_dir ?: '03_bowtie2'
params.bam_dir = params.bam_dir ?: '04_bam'
params.vcf_dir = params.vcf_dir ?: '05_vcf'
params.consen_dir = params.consen_dir ?: '06_consen'
params.macse_dir = params.macse_dir ?: '07_macse'
params.trimal_dir = params.trimal_dir ?: '08_trimal'
params.fastp_paras = params.fastp_paras ?: ''
params.spades_paras = params.spades_paras ?: '--careful --phred-offset 33'
params.diamond_makedb_paras = params.diamond_makedb_paras ?: ''
params.diamond_blastx_paras = params.diamond_blastx_paras ?: '--ultra-sensitive'
params.bowtie2_build_paras = params.bowtie2_build_paras ?: ''
params.bowtie2_paras = params.bowtie2_paras ?: ''
params.samtools_view_paras = params.samtools_view_paras ?: ''
params.samtools_fixmate_paras = params.samtools_fixmate_paras ?: ''
params.samtools_sort_paras = params.samtools_sort_paras ?: ''
params.samtools_markdup_paras = params.samtools_markdup_paras ?: ''
params.samtools_index_paras = params.samtools_index_paras ?: ''
params.bcftools_mpileup_paras = params.bcftools_mpileup_paras ?: ''
params.bcftools_call_paras = params.bcftools_call_paras ?: '-mv'
params.bcftools_norm_paras = params.bcftools_norm_paras ?: '--check-ref s'
params.bcftools_filter_paras = params.bcftools_filter_paras ?: ''
params.bcftools_index_paras = params.bcftools_index_paras ?: ''
params.bcftools_consensus_paras = params.bcftools_consensus_paras ?: ''
params.exonerate_paras = params.exonerate_paras ?: ''
params.macse_paras = params.macse_paras ?: ''
params.trimal_paras = params.trimal_paras ?: '-gt 0.7'

def helpMessage() {
    return """
    RGBEPP Nextflow workflow

    General options:
      --codon         Run the optional codon extraction step before alignment
      --config        RGBEPP config; values override conflicting tool/directory options
      --function      Run through this RGBEPP function (default: all)
      --genes         Gene list file; generated from --reference if omitted
      --list          Sample list file (default: list)
      --log           Log filename/prefix (default: rgbepp.log; one file per stage)
      --memory        Memory in GB for RGBEPP/SPAdes (default: 16)
      --outdir        Published results directory (default: results)
      --reference     Reference amino-acid FASTA file (required)
      --rgbepp        RGBEPP executable command (default: RGBEPP)
      --threads       Threads per RGBEPP step (default: 8)

    Tool paths (command names are resolved through PATH):
      --bcftools      bcftools path
      --bowtie2       bowtie2 path
      --concataln     concataln path
      --delstop       delstop path
      --deltaxa       deltaxa path
      --diamond       diamond path
      --exonerate     exonerate path
      --fastp         fastp path
      --macse         MACSE JAR path
      --samtools      samtools path
      --sortdiamond   sortdiamond path
      --spades        spades.py path
      --trimal        trimal path

    Directories (config.example):
      --raw_dir       Directory containing paired FASTQ files (default: 00_raw)
      --fastp_dir     RGBEPP fastp output directory (default: 01_fastp)
      --spades_dir    RGBEPP SPAdes output directory (default: 02_spades)
      --bowtie2_dir   RGBEPP Bowtie2 output directory (default: 03_bowtie2)
      --bam_dir       RGBEPP BAM output directory (default: 04_bam)
      --vcf_dir       RGBEPP VCF output directory (default: 05_vcf)
      --consen_dir    RGBEPP consensus output directory (default: 06_consen)
      --macse_dir     RGBEPP MACSE output directory (default: 07_macse)
      --trimal_dir    RGBEPP trim output directory (default: 08_trimal)
      Tool-specific *_paras values are configured in the config file.

    Example:
      nextflow run starsareintherose/RGBEPP --reference reference.aa.fasta --list list --raw_dir 00_raw --rgbepp RGBEPP
    """.stripIndent()
}

def cmdPath(value) {
    return value && value.contains('/') ? file(value).toAbsolutePath() : value
}

def shellQuote(value) {
    return "'${value.toString().replace("'", "'\"'\"'")}'"
}

def logName(stage) {
    def requested = params.log.toString()
    if (!requested || requested.contains('/') || requested.contains('\\')) {
        error "--log must be a filename or prefix, not a path: ${requested}"
    }
    def dot = requested.lastIndexOf('.')
    return dot > 0
        ? "${requested.substring(0, dot)}.${stage}${requested.substring(dot)}"
        : "${requested}.${stage}.log"
}

def logArg(stage) {
    return "--log ${shellQuote(logName(stage))}"
}

def functionStages() {
    return ['clean', 'assembly', 'map', 'postmap', 'varcall', 'consen', 'codon', 'align', 'ortholog', 'trim', 'concat']
}

def functionEnabled(stage) {
    return params.function == 'all' || functionStages().indexOf(stage) <= functionStages().indexOf(params.function)
}

def codonEnabled() {
    return params.codon || params.function == 'codon'
}

def validateFunction() {
    def valid = functionStages() + ['all']
    if (!valid.contains(params.function)) {
        error "Invalid --function '${params.function}'. Valid values: ${valid.join(', ')}"
    }
}

def validateRelativeDirs(directoryValues) {
    def dirParams = ['raw_dir', 'fastp_dir', 'spades_dir', 'bowtie2_dir', 'bam_dir', 'vcf_dir', 'consen_dir', 'macse_dir', 'trimal_dir']
    dirParams.each { name ->
        def value = directoryValues[name]?.toString()
        if (!value) {
            error "Parameter --${name} cannot be empty"
        }
        if (value.startsWith('/')) {
            error "Parameter --${name} must be a relative path inside the Nextflow work directory: ${value}"
        }
        if (value == '..' || value.startsWith('../') || value.endsWith('/..') || value.contains('/../')) {
            error "Parameter --${name} must not contain '..': ${value}"
        }
    }
}

def rgbeppCmd() {
    return params.rgbepp.contains('/') ? file(params.rgbepp).toAbsolutePath() : params.rgbepp
}

def suppliedConfig() {
    return params.config ? file(params.config).toAbsolutePath() : null
}

def externalConfigValues() {
    def values = [:]
    def configFile = suppliedConfig()
    if (!configFile) return values
    if (!configFile.isFile()) error "RGBEPP config file does not exist: ${configFile}"
    configFile.eachLine { line ->
        def text = line.trim()
        if (!text || text.startsWith('#')) return
        int separator = text.indexOf('=')
        if (separator > 0) {
            values[text.substring(0, separator).trim()] = text.substring(separator + 1).trim()
        }
    }
    return values
}

def applyExternalConfigOverrides() {
    if (!params.config) return
    def values = externalConfigValues()
    def toolDefaults = [
        fastp: 'fastp', spades: 'spades.py', diamond: 'diamond', sortdiamond: 'sortdiamond',
        bowtie2: 'bowtie2', samtools: 'samtools', bcftools: 'bcftools', exonerate: 'exonerate',
        macse: '/usr/share/java/macse.jar', delstop: 'delstop', deltaxa: 'deltaxa',
        trimal: 'trimal', concataln: 'concataln'
    ]
    toolDefaults.each { key, defaultValue ->
        def configValue = values[key]
        def parameterValue = params[key]?.toString()
        if (configValue && parameterValue && parameterValue != defaultValue && parameterValue != configValue) {
            log.warn "--${key}='${parameterValue}' conflicts with config value '${configValue}'; the config value takes precedence"
        }
    }

}

def effectiveDirectoryValues() {
    def defaults = [
        raw_dir: '00_raw', fastp_dir: '01_fastp', spades_dir: '02_spades',
        bowtie2_dir: '03_bowtie2', bam_dir: '04_bam', vcf_dir: '05_vcf',
        consen_dir: '06_consen', macse_dir: '07_macse', trimal_dir: '08_trimal'
    ]
    def configValues = params.config ? externalConfigValues() : [:]
    def result = [:]
    defaults.each { key, defaultValue ->
        def parameterValue = params[key]?.toString() ?: defaultValue
        def configValue = configValues[key]
        if (configValue) {
            if (parameterValue != defaultValue && parameterValue != configValue) {
                log.warn "--${key}='${parameterValue}' conflicts with config value '${configValue}'; the config value takes precedence"
            }
            result[key] = configValue
        } else {
            result[key] = parameterValue
        }
    }
    return result
}

def configArg() {
    return suppliedConfig() ? "-c ${suppliedConfig()}" : '-c rgbepp.nextflow.config'
}

def commonArgs() {
    return "-t ${params.threads} -m ${params.memory} ${configArg()}".trim()
}

def toolConfigLine(key, value) {
    return "write_tool ${shellQuote(key)} ${shellQuote(cmdPath(value))}"
}

def valueConfigLine(key, value) {
    return "write_value ${shellQuote(key)} ${shellQuote(value)}"
}

def generatedConfigCommands() {
    return [
        toolConfigLine('fastp', params.fastp),
        toolConfigLine('spades', params.spades),
        toolConfigLine('diamond', params.diamond),
        toolConfigLine('sortdiamond', params.sortdiamond),
        toolConfigLine('bowtie2', params.bowtie2),
        toolConfigLine('samtools', params.samtools),
        toolConfigLine('bcftools', params.bcftools),
        toolConfigLine('exonerate', params.exonerate),
        toolConfigLine('macse', params.macse),
        toolConfigLine('delstop', params.delstop),
        toolConfigLine('deltaxa', params.deltaxa),
        toolConfigLine('trimal', params.trimal),
        toolConfigLine('concataln', params.concataln),
        valueConfigLine('raw_dir', params.raw_dir),
        valueConfigLine('fastp_dir', params.fastp_dir),
        valueConfigLine('spades_dir', params.spades_dir),
        valueConfigLine('bowtie2_dir', params.bowtie2_dir),
        valueConfigLine('bam_dir', params.bam_dir),
        valueConfigLine('vcf_dir', params.vcf_dir),
        valueConfigLine('consen_dir', params.consen_dir),
        valueConfigLine('macse_dir', params.macse_dir),
        valueConfigLine('trimal_dir', params.trimal_dir),
        valueConfigLine('fastp_paras', params.fastp_paras),
        valueConfigLine('spades_paras', params.spades_paras),
        valueConfigLine('diamond_makedb_paras', params.diamond_makedb_paras),
        valueConfigLine('diamond_blastx_paras', params.diamond_blastx_paras),
        valueConfigLine('bowtie2_build_paras', params.bowtie2_build_paras),
        valueConfigLine('bowtie2_paras', params.bowtie2_paras),
        valueConfigLine('samtools_view_paras', params.samtools_view_paras),
        valueConfigLine('samtools_fixmate_paras', params.samtools_fixmate_paras),
        valueConfigLine('samtools_sort_paras', params.samtools_sort_paras),
        valueConfigLine('samtools_markdup_paras', params.samtools_markdup_paras),
        valueConfigLine('samtools_index_paras', params.samtools_index_paras),
        valueConfigLine('bcftools_mpileup_paras', params.bcftools_mpileup_paras),
        valueConfigLine('bcftools_call_paras', params.bcftools_call_paras),
        valueConfigLine('bcftools_norm_paras', params.bcftools_norm_paras),
        valueConfigLine('bcftools_filter_paras', params.bcftools_filter_paras),
        valueConfigLine('bcftools_index_paras', params.bcftools_index_paras),
        valueConfigLine('bcftools_consensus_paras', params.bcftools_consensus_paras),
        valueConfigLine('exonerate_paras', params.exonerate_paras),
        valueConfigLine('macse_paras', params.macse_paras),
        valueConfigLine('trimal_paras', params.trimal_paras)
    ].join('\n    ')
}

def configSetup() {
    return suppliedConfig() ? '' : """
resolve_tool() {
    if [[ "\$1" == */* ]]; then
        [[ -e "\$1" ]] || { echo "Tool path does not exist: \$1" >&2; return 1; }
        printf '%s' "\$1"
    else
        command -v "\$1" || { echo "Tool not found in PATH: \$1" >&2; return 1; }
    fi
}

write_tool() {
    local key="\$1"
    local value
    value="\$(resolve_tool "\$2")"
    printf '%s = %s\n' "\$key" "\$value"
}

write_value() {
    printf '%s = %s\n' "\$1" "\$2"
}

{
    ${generatedConfigCommands()}
} > rgbepp.nextflow.config
""".stripIndent().trim()
}

process PREPARE_GENES {
    tag 'prepare_genes'
    publishDir params.outdir, mode: params.publish_mode, overwrite: true

    input:
    path reference

    output:
    path 'genes'

    script:
    """
    set -euo pipefail
    awk '/^>/ { sub(/^>/, ""); print }' ${reference} > genes
    """
}

process CLEAN {
    tag 'clean'
    cpus params.threads
    memory "${params.memory} GB"
    publishDir params.outdir, mode: params.publish_mode, overwrite: true

    input:
    path sample_list, name: 'list'
    path raw_reads

    val directories

    output:
    path "${directories.fastp_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    if [[ "${directories.raw_dir}" != /* ]]; then
        mkdir -p "${directories.raw_dir}"
        cp -L ${raw_reads.join(' ')} "${directories.raw_dir}/"
    fi
    ${rgbeppCmd()} -f clean -l list ${commonArgs()} ${logArg('clean')}
    cp ${shellQuote(logName('clean'))} "${directories.fastp_dir}/"
    compgen -G "${directories.fastp_dir}/*_R1.fastq.gz" > /dev/null || { echo "No R1 cleaned reads were produced in ${directories.fastp_dir}" >&2; exit 1; }
    compgen -G "${directories.fastp_dir}/*_R2.fastq.gz" > /dev/null || { echo "No R2 cleaned reads were produced in ${directories.fastp_dir}" >&2; exit 1; }
    """
}

process ASSEMBLY {
    tag 'assembly'
    cpus params.threads
    memory "${params.memory} GB"
    publishDir params.outdir, mode: params.publish_mode, overwrite: true, enabled: params.function == 'assembly'

    input:
    path sample_list, name: 'list'
    path fastp_dir, name: 'input_01_fastp'

    val directories

    output:
    path "${directories.spades_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_01_fastp "${directories.fastp_dir}"
    ${rgbeppCmd()} -f assembly -l list ${commonArgs()} ${logArg('assembly')}
    cp ${shellQuote(logName('assembly'))} "${directories.spades_dir}/"
    compgen -G "${directories.spades_dir}/scaffolds/*.fasta" > /dev/null || { echo "No scaffold FASTA files were produced in ${directories.spades_dir}/scaffolds" >&2; exit 1; }
    """
}

process MAP {
    tag 'map'
    cpus params.threads
    memory "${params.memory} GB"
    publishDir params.outdir, mode: params.publish_mode, overwrite: true

    input:
    path sample_list, name: 'list'
    path reference
    path fastp_dir, name: 'input_01_fastp'
    path assembly_dir, name: 'input_02_spades'

    val directories

    output:
    path "${directories.spades_dir}", emit: assembly
    path "${directories.bowtie2_dir}", emit: map

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_01_fastp "${directories.fastp_dir}"
    cp -aL input_02_spades "${directories.spades_dir}"
    ${rgbeppCmd()} -f map -l list -r ${reference} ${commonArgs()} ${logArg('map')}
    cp ${shellQuote(logName('map'))} "${directories.bowtie2_dir}/"
    compgen -G "${directories.spades_dir}/fasta/*.fasta" > /dev/null || { echo "No mapped reference FASTA files were produced in ${directories.spades_dir}/fasta" >&2; exit 1; }
    compgen -G "${directories.bowtie2_dir}/*.bam" > /dev/null || { echo "No mapping BAM files were produced in ${directories.bowtie2_dir}" >&2; exit 1; }
    """
}

process POSTMAP {
    tag 'postmap'
    cpus params.threads
    memory "${params.memory} GB"
    publishDir params.outdir, mode: params.publish_mode, overwrite: true

    input:
    path sample_list, name: 'list'
    path map_dir, name: 'input_03_bowtie2'

    val directories

    output:
    path "${directories.bam_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_03_bowtie2 "${directories.bowtie2_dir}"
    ${rgbeppCmd()} -f postmap -l list ${commonArgs()} ${logArg('postmap')}
    cp ${shellQuote(logName('postmap'))} "${directories.bam_dir}/"
    compgen -G "${directories.bam_dir}/*.bam" > /dev/null || { echo "No post-mapping BAM files were produced in ${directories.bam_dir}" >&2; exit 1; }
    compgen -G "${directories.bam_dir}/*.bam.bai" > /dev/null || { echo "No BAM index files were produced in ${directories.bam_dir}" >&2; exit 1; }
    """
}

process VARCALL {
    tag 'varcall'
    cpus params.threads
    memory "${params.memory} GB"
    publishDir params.outdir, mode: params.publish_mode, overwrite: true

    input:
    path sample_list, name: 'list'
    path assembly_dir, name: 'input_02_spades'
    path bam_dir, name: 'input_04_bam'

    val directories

    output:
    path "${directories.vcf_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_02_spades "${directories.spades_dir}"
    cp -aL input_04_bam "${directories.bam_dir}"
    ${rgbeppCmd()} -f varcall -l list ${commonArgs()} ${logArg('varcall')}
    cp ${shellQuote(logName('varcall'))} "${directories.vcf_dir}/"
    compgen -G "${directories.vcf_dir}/*.vcf.gz" > /dev/null || { echo "No VCF files were produced in ${directories.vcf_dir}" >&2; exit 1; }
    """
}

process CONSEN {
    tag 'consen'
    cpus params.threads
    memory "${params.memory} GB"
    publishDir params.outdir, mode: params.publish_mode, overwrite: true, enabled: !(params.codon || params.function == 'codon')

    input:
    path sample_list, name: 'list'
    path genes, name: 'genes'
    path assembly_dir, name: 'input_02_spades'
    path vcf_dir, name: 'input_05_vcf'

    val directories

    output:
    path "${directories.consen_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_02_spades "${directories.spades_dir}"
    cp -aL input_05_vcf "${directories.vcf_dir}"
    ${rgbeppCmd()} -f consen -l list -g genes ${commonArgs()} ${logArg('consen')}
    cp ${shellQuote(logName('consen'))} "${directories.consen_dir}/"
    compgen -G "${directories.consen_dir}/taxa/*.fasta" > /dev/null || { echo "No consensus taxa FASTA files were produced in ${directories.consen_dir}/taxa" >&2; exit 1; }
    compgen -G "${directories.consen_dir}/gene/*.fasta" > /dev/null || { echo "No consensus gene FASTA files were produced in ${directories.consen_dir}/gene" >&2; exit 1; }
    """
}

process CODON {
    tag 'codon'
    cpus params.threads
    memory "${params.memory} GB"
    publishDir params.outdir, mode: params.publish_mode, overwrite: true

    input:
    path genes, name: 'genes'
    path reference
    path consensus_dir, name: 'input_06_consen'

    val directories

    output:
    path "${directories.consen_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_06_consen "${directories.consen_dir}"
    ${rgbeppCmd()} -f codon -g genes -r ${reference} ${commonArgs()} ${logArg('codon')}
    cp ${shellQuote(logName('codon'))} "${directories.consen_dir}/"
    compgen -G "${directories.consen_dir}/gene/*.fasta" > /dev/null || { echo "No codon-filtered gene FASTA files were produced in ${directories.consen_dir}/gene" >&2; exit 1; }
    """
}

process ALIGN {
    tag 'align'
    cpus params.threads
    memory "${params.memory} GB"
    publishDir params.outdir, mode: params.publish_mode, overwrite: true, enabled: params.function == 'align'

    input:
    path genes, name: 'genes'
    path consensus_dir, name: 'input_06_consen'

    val directories

    output:
    path "${directories.macse_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_06_consen "${directories.consen_dir}"
    ${rgbeppCmd()} -f align -g genes ${commonArgs()} ${logArg('align')}
    cp ${shellQuote(logName('align'))} "${directories.macse_dir}/"
    compgen -G "${directories.macse_dir}/AA/*.fasta" > /dev/null || { echo "No MACSE AA alignments were produced in ${directories.macse_dir}/AA" >&2; exit 1; }
    compgen -G "${directories.macse_dir}/NT/*.fasta" > /dev/null || { echo "No MACSE NT alignments were produced in ${directories.macse_dir}/NT" >&2; exit 1; }
    """
}

process ORTHOLOG {
    tag 'ortholog'
    cpus params.threads
    memory "${params.memory} GB"
    publishDir params.outdir, mode: params.publish_mode, overwrite: true

    input:
    path sample_list, name: 'list'
    path genes, name: 'genes'
    path reference
    path align_dir, name: 'input_07_macse'

    val directories

    output:
    path "${directories.macse_dir}", emit: align
    path 'paralog.csv', emit: paralog

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_07_macse "${directories.macse_dir}"
    ${rgbeppCmd()} -f ortholog -l list -g genes -r ${reference} ${commonArgs()} ${logArg('ortholog')}
    cp ${shellQuote(logName('ortholog'))} "${directories.macse_dir}/"
    [[ -e paralog.csv ]] || { echo "Paralog report was not produced: paralog.csv" >&2; exit 1; }
    compgen -G "${directories.macse_dir}/taxa/*.tsv" > /dev/null || { echo "No ortholog TSV files were produced in ${directories.macse_dir}/taxa" >&2; exit 1; }
    """
}

process TRIM {
    tag 'trim'
    cpus params.threads
    memory "${params.memory} GB"
    publishDir params.outdir, mode: params.publish_mode, overwrite: true, enabled: params.function == 'trim'

    input:
    path genes, name: 'genes'
    path align_dir, name: 'input_07_macse'

    val directories

    output:
    path "${directories.trimal_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_07_macse "${directories.macse_dir}"
    ${rgbeppCmd()} -f trim -g genes ${commonArgs()} ${logArg('trim')}
    cp ${shellQuote(logName('trim'))} "${directories.trimal_dir}/"
    compgen -G "${directories.trimal_dir}/NT/*.fasta" > /dev/null || { echo "No trimmed FASTA files were produced in ${directories.trimal_dir}/NT" >&2; exit 1; }
    """
}

process CONCAT {
    tag 'concat'
    cpus params.threads
    memory "${params.memory} GB"
    publishDir params.outdir, mode: params.publish_mode, overwrite: true

    input:
    path genes, name: 'genes'
    path trim_dir, name: 'input_08_trimal'

    val directories

    output:
    path "${directories.trimal_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_08_trimal "${directories.trimal_dir}"
    ${rgbeppCmd()} -f concat -g genes ${commonArgs()} ${logArg('concat')}
    cp ${shellQuote(logName('concat'))} "${directories.trimal_dir}/"
    [[ -s "${directories.trimal_dir}/concat.fasta.fasta" ]] || { echo "Concatenated FASTA was not produced: ${directories.trimal_dir}/concat.fasta.fasta" >&2; exit 1; }
    """
}

workflow {
    if (params.help) {
        log.info helpMessage()
        System.exit(0)
    }

    if (!params.reference) {
        error 'Missing required parameter: --reference. Run `nextflow run starsareintherose/RGBEPP --help` for usage.'
    }

    applyExternalConfigOverrides()
    directory_values = effectiveDirectoryValues()
    validateFunction()
    validateRelativeDirs(directory_values)
    directory_values_ch = Channel.value(directory_values)

    sample_list_ch = Channel.fromPath(params.list, checkIfExists: true)
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
    raw_reads_ch = Channel.fromPath("${directory_values.raw_dir}/*_R[12].fastq.gz", checkIfExists: true).collect()

    if (params.genes) {
        genes_ch = Channel.fromPath(params.genes, checkIfExists: true)
    } else {
        PREPARE_GENES(reference_ch)
        genes_ch = PREPARE_GENES.out
    }

    CLEAN(sample_list_ch, raw_reads_ch, directory_values_ch)

    if (functionEnabled('assembly')) {
        ASSEMBLY(sample_list_ch, CLEAN.out, directory_values_ch)
    }

    if (functionEnabled('map')) {
        MAP(sample_list_ch, reference_ch, CLEAN.out, ASSEMBLY.out, directory_values_ch)
    }

    if (functionEnabled('postmap')) {
        POSTMAP(sample_list_ch, MAP.out.map, directory_values_ch)
    }

    if (functionEnabled('varcall')) {
        VARCALL(sample_list_ch, MAP.out.assembly, POSTMAP.out, directory_values_ch)
    }

    if (functionEnabled('consen')) {
        CONSEN(sample_list_ch, genes_ch, MAP.out.assembly, VARCALL.out, directory_values_ch)
    }

    if (functionEnabled('codon') && codonEnabled()) {
        CODON(genes_ch, reference_ch, CONSEN.out, directory_values_ch)
    }

    if (functionEnabled('align')) {
        if (codonEnabled()) {
            ALIGN(genes_ch, CODON.out, directory_values_ch)
        } else {
            ALIGN(genes_ch, CONSEN.out, directory_values_ch)
        }
    }

    if (functionEnabled('ortholog')) {
        ORTHOLOG(sample_list_ch, genes_ch, reference_ch, ALIGN.out, directory_values_ch)
    }

    if (functionEnabled('trim')) {
        TRIM(genes_ch, ORTHOLOG.out.align, directory_values_ch)
    }

    if (functionEnabled('concat')) {
        CONCAT(genes_ch, TRIM.out, directory_values_ch)
    }
}
