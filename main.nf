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

    Required:
      --reference     Reference amino-acid FASTA file

    Common options:
      --list          Sample list file (default: list)
      --genes         Gene list file; generated from --reference if omitted
      --raw_dir       Directory containing *_R1.fastq.gz and *_R2.fastq.gz files (default: 00_raw)
      --rgbepp        RGBEPP executable command (default: RGBEPP)
      --config        RGBEPP config file for tool paths and parameters
      --threads       Threads per RGBEPP step (default: 8)
      --memory        Memory in GB for RGBEPP/SPAdes (default: 16)
      --codon         Run the optional codon extraction step before alignment
      --function      Run through this RGBEPP function (default: all)
      --outdir        Published results directory (default: results)

    RGBEPP config keys:
      --fastp         fastp executable path
      --spades        spades.py executable path
      --diamond       diamond executable path
      --sortdiamond   sortdiamond executable path
      --bowtie2       bowtie2 executable path
      --samtools      samtools executable path
      --bcftools      bcftools executable path
      --exonerate     exonerate executable path
      --macse         macse jar path
      --delstop       delstop executable path
      --deltaxa       deltaxa executable path
      --trimal        trimal executable path
      --concataln     concataln executable path
      --fastp_dir     RGBEPP fastp output directory (default: 01_fastp)
      --spades_dir    RGBEPP spades output directory (default: 02_spades)
      --bowtie2_dir   RGBEPP bowtie2 output directory (default: 03_bowtie2)
      --bam_dir       RGBEPP BAM output directory (default: 04_bam)
      --vcf_dir       RGBEPP VCF output directory (default: 05_vcf)
      --consen_dir    RGBEPP consensus output directory (default: 06_consen)
      --macse_dir     RGBEPP MACSE output directory (default: 07_macse)
      --trimal_dir    RGBEPP trim output directory (default: 08_trimal)
      --*_paras       Tool-specific RGBEPP config parameters, e.g. --spades_paras "--careful --phred-offset 33"

    Example:
      nextflow run . --reference reference.aa.fasta --list list --raw_dir 00_raw --rgbepp RGBEPP
    """.stripIndent()
}

def cmdPath(value) {
    return value && value.contains('/') ? file(value).toAbsolutePath() : value
}

def shellQuote(value) {
    return "'${value.toString().replace("'", "'\"'\"'")}'"
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

def validateRelativeDirs() {
    def dirParams = ['raw_dir', 'fastp_dir', 'spades_dir', 'bowtie2_dir', 'bam_dir', 'vcf_dir', 'consen_dir', 'macse_dir', 'trimal_dir']
    dirParams.each { name ->
        def value = params[name]?.toString()
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

    output:
    path "${params.fastp_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    if [[ "${params.raw_dir}" != /* ]]; then
        mkdir -p "${params.raw_dir}"
        cp -L ${raw_reads.join(' ')} "${params.raw_dir}/"
    fi
    ${rgbeppCmd()} -f clean -l list ${commonArgs()}
    compgen -G "${params.fastp_dir}/*_R1.fastq.gz" > /dev/null || { echo "No R1 cleaned reads were produced in ${params.fastp_dir}" >&2; exit 1; }
    compgen -G "${params.fastp_dir}/*_R2.fastq.gz" > /dev/null || { echo "No R2 cleaned reads were produced in ${params.fastp_dir}" >&2; exit 1; }
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

    output:
    path "${params.spades_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_01_fastp "${params.fastp_dir}"
    ${rgbeppCmd()} -f assembly -l list ${commonArgs()}
    compgen -G "${params.spades_dir}/scaffolds/*.fasta" > /dev/null || { echo "No scaffold FASTA files were produced in ${params.spades_dir}/scaffolds" >&2; exit 1; }
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

    output:
    path "${params.spades_dir}", emit: assembly
    path "${params.bowtie2_dir}", emit: map

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_01_fastp "${params.fastp_dir}"
    cp -aL input_02_spades "${params.spades_dir}"
    ${rgbeppCmd()} -f map -l list -r ${reference} ${commonArgs()}
    compgen -G "${params.spades_dir}/fasta/*.fasta" > /dev/null || { echo "No mapped reference FASTA files were produced in ${params.spades_dir}/fasta" >&2; exit 1; }
    compgen -G "${params.bowtie2_dir}/*.bam" > /dev/null || { echo "No mapping BAM files were produced in ${params.bowtie2_dir}" >&2; exit 1; }
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

    output:
    path "${params.bam_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_03_bowtie2 "${params.bowtie2_dir}"
    ${rgbeppCmd()} -f postmap -l list ${commonArgs()}
    compgen -G "${params.bam_dir}/*.bam" > /dev/null || { echo "No post-mapping BAM files were produced in ${params.bam_dir}" >&2; exit 1; }
    compgen -G "${params.bam_dir}/*.bam.bai" > /dev/null || { echo "No BAM index files were produced in ${params.bam_dir}" >&2; exit 1; }
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

    output:
    path "${params.vcf_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_02_spades "${params.spades_dir}"
    cp -aL input_04_bam "${params.bam_dir}"
    ${rgbeppCmd()} -f varcall -l list ${commonArgs()}
    compgen -G "${params.vcf_dir}/*.vcf.gz" > /dev/null || { echo "No VCF files were produced in ${params.vcf_dir}" >&2; exit 1; }
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

    output:
    path "${params.consen_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_02_spades "${params.spades_dir}"
    cp -aL input_05_vcf "${params.vcf_dir}"
    ${rgbeppCmd()} -f consen -l list -g genes ${commonArgs()}
    compgen -G "${params.consen_dir}/taxa/*.fasta" > /dev/null || { echo "No consensus taxa FASTA files were produced in ${params.consen_dir}/taxa" >&2; exit 1; }
    compgen -G "${params.consen_dir}/gene/*.fasta" > /dev/null || { echo "No consensus gene FASTA files were produced in ${params.consen_dir}/gene" >&2; exit 1; }
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

    output:
    path "${params.consen_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_06_consen "${params.consen_dir}"
    ${rgbeppCmd()} -f codon -g genes -r ${reference} ${commonArgs()}
    compgen -G "${params.consen_dir}/gene/*.fasta" > /dev/null || { echo "No codon-filtered gene FASTA files were produced in ${params.consen_dir}/gene" >&2; exit 1; }
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

    output:
    path "${params.macse_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_06_consen "${params.consen_dir}"
    ${rgbeppCmd()} -f align -g genes ${commonArgs()}
    compgen -G "${params.macse_dir}/AA/*.fasta" > /dev/null || { echo "No MACSE AA alignments were produced in ${params.macse_dir}/AA" >&2; exit 1; }
    compgen -G "${params.macse_dir}/NT/*.fasta" > /dev/null || { echo "No MACSE NT alignments were produced in ${params.macse_dir}/NT" >&2; exit 1; }
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

    output:
    path "${params.macse_dir}", emit: align
    path 'paralog.csv', emit: paralog

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_07_macse "${params.macse_dir}"
    ${rgbeppCmd()} -f ortholog -l list -g genes -r ${reference} ${commonArgs()}
    [[ -e paralog.csv ]] || { echo "Paralog report was not produced: paralog.csv" >&2; exit 1; }
    compgen -G "${params.macse_dir}/taxa/*.tsv" > /dev/null || { echo "No ortholog TSV files were produced in ${params.macse_dir}/taxa" >&2; exit 1; }
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

    output:
    path "${params.trimal_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_07_macse "${params.macse_dir}"
    ${rgbeppCmd()} -f trim -g genes ${commonArgs()}
    compgen -G "${params.trimal_dir}/NT/*.fasta" > /dev/null || { echo "No trimmed FASTA files were produced in ${params.trimal_dir}/NT" >&2; exit 1; }
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

    output:
    path "${params.trimal_dir}"

    script:
    """
    set -euo pipefail
    ${configSetup()}
    cp -aL input_08_trimal "${params.trimal_dir}"
    ${rgbeppCmd()} -f concat -g genes ${commonArgs()}
    [[ -s "${params.trimal_dir}/concat.fasta.fasta" ]] || { echo "Concatenated FASTA was not produced: ${params.trimal_dir}/concat.fasta.fasta" >&2; exit 1; }
    """
}

workflow {
    if (params.help) {
        log.info helpMessage()
        System.exit(0)
    }

    if (!params.reference) {
        error 'Missing required parameter: --reference. Run `nextflow run . --help` for usage.'
    }

    validateFunction()
    validateRelativeDirs()

    sample_list_ch = Channel.fromPath(params.list, checkIfExists: true)
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
    raw_reads_ch = Channel.fromPath("${params.raw_dir}/*_R[12].fastq.gz", checkIfExists: true).collect()

    if (params.genes) {
        genes_ch = Channel.fromPath(params.genes, checkIfExists: true)
    } else {
        PREPARE_GENES(reference_ch)
        genes_ch = PREPARE_GENES.out
    }

    CLEAN(sample_list_ch, raw_reads_ch)

    if (functionEnabled('assembly')) {
        ASSEMBLY(sample_list_ch, CLEAN.out)
    }

    if (functionEnabled('map')) {
        MAP(sample_list_ch, reference_ch, CLEAN.out, ASSEMBLY.out)
    }

    if (functionEnabled('postmap')) {
        POSTMAP(sample_list_ch, MAP.out.map)
    }

    if (functionEnabled('varcall')) {
        VARCALL(sample_list_ch, MAP.out.assembly, POSTMAP.out)
    }

    if (functionEnabled('consen')) {
        CONSEN(sample_list_ch, genes_ch, MAP.out.assembly, VARCALL.out)
    }

    if (functionEnabled('codon') && codonEnabled()) {
        CODON(genes_ch, reference_ch, CONSEN.out)
    }

    if (functionEnabled('align')) {
        if (codonEnabled()) {
            ALIGN(genes_ch, CODON.out)
        } else {
            ALIGN(genes_ch, CONSEN.out)
        }
    }

    if (functionEnabled('ortholog')) {
        ORTHOLOG(sample_list_ch, genes_ch, reference_ch, ALIGN.out)
    }

    if (functionEnabled('trim')) {
        TRIM(genes_ch, ORTHOLOG.out.align)
    }

    if (functionEnabled('concat')) {
        CONCAT(genes_ch, TRIM.out)
    }
}
