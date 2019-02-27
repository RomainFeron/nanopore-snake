configfile: 'config.yaml'

wildcard_constraints:
    racon_round = '\d{1,2}',  # racon_round should be a number between 1 and 99
    pilon_round = '\d{1,2}',  # pilon_round should be a number between 1 and 99
    extension = '.(fa|fna|fasta)',  # extension is either .fa, .fna, or .fasta

# Get base filename from assembly file by stripping extension and folders
FINAL_BASE, EXTENSION = os.path.splitext(os.path.split(config['assembly'])[-1])
if config['racon_rounds'] > 0:
    FINAL_BASE += '_racon%d' % config['racon_rounds']  # Add racon suffix to pilon base if racon is run
if config['pilon_rounds'] > 0:
    FINAL_BASE += '_pilon%d' % config['pilon_rounds']  # Add pilon suffix to final base if pilon is run

# Final assembly file
FINAL_ASSEMBLY = 'output/%s%s' % (FINAL_BASE, EXTENSION)


def racon_input(wildcards):
    '''
    Input file for Racon.
    '''
    if wildcards.racon_round == '1':  # If first round of Racon, input is the original assembly file
        return config['assembly']
    else:  # If not first round of Racon, input is the output of the previous round
        return 'output/%s_racon%d%s' % (wildcards.racon_base, int(wildcards.racon_round) - 1, wildcards.extension)


def pilon_input(wildcards):
    '''
    Input file for Pilon.
    '''
    if wildcards.pilon_round == '1':  # First round of Pilon
        if config['racon_rounds'] == 0:  # If no Racon rounds, input file is the original assembly file
            return config['assembly']
        else:  # If Racon rounds, input file is the final output from Racon, taken from the
            return 'output/%s%s' % (wildcards.pilon_base, wildcards.extension)
    else:
        return 'output/%s_pilon%d%s' % (wildcards.pilon_base, int(wildcards.pilon_round) - 1, wildcards.extension)


rule all:
    input:
        FINAL_ASSEMBLY


rule racon:
    input:
        racon_input
    output:
        'output/{racon_base}_racon{racon_round}{extension}'
    threads:
        config['threads']
    params:
        racon = config['racon_path'],
        minimap2 = config['minimap2_path'],
        reads = config['nanopore'],
        minimap2_output = 'output/racon_minimap2_round{racon_round}.paf'
    shell:
        '{params.minimap2} -t {threads} -x map-ont {input} {params.reads} -o {params.minimap2_output};'
        '{params.racon} -t {threads} {params.reads} {params.minimap2_output} {input} > {output};'
        'rm -f {params.minimap2_output}'


rule pilon:
    input:
        pilon_input
    output:
        'output/{pilon_base}_pilon{pilon_round}{extension}'
    threads:
        config['threads']
    params:
        max_mem = config['max_mem'],
        pilon = config['racon_path'],
        reads_r1 = config['illumina']['R1'],
        reads_r2 = config['illumina']['R2'],
        load_bwa = config['bwa_module'],
        load_samtools = config['samtools_module'],
        load_java = config['java_module'],
        bwa_output = 'output/pilon_bwa_round{pilon_round}.bam',
        output_dir = config['pilon_output_dir'],
        output_prefix = 'output/{pilon_base}_pilon{pilon_round}'
    shell:
        '{params.load_bwa}; {params.load_samtools}; {params.load_java};'
        'bwa index {input};'
        'bwa mem -t {threads} {input} {params.reads_r1} {params.reads_r2} | '
        'samtools sort -@ {threads} - > {params.bwa_output};'
        'java -Xmx{max_mem} -jar {params.pilon} --threads {threads} --changes --fix all '
        '--genome {input} --bam {params.bwa_output} --outdir {params.output_dir} --output {params.output_prefix}'
