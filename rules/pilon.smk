

def pilon_input(wildcards):
    '''
    Input file for Pilon.
    '''
    if wildcards.pilon_round == '1':  # First round of Pilon
        if config['racon_rounds'] == 0:  # If no Racon rounds, input file is the original assembly file
            return config['assembly']
        else:  # If Racon rounds, input file is the final output from Racon, taken from the
            return 'output/%s.%s' % (wildcards.pilon_base, wildcards.extension)
    else:
        return 'output/%s_pilon%d.%s' % (wildcards.pilon_base, int(wildcards.pilon_round) - 1, wildcards.extension)


rule bwa:
    input:
        assembly = pilon_input,
        reads_r1 = config['illumina']['R1'],
        reads_r2 = config['illumina']['R2']
    output:
        temp('output/{pilon_base}_{extension}_bwa_round{pilon_round}.bam')
    benchmark:
        'benchmarks/{pilon_base}_{extension}_bwa_round{pilon_round}.tsv'
    log:
        'logs/{pilon_base}_{extension}_bwa_round{pilon_round}.txt'
    conda:
        '../envs/polishing.yaml'
    threads:
        config['bwa']['threads']
    resources:
        memory = config['bwa']['memory']
    params:
        runtime = config['bwa']['runtime']
    shell:
        'bwa index {input.assembly} 2>> {log};'
        'bwa mem -t {threads} {input.assembly} {input.reads_r1} {input.reads_r2}  2>> {log} | '
        'samtools sort - > {output} 2>> {log};'
        'samtools index -@ {threads} {output} 2>> {log};'


rule pilon:
    input:
        assembly = pilon_input,
        bam = 'output/{pilon_base}_{extension}_bwa_round{pilon_round}.bam'
    output:
        'output/{pilon_base}_pilon{pilon_round}.{extension}'
    benchmark:
        'benchmarks/{pilon_base}_{extension}_pilon_round{pilon_round}.tsv'
    log:
        'logs/{pilon_base}_{extension}_pilon_round{pilon_round}.txt'
    conda:
        '../envs/polishing.yaml'
    threads:
        config['pilon']['threads']
    resources:
        memory = config['pilon']['memory']
    params:
        runtime = config['pilon']['runtime'],
        output_dir = 'tmp',
        output_prefix = 'output/{pilon_base}_pilon{pilon_round}',
        max_mem = f'{int(config["pilon"]["memory"]) / 1000}G'
    shell:
        'java --changes --fix all -Xmx{params.max_mem} -jar pilon.jar --threads {threads} '
        '--genome {input.assembly} --bam {input.bam} --outdir {params.output_dir} '
        '--output {params.output_prefix} 2>> {log}'
