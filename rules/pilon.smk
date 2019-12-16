

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
        'rm *.{wildcards.extension}.*'


checkpoint get_contigs:
    input:
        assembly = pilon_input
    output:
        directory('output/.{pilon_base}_pilon{pilon_round}_{extension}/')
    benchmark:
        'benchmarks/get_contigs_{pilon_base}_pilon{pilon_round}_{extension}.tsv'
    log:
        'logs/get_contigs_{pilon_base}_pilon{pilon_round}_{extension}.txt'
    params:
        output_dir = lambda wildcards, output: output[0].replace('/', '\\/')
    shell:
        'mkdir -p {output};'
        'grep ">" {input} 2> {log} | sed "s/>/{params.output_dir}/g" 2>> {log} | xargs touch 2>> {log}'


rule pilon_contig:
    input:
        assembly = pilon_input,
        bam = 'output/{pilon_base}_{extension}_bwa_round{pilon_round}.bam',
        contig = 'output/.{pilon_base}_pilon{pilon_round}_{extension}/{contig}'
    output:
        temp('output/pilon/{pilon_base}_pilon{pilon_round}_{contig}.{extension}')
    benchmark:
        'benchmarks/{pilon_base}_{extension}_pilon_round{pilon_round}/{contig}.tsv'
    log:
        'logs/{pilon_base}_{extension}_pilon_round{pilon_round}/{contig}.txt'
    conda:
        '../envs/polishing.yaml'
    threads:
        config['pilon']['threads']
    resources:
        memory = lambda wildcards, attempt: config['pilon']['memory'] * attempt
    params:
        runtime = config['pilon']['runtime'],
        output_dir = 'output/pilon',
        output_prefix = '{pilon_base}_pilon{pilon_round}_{contig}',
        max_mem = lambda wildcards, resources: f'{int(resources.memory / 1000) * attempt}G'
    shell:
        'pilon -Xmx{params.max_mem} --genome {input.assembly} --bam {input.bam} --output {params.output_prefix} --outdir {params.output_dir} '
        '--changes --fix all --threads {threads} --targets {wildcards.contig} 2> {log}'


def pilon_all_contigs(wildcards):
    '''
    '''
    pilon_contig_output = checkpoints.get_contigs.get(**wildcards).output[0]
    pilon_base = wildcards.pilon_base
    pilon_round = wildcards.pilon_round
    extension = wildcards.extension
    batch_list = expand('output/pilon/{pilon_base}_pilon{pilon_round}_{contig}.{extension}',
                        pilon_base=pilon_base,
                        pilon_round=pilon_round,
                        extension=extension,
                        contig=glob_wildcards(f'output/.{pilon_base}_pilon{pilon_round}_{extension}/{{contig}}').contig)
    return batch_list


rule pilon_round:
    input:
        pilon_all_contigs
    output:
        'output/{pilon_base}_pilon{pilon_round}.{extension}'
    benchmark:
        'benchmarks/{pilon_base}_{extension}_pilon_round{pilon_round}.tsv'
    log:
        'logs/{pilon_base}_{extension}_pilon_round{pilon_round}.txt'
    shell:
        'cat {input} > {output} 2> {log}'
