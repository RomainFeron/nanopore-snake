

def racon_input(wildcards):
    '''
    Input file for Racon.
    '''
    if wildcards.racon_round == '1':  # If first round of Racon, input is the original assembly file
        return config['assembly']
    else:  # If not first round of Racon, input is the output of the previous round
        return 'output/%s_racon%d.%s' % (wildcards.racon_base, int(wildcards.racon_round) - 1, wildcards.extension)


rule minimap2:
    input:
        assembly = racon_input,
        reads = config['nanopore']
    output:
        temp('output/{racon_base}_{extension}_round{racon_round}.paf')
    benchmark:
        'benchmarks/{racon_base}_{extension}_minimap_round{racon_round}.tsv'
    log:
        'logs/{racon_base}_{extension}_minimap_round{racon_round}.txt'
    conda:
        '../envs/polishing.yaml'
    threads:
        config['minimap2']['threads']
    resources:
        memory = config['minimap2']['memory']
    params:
        runtime = config['minimap2']['runtime']
    shell:
        'minimap2 -x map-ont -t {threads} -o {output} {input.assembly} {input.reads} 2> {log}'


rule racon:
    input:
        assembly = racon_input,
        overlaps = 'output/{racon_base}_{extension}_round{racon_round}.paf',
        reads = config['nanopore']
    output:
        'output/{racon_base}_racon{racon_round}.{extension}'
    benchmark:
        'benchmarks/{racon_base}_{extension}_racon_round{racon_round}.tsv'
    log:
        'logs/{racon_base}_{extension}_racon_round{racon_round}.txt'
    conda:
        '../envs/polishing.yaml'
    threads:
        config['racon']['threads']
    resources:
        memory = config['racon']['memory']
    params:
        runtime = config['racon']['runtime']
    shell:
        'racon -t {threads} {input.reads} {input.overlaps} {input.assembly} 1> {output} 2> {log}'
