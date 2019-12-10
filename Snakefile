configfile: 'config.yaml'

include: 'rules/racon.smk'
include: 'rules/pilon.smk'

wildcard_constraints:
    racon_round = '\d{1,2}',  # racon_round should be a number between 1 and 99
    pilon_round = '\d{1,2}',  # pilon_round should be a number between 1 and 99
    extension = '(fa|fna|fasta)',  # extension is either .fa, .fna, or .fasta

# Get base filename from assembly file by stripping extension and folders
FINAL_BASE, EXTENSION = os.path.splitext(os.path.split(config['assembly'])[-1])
if config['racon_rounds'] > 0:
    # Add racon suffix to pilon base if racon is run
    FINAL_BASE += '_racon%d' % config['racon_rounds']
if config['pilon_rounds'] > 0:
    # Add pilon suffix to final base if pilon is run
    FINAL_BASE += '_pilon%d' % config['pilon_rounds']

# Final assembly file
FINAL_ASSEMBLY = 'output/%s%s' % (FINAL_BASE, EXTENSION)

rule all:
    input:
        FINAL_ASSEMBLY
