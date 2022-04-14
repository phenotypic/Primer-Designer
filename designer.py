import ensembl_rest
import argparse
import subprocess
import string
from prettytable import PrettyTable

parser = argparse.ArgumentParser()
parser.add_argument('-s', help='species')
parser.add_argument('-g', help='gene')
parser.add_argument('-c', help='config path')
parser.add_argument('-v', help='verbose output', action='store_true')  # verbose
args = parser.parse_args()

if args.s is None:
    input_species = input('\nInput species name: ')
else:
    input_species = args.s

if args.g is None:
    input_gene = input('\nInput gene name: ')
else:
    input_gene = args.g

if args.c is None:
    config_location = '/opt/homebrew/Cellar/primer3/2.4.0/share/primer3/primer3_config/'
else:
    config_location = args.c

gene = ensembl_rest.symbol_lookup(input_species, input_gene, params={'expand': True})

splice_count = len(gene['Transcript'])

print('\nGene name:', gene['display_name'])
print('Species:', gene['species'])
print('Gene ID:', gene['id'])
print('Description:', gene['description'])
print('Genome:', gene['assembly_name'])
print('Strand direction:', gene['strand'])
print('Splice varaints:', splice_count)

exon_variants = []
for sv in range(splice_count):
    summary = gene['Transcript'][sv]
    exon_count = len(summary['Exon'])
    print('\nSplice variant', sv + 1, '(' + summary['display_name'] + ') has', exon_count, 'exons')

    summary['location'] = []
    for x in range(exon_count):
        variant = summary['Exon'][x]
        summary['location'].append((variant['start'] - gene['start'], variant['end'] - gene['start']))

    exon_variants.append(summary)


def exon_values(exon):
    values = set()
    for i, item in enumerate(exon['location']):
        for i in range(item[0], item[1] + 1):
            values.add(i)
    return values


whole_sequence = ensembl_rest.sequence_id(gene['id'])['seq']
if gene['strand'] == -1:
    whole_sequence = whole_sequence[::-1]

if splice_count > 1:
    exons = iter(exon_variants)
    values = exon_values(next(exons))
    for exon in exons:
        values.intersection_update(exon_values(exon))

    sorted_values = sorted(list(values))
    ranges = []
    for value in sorted_values:
        if ranges and value == ranges[-1][1] + 1:
            ranges[-1] = (ranges[-1][0], value)
        else:
            ranges.append((value, value))
else:
    ranges = exon_variants[0]['location']

search_sequence = []
for i in range(len(ranges)):
    start = ranges[i][0]
    end = ranges[i][1]
    if gene['strand'] == 1:
        search_sequence.append(whole_sequence[start:end + 1])
    else:
        search_sequence.append(whole_sequence[start:end + 1][::-1])
        if splice_count > 1:
            search_sequence.reverse()

if splice_count > 1:
    print('\nFound', len(search_sequence), 'overlapping sequences to search')

for i, item in enumerate(search_sequence):
    print('\nSequence', i + 1)
    print(item)

primer_params = {
    'Primer min size': 20,
    'Primer opt size': 23,
    'Primer max size': 26,
    'Primer min Tm': 56,
    'Primer opt Tm': 60,
    'Primer max Tm': 64,
    'Product size range': '300-500',
}


def renderList():
    list = PrettyTable(['Index', 'Parameter', 'Value'])

    for index, (key, value) in enumerate(primer_params.items()):
        list.add_row([index + 1, key, value])
    print('\nPrimer search parameters:')
    print(list)


def selectionMenu():
    response = input('\nEnter index of parameter to edit or press return: ').lower()
    if response in string.ascii_lowercase:
        return
    else:
        field = list(primer_params)[int(response) - 1]
        value = input(field + ': ')
        primer_params[field] = value
        renderList()
        selectionMenu()


renderList()
selectionMenu()
in_args = list(primer_params.values())

in_file = '''SEQUENCE_ID=ensembl
SEQUENCE_TEMPLATE={0}
PRIMER_FIRST_BASE_INDEX=1
PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1
PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_LIBERAL_BASE=1
PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0
PRIMER_LOWERCASE_MASKING=0
PRIMER_PICK_ANYWAY=0
PRIMER_EXPLAIN_FLAG=1
PRIMER_MASK_TEMPLATE=0
PRIMER_TASK=generic
PRIMER_MASK_FAILURE_RATE=0.1
PRIMER_MASK_5P_DIRECTION=1
PRIMER_MASK_3P_DIRECTION=0
PRIMER_MIN_QUALITY=0
PRIMER_MIN_END_QUALITY=0
PRIMER_QUALITY_RANGE_MIN=0
PRIMER_QUALITY_RANGE_MAX=100
PRIMER_MIN_SIZE={1}
PRIMER_OPT_SIZE={2}
PRIMER_MAX_SIZE={3}
PRIMER_MIN_TM={4}
PRIMER_OPT_TM={5}
PRIMER_MAX_TM={6}
PRIMER_PAIR_MAX_DIFF_TM=5.0
PRIMER_TM_FORMULA=1
PRIMER_PRODUCT_MIN_TM=-1000000.0
PRIMER_PRODUCT_OPT_TM=0.0
PRIMER_PRODUCT_MAX_TM=1000000.0
PRIMER_MIN_GC=30.0
PRIMER_OPT_GC_PERCENT=50.0
PRIMER_MAX_GC=70.0
PRIMER_PRODUCT_SIZE_RANGE={7}
PRIMER_NUM_RETURN=5
PRIMER_MAX_END_STABILITY=9.0
PRIMER_MAX_LIBRARY_MISPRIMING=12.00
PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00
PRIMER_MAX_SELF_ANY_TH=45.0
PRIMER_MAX_SELF_END_TH=35.0
PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0
PRIMER_PAIR_MAX_COMPL_END_TH=35.0
PRIMER_MAX_HAIRPIN_TH=24.0
PRIMER_MAX_SELF_ANY=8.00
PRIMER_MAX_SELF_END=3.00
PRIMER_PAIR_MAX_COMPL_ANY=8.00
PRIMER_PAIR_MAX_COMPL_END=3.00
PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00
PRIMER_MAX_TEMPLATE_MISPRIMING=12.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_MAX_POLY_X=4
PRIMER_INSIDE_PENALTY=-1.0
PRIMER_OUTSIDE_PENALTY=0
PRIMER_GC_CLAMP=0
PRIMER_MAX_END_GC=5
PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=3
PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=3
PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION=7
PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=4
PRIMER_SALT_MONOVALENT=50.0
PRIMER_SALT_CORRECTIONS=1
PRIMER_SALT_DIVALENT=1.5
PRIMER_DNTP_CONC=0.6
PRIMER_DNA_CONC=50.0
PRIMER_SEQUENCING_SPACING=500
PRIMER_SEQUENCING_INTERVAL=250
PRIMER_SEQUENCING_LEAD=50
PRIMER_SEQUENCING_ACCURACY=20
PRIMER_WT_SIZE_LT=1.0
PRIMER_WT_SIZE_GT=1.0
PRIMER_WT_TM_LT=1.0
PRIMER_WT_TM_GT=1.0
PRIMER_WT_GC_PERCENT_LT=0.0
PRIMER_WT_GC_PERCENT_GT=0.0
PRIMER_WT_SELF_ANY_TH=0.0
PRIMER_WT_SELF_END_TH=0.0
PRIMER_WT_HAIRPIN_TH=0.0
PRIMER_WT_TEMPLATE_MISPRIMING_TH=0.0
PRIMER_WT_SELF_ANY=0.0
PRIMER_WT_SELF_END=0.0
PRIMER_WT_TEMPLATE_MISPRIMING=0.0
PRIMER_WT_NUM_NS=0.0
PRIMER_WT_LIBRARY_MISPRIMING=0.0
PRIMER_WT_SEQ_QUAL=0.0
PRIMER_WT_END_QUAL=0.0
PRIMER_WT_POS_PENALTY=0.0
PRIMER_WT_END_STABILITY=0.0
PRIMER_WT_MASK_FAILURE_RATE=0.0
PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0
PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0
PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0
PRIMER_PAIR_WT_PRODUCT_TM_GT=0.0
PRIMER_PAIR_WT_COMPL_ANY_TH=0.0
PRIMER_PAIR_WT_COMPL_END_TH=0.0
PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH=0.0
PRIMER_PAIR_WT_COMPL_ANY=0.0
PRIMER_PAIR_WT_COMPL_END=0.0
PRIMER_PAIR_WT_TEMPLATE_MISPRIMING=0.0
PRIMER_PAIR_WT_DIFF_TM=0.0
PRIMER_PAIR_WT_LIBRARY_MISPRIMING=0.0
PRIMER_PAIR_WT_PR_PENALTY=1.0
PRIMER_PAIR_WT_IO_PENALTY=0.0
PRIMER_INTERNAL_MIN_SIZE=18
PRIMER_INTERNAL_OPT_SIZE=20
PRIMER_INTERNAL_MAX_SIZE=27
PRIMER_INTERNAL_MIN_TM=57.0
PRIMER_INTERNAL_OPT_TM=60.0
PRIMER_INTERNAL_MAX_TM=63.0
PRIMER_INTERNAL_MIN_GC=20.0
PRIMER_INTERNAL_OPT_GC_PERCENT=50.0
PRIMER_INTERNAL_MAX_GC=80.0
PRIMER_INTERNAL_MAX_SELF_ANY_TH=47.00
PRIMER_INTERNAL_MAX_SELF_END_TH=47.00
PRIMER_INTERNAL_MAX_HAIRPIN_TH=47.00
PRIMER_INTERNAL_MAX_SELF_ANY=12.00
PRIMER_INTERNAL_MAX_SELF_END=12.00
PRIMER_INTERNAL_MIN_QUALITY=0
PRIMER_INTERNAL_MAX_NS_ACCEPTED=0
PRIMER_INTERNAL_MAX_POLY_X=5
PRIMER_INTERNAL_MAX_LIBRARY_MISHYB=12.00
PRIMER_INTERNAL_SALT_MONOVALENT=50.0
PRIMER_INTERNAL_DNA_CONC=50.0
PRIMER_INTERNAL_SALT_DIVALENT=1.5
PRIMER_INTERNAL_DNTP_CONC=0.0
PRIMER_INTERNAL_WT_SIZE_LT=1.0
PRIMER_INTERNAL_WT_SIZE_GT=1.0
PRIMER_INTERNAL_WT_TM_LT=1.0
PRIMER_INTERNAL_WT_TM_GT=1.0
PRIMER_INTERNAL_WT_GC_PERCENT_LT=0.0
PRIMER_INTERNAL_WT_GC_PERCENT_GT=0.0
PRIMER_INTERNAL_WT_SELF_ANY_TH=0.0
PRIMER_INTERNAL_WT_SELF_END_TH=0.0
PRIMER_INTERNAL_WT_HAIRPIN_TH=0.0
PRIMER_INTERNAL_WT_SELF_ANY=0.0
PRIMER_INTERNAL_WT_SELF_END=0.0
PRIMER_INTERNAL_WT_NUM_NS=0.0
PRIMER_INTERNAL_WT_LIBRARY_MISHYB=0.0
PRIMER_INTERNAL_WT_SEQ_QUAL=0.0
PRIMER_INTERNAL_WT_END_QUAL=0.0
PRIMER_THERMODYNAMIC_PARAMETERS_PATH={8}
='''.format(''.join(search_sequence), *in_args, config_location)

output = subprocess.run(['primer3_core'], stdout=subprocess.PIPE, input=in_file, encoding='ascii')
output = output.stdout.replace('\n', '=').split('=')

pair_count = int(output[output.index('PRIMER_PAIR_NUM_RETURNED') + 1])

print('\nFound', pair_count, 'primer pairs')


def result_parser(side, i, s_string):
    if side == 1:
        return output[output.index('PRIMER_PAIR_' + str(i) + s_string) + 1]
    else:
        return output[output.index('PRIMER_' + side + '_' + str(i) + s_string) + 1]


pairs = []
for i in range(pair_count):
    primers = {}
    primers['product_size'] = result_parser(1, i, '_PRODUCT_SIZE')
    primers['compl_any'] = result_parser(1, i, '_COMPL_ANY_TH')
    primers['compl_end'] = result_parser(1, i, '_COMPL_END_TH')
    for index, side in enumerate(['LEFT', 'RIGHT']):
        primers[index] = {}
        position = result_parser(side, i, '').split(',')
        primers[index]['start'] = position[0]
        primers[index]['length'] = position[1]
        primers[index]['sequence'] = result_parser(side, i, '_SEQUENCE')
        primers[index]['tm'] = result_parser(side, i, '_TM')
        primers[index]['self_binding_any'] = result_parser(side, i, '_SELF_ANY_TH')
        primers[index]['any_hairpin'] = result_parser(side, i, '_HAIRPIN_TH')
        primers[index]['end_stability'] = result_parser(side, i, '_END_STABILITY')
        primers[index]['gc_percent'] = result_parser(side, i, '_GC_PERCENT')
    pairs.append(primers)

for index in reversed(range(len(pairs))):
    print('\n\nPrimer pair', index + 1)
    print('---------------------')
    print('Product size:', pairs[index]['product_size'])
    print('Any compl   :', pairs[index]['compl_any'])
    print('End compl   :', pairs[index]['compl_end'])
    list = PrettyTable(['Direction', 'Length', 'Tm', 'GC%', 'Self-binding', 'Hairpin', 'End Stability', 'Sequence'])
    for i, side in enumerate(['Left', 'Right']):
        dict = pairs[index][i]
        list.add_row([side, dict['length'], dict['tm'], dict['gc_percent'], dict['self_binding_any'], dict['any_hairpin'], dict['end_stability'], dict['sequence']])
    print(list)
