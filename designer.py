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

with open('parameters.txt', 'r') as file:
    in_file = file.read().format(''.join(search_sequence), *in_args, config_location)

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
    print('Gene name:', gene['display_name'])
    print('Species:', gene['species'])
    print('Product size:', pairs[index]['product_size'])
    print('Any compl   :', pairs[index]['compl_any'])
    print('End compl   :', pairs[index]['compl_end'])
    list = PrettyTable(['Direction', 'Length', 'Tm', 'GC%', 'Self-binding', 'Hairpin', 'End Stability', 'Sequence'])
    for i, side in enumerate(['Left', 'Right']):
        dict = pairs[index][i]
        list.add_row([side, dict['length'], dict['tm'], dict['gc_percent'], dict['self_binding_any'], dict['any_hairpin'], dict['end_stability'], dict['sequence']])
    print(list)
