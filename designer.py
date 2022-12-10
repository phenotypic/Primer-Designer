import ensembl_rest
import argparse
import subprocess
import string
from prettytable import PrettyTable

parser = argparse.ArgumentParser()
parser.add_argument('-s', help='species')
parser.add_argument('-g', help='gene')
parser.add_argument('-p', help='parameter file name')
parser.add_argument('-t', help='thermal path')
parser.add_argument('-v', help='verbose output', action='store_true')
args = parser.parse_args()

if args.s is None:
    input_species = input('\nInput species name (e.g. zebrafish): ')
else:
    input_species = args.s

if args.g is None:
    input_gene = input('\nInput gene name (e.g. wnt10a): ')
else:
    input_gene = args.g

if args.p is None:
    parameter_file = 'parameters.txt'
else:
    parameter_file = args.p

gene = ensembl_rest.symbol_lookup(input_species, input_gene, params={'expand': True})

splice_count = len(gene['Transcript'])

print('Species:', gene['species'])
print('Genome:', gene['assembly_name'])
print('\nGene name:', gene['display_name'])
print('Description:', gene['description'])
print('Gene ID:', gene['id'])
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
        if args.v:
            print('Exon', str(x + 1) + ':', variant['start'] - gene['start'], '-', variant['end'] - gene['start'])
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


def find_overlaps(variants):
    if splice_count > 1:
        exons = iter(variants)
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
        ranges = variants[0]['location']

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
        if len(search_sequence) > 0:
            print('\nFound', len(search_sequence), 'overlapping sequences')
        else:
            print('\nError: found 0 overlapping sequences')
            indices = [int(i) - 1 for i in input('Manually select splice variants to find overlaps for (e.g. 1,5,12): ').replace(' ', '').split(',')]
            selected_elements = [exon_variants[index] for index in indices]
            search_sequence = find_overlaps(selected_elements)
    return search_sequence


search_sequence = find_overlaps(exon_variants)

for i, item in enumerate(search_sequence):
    print('\nExon', i + 1)
    print(item)

with open(parameter_file, 'r') as file:
    in_file = file.read()


primer_params = {
    'Primer min size': ['PRIMER_MIN_SIZE='],
    'Primer opt size': ['PRIMER_OPT_SIZE='],
    'Primer max size': ['PRIMER_MAX_SIZE='],
    'Primer min Tm': ['PRIMER_MIN_TM='],
    'Primer opt Tm': ['PRIMER_OPT_TM='],
    'Primer max Tm': ['PRIMER_MAX_TM='],
    'Product size range': ['PRIMER_PRODUCT_SIZE_RANGE='],
}

for key, value in primer_params.items():
    original_value = in_file.split(value[0])[1].split('\n')[0]
    primer_params[key].extend((original_value, original_value))


def selectionMenu():
    list = PrettyTable(['Index', 'Parameter', 'Value'])

    for index, (key, value) in enumerate(primer_params.items()):
        list.add_row([index + 1, key, value[1]])
    print('\nPrimer search parameters:')
    print(list)

    response = input('\nEnter index to edit parameter (or press return): ').lower()
    if response in string.ascii_lowercase:
        return
    else:
        field = list(primer_params)[int(response) - 1]
        value = input(field + ': ')
        primer_params[field][1] = value
        selectionMenu()


selectionMenu()

joined_sequence = ''.join(search_sequence)
if args.v:
    print('\nSequence sent to primer3:')
    print(joined_sequence)

for value in primer_params.values():
    in_file = in_file.replace(value[0] + value[2], value[0] + value[1])

in_file = in_file.replace('SEQUENCE_TEMPLATE=', 'SEQUENCE_TEMPLATE=' + joined_sequence)

if args.t is not None:
    config_text = 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH='
    original_value = in_file.split(config_text)[1].split('\n')[0]
    in_file = in_file.replace(config_text + original_value, config_text + args.t)

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

for index in reversed(range(pair_count)):
    print('\n\nPrimer pair', index + 1)
    print('---------------------')
    print('Gene name:', gene['display_name'])
    print('Species:', gene['species'])
    print('Product size:', pairs[index]['product_size'])
    print('Any compl   :', pairs[index]['compl_any'])
    print('End compl   :', pairs[index]['compl_end'])
    list = PrettyTable(['Primer', 'Length', 'Tm', 'GC%', 'Self-binding', 'Hairpin', 'End Stability', 'Sequence'])
    for i, side in enumerate(['Left', 'Right']):
        dict = pairs[index][i]
        list.add_row([side, dict['length'], dict['tm'], dict['gc_percent'], dict['self_binding_any'], dict['any_hairpin'], dict['end_stability'], dict['sequence']])
    print(list)
