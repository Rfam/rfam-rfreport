#!/usr/bin/env python

# -*- coding: utf-8 -*-

import datetime
import os
import re
import sys

from collections import defaultdict

# activate virtualenv
folder_of_script = os.path.dirname(os.path.realpath(__file__))
activate_script = os.path.join(folder_of_script, 'env', 'bin', 'activate_this.py')
execfile(activate_script, dict(__file__=activate_script))

import click
import emoji
from jinja2 import Environment, FileSystemLoader
import requests


rnacentral_metadata = {}


def parse_species(filename):
    species = []
    ga_threshold = None
    best_reversed = None
    with open(filename, 'r') as f_in:
        for line in f_in:
            m = re.search(r'CURRENT GA THRESHOLD: (.+) BITS', line)
            if m:
                ga_threshold = m.group(1)
                species.append(m.group(0))
                continue
            m = re.search(r'BEST REVERSED HIT E-VALUE: (.+) ', line)
            if m:
                best_reversed = m.group(1)
                species.append(m.group(0))
                continue
            if line.startswith('#'):
                continue
            tabbed_line = re.sub(r'\s{2,}', '\t', line.strip())
            bits, evalue, seqLabel, name, overlap, ncbiId, species_name, extra, taxString = re.split('\t', tabbed_line)
            if seqLabel == 'FULL-SEED':
                continue
            if name.startswith('URS00'):
                seq_name = name
                urs_taxid = re.sub(r'\/.+', '', name) if name.startswith('URS00') else ''
                taxString, _, species_name = get_rnacentral_metadata(urs_taxid)
                ncbiId = urs_taxid.split('_')[1]
            species.append({
                'bits': bits,
                'evalue': evalue,
                'seqLabel': seqLabel,
                'name': name,
                'overlap': overlap,
                'ncbiId': ncbiId,
                'species': species_name,
                'extra': extra,
                'taxString': taxString,
                'tax_string_split': [x.strip().replace('.', '') for x in taxString.split(';')],
            })
    return species, ga_threshold, best_reversed


def parse_align(filename):
    align = []
    with open(filename, 'r') as f_in:
        for line in f_in:
            if line.startswith('#=GC SS_cons'):
                _, _, ss_cons = re.split(r'\s+', line.strip())
            if line.startswith('#') or len(line) < 50:
                continue
            name, sequence = re.split(r'\s+', line.strip())
            align.append({
                'name': name,
                'sequence': sequence,
                'sequence_split': list(sequence),
            })
    return align, ss_cons


def parse_align_with_seed(data_path, threshold):
    ss_cons = ''
    rf_line = ''
    align = os.path.join(data_path, 'align-{}'.format(threshold))
    align_with_seed = os.path.join(data_path, 'align-with-seed-{}'.format(threshold))
    align_with_seed_pfam = os.path.join(data_path, 'align-with-seed-pfam-{}'.format(threshold))

    if not os.path.exists(align):
        cmd = 'cd {} && rfmake.pl -t {} -local -a -forcethr -relax && cp align {}'.format(data_path, threshold, align)
        os.system(cmd)
    if not os.path.exists(align_with_seed) or os.stat(align_with_seed).st_size == 0:
        cmd = 'cd {} && esl-reformat fasta {} | cmalign --mapali SEED CM - > {}'
        os.system(cmd.format(data_path, align, align_with_seed))
    if not os.path.exists(align_with_seed_pfam) or os.stat(align_with_seed_pfam).st_size == 0:
        cmd = 'esl-reformat pfam {} > {}'.format(align_with_seed, align_with_seed_pfam)
        os.system(cmd)

    align = dict()
    with open(align_with_seed_pfam, 'r') as f_in:
        for line in f_in:
            if line.startswith('#=GC SS_cons'):
                _, _, ss_cons = re.split(r'\s+', line.strip())
            if line.startswith('#=GC RF'):
                _, _, rf_line = re.split(r'\s+', line.strip())
            if line.startswith('#') or len(line) < 50:
                continue
            name, sequence = re.split(r'\s+', line.strip())
            if re.search(r'^\d+\|', name):
                name = re.sub(r'^\d+\|', '', name)
            align[name] = {
                'sequence': sequence,
                'sequence_split': list(sequence),
            }
    return align, ss_cons, rf_line


def parse_outlist(filename, maxhits):
    """
    # bits  evalue   seqLabel  name                      overlap  start      end        str  qstart  qend  trunc  species                            extra                  description
    #=====  =======  ========  ========================  =======  =========  =========  ===  ======  ====  =====  =================================  =====================  ==========================================================================================================
      98.3  1.1e-17  FULL      CM000306.1                      -   15357359   15357425    +       1    67     no  Macaca_mulatta_(Rhesus_..[9544]    GA:A;RV:A;SO:N[0.000]  Macaca mulatta chromosome 20, whole genome shotgun sequence.
      98.3  1.1e-17  FULL      CM000330.3                      -   15864341   15864407    +       1    67     no  Pan_troglodytes_(chimpa..[9598]    GA:A;RV:A;SO:N[0.000]  Pan troglodytes isolate Yerkes chimp pedigree #C0471 (Clint) chromosome 16, whole genome shotgun sequence.
    """
    outlist = []
    found_hits_below_reversed = False
    num_hits_below_reversed = 0
    with open(filename, 'r') as f_in:
        for line in f_in:
            m = re.search(r'CURRENT GA THRESHOLD: (.+) BITS', line)
            if m:
                outlist.append(m.group(0))
                continue
            m = re.search(r'BEST REVERSED HIT E-VALUE: (.+) ', line)
            if m:
                outlist.append(m.group(0))
                found_hits_below_reversed = True
                continue
            if line.startswith('#'):
                continue
            tabbed_line = re.sub(r'\s{2,}', '\t', line.strip())
            try:
                bits, evalue, seqLabel, name, overlap, start, end, str, qstart, qend, trunc, species, extra, description = re.split('\t', tabbed_line)
            except:
                bits, evalue, seqLabel, name, overlap, start, end, str, qstart, qend, trunc, species, extra = re.split('\t', tabbed_line)
                description = ''
            if seqLabel == 'FULL-SEED':
                continue
            urs_taxid = re.sub(r'\/.+', '', name) if name.startswith('URS00') else ''
            if urs_taxid:
                seq_name = name
                _, description, species = get_rnacentral_metadata(urs_taxid)
            elif seqLabel == 'SEED':
                seq_name = name
            elif not re.match(r'^(\S+)\/(\d+)\-(\d+)\s*', name):
                seq_name = '{}/{}-{}'.format(name, start, end)
            else:
                seq_name = name
            if found_hits_below_reversed:
                num_hits_below_reversed += 1
            outlist.append({
                'bits': bits,
                'evalue': evalue,
                'seqLabel': seqLabel,
                'name': name,
                'overlap': overlap,
                'start': start,
                'end': end,
                'str': str,
                'qstart': qstart,
                'qend': qend,
                'trunc': trunc,
                'species': species,
                'extra': extra,
                'description': description,
                'seq_name': seq_name,
                'urs_taxid': urs_taxid,
            })
    print('Found {} hits below reversed'.format(num_hits_below_reversed))
    return outlist[:maxhits], num_hits_below_reversed


def parse_mature_mirna_file(filename):
    """
    Transform:
    URS000000076D_6239	15	37	53	74
    into:
    [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74]
    """
    data = dict()
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), filename), 'r') as f_in:
        for line in f_in:
            fields = line.strip().split('\t')
            if (len(fields) - 1) % 2 != 0:
                print('Unusual number of tabs found in line {}'.format(line))
                continue
            data[fields[0]] = []
            flanks = [int(x) for x in fields[1:]]
            while flanks:
                right = flanks.pop()
                left = flanks.pop()
                data[fields[0]] += range(left, right + 1)
    return data


def get_mature_mirna_locations(mature_mirna, outlist, align):
    data = dict()
    for row in outlist:
        if not isinstance(row, dict) or not row['urs_taxid']:
            continue
        if row['urs_taxid'] not in mature_mirna:
            continue
        aligned_sequence = align[row['seq_name']]['sequence']
        mature_ids = mature_mirna[row['urs_taxid']]
        matures = len(aligned_sequence) * [0]
        m = re.match(r'URS\w{10}_\d+\/(\d+)-\d+', row['seq_name'])
        if m:
            urs_taxid_start = int(m.group(1))
        else:
            import pdb; pdb.set_trace()
        seq_id = -1
        for i, nt in enumerate(aligned_sequence):
            if nt.upper() in ['A', 'C', 'G', 'U']:
                seq_id += 1
                if seq_id + urs_taxid_start - 1 in mature_ids:
                    matures[i] = 1
        data[row['urs_taxid']] = matures
    return data


def get_rnacentral_metadata(urs_taxid):
    url = 'http://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral?query={}&fields=description,tax_string,species&format=json'
    if urs_taxid in rnacentral_metadata:
        return rnacentral_metadata[urs_taxid]
    tax_string = ''
    description = ''
    species = ''
    try:
        data = requests.get(url.format(urs_taxid))
        print('fetching {}'.format(urs_taxid))
        if data.json()['hitCount'] == 1:
            tax_string = data.json()['entries'][0]['fields']['tax_string'][0]
            description = data.json()['entries'][0]['fields']['description'][0]
            species = data.json()['entries'][0]['fields']['species'][0]
            rnacentral_metadata[urs_taxid] = (tax_string, description, species)
    except:
        print('Error fetching metadata for {}'.format(urs_taxid))
    return (tax_string, description, species)


def get_emoji(tax_string):
    mapping = {
        'Primates': ':monkey_face:',
        'Viridiplantae': ':herb:',
        'Mollusca': ':oyster:',
        'Suidae': ':pig:',
        'Camelidae': ':camel:',
        'Bovinae': ':cow_face:',
        'Equidae': ':horse_face:',
        'Canidae': ':dog_face:',
        'Rodentia': ':mouse:',
        'Erinaceidae': ':hedgehog:',
        'Chiroptera': ':bat:',
        'Felinae': ':cat_face:',
        'Bacteria': ':microbe:',
        'Ailuropoda': ':panda_face:',
        'Proboscidea': ':elephant:',
        'Ovis': ':ewe:',
        'Xenopodinae': ':frog_face:',
        'Aves': ':bird:',
        'Insecta': ':cricket:',
    }
    found = False
    for taxon, emoji_string in mapping.iteritems():
        if taxon in tax_string:
            return emoji.emojize(emoji_string)
    if not found:
        return ''


def detect_bit_score_drops(outlist):
    big_drop = [False] * len(outlist)
    previous = float(outlist[0]['bits'])
    for i, entry in enumerate(outlist):
        if 'bits' in entry:
            current = float(entry['bits'])
        if previous - current > 10:
            big_drop[i] = True
        previous = current
    return big_drop


def get_seed_nts(align, outlist):
    seed_nts = defaultdict(set)
    for row in outlist:
        if 'seqLabel' in row and row['seqLabel'] == 'SEED':
            try:
                sequence = align[row['seq_name']]['sequence_split']
            except:
                print('Error: {} not found in align data'.format(row['seq_name']))
                continue
            for i, nt in enumerate(sequence):
                seed_nts[i].add(nt)
    return seed_nts


def process_tax_string(species):
    seed_taxa = set()
    for row in species:
        if 'seqLabel' in row and row['seqLabel'] == 'SEED':
            fields = row['taxString'].split(';')
            [seed_taxa.add(x.strip().replace('.', '')) for x in fields]
    return seed_taxa


def process_large_outlist(outlist, num_hits_below_reversed):
    outlist_skip = [0] * len(outlist)
    if num_hits_below_reversed < 50:
        return outlist_skip
    species = set()
    below_reversed = False
    for i, row in enumerate(outlist):
        if isinstance(row, basestring):
            if row.startswith('BEST REVERSED'):
                below_reversed = True
            continue
        if not below_reversed:
            continue
        elif row['seqLabel'] == 'SEED':
            continue
        else:
            if row['species'] not in species:
                species.add(row['species'])
            else:
                outlist_skip[i] = 1
    print('Hiding {} hits below reversed'.format(outlist_skip.count(1)))
    return outlist_skip


def write_html(output_path, species, align, ss_cons, rf_line, outlist, family, ga_threshold, best_reversed, big_drops, outlist_skip, seed_nts, mature_mirnas, seed_taxa):
    env = Environment(
        loader=FileSystemLoader(os.path.dirname(os.path.realpath(__file__)))
    )
    env.globals.update(zip=zip)
    env.globals['get_emoji'] = get_emoji
    template = env.get_template('template.html')
    output_file = os.path.join(output_path, '{}.html'.format(family))
    with open(output_file, 'w') as f_out:
        output = template.render(species=species, outlist=outlist, align=align, ss_cons=ss_cons, ss_cons_split=list(ss_cons), rf_line=rf_line, family=family, big_drops=big_drops, outlist_skip=outlist_skip, seed_nts=seed_nts, mature_mirnas=mature_mirnas, seed_taxa=seed_taxa)
        f_out.write(output.encode('utf-8'))
    print('Created file {}'.format(output_file))
    return output_file


def minify_html(html_file):
    if os.stat(html_file).st_size < 1048576: # 1Mb
        return
    unminified_html_file = html_file.replace('.html', '-unminified.html')
    minified_html_file = html_file
    os.rename(html_file, unminified_html_file)
    cmd = 'htmlmin {} {}'.format(unminified_html_file, minified_html_file)
    try:
        print('Minifying html {}'.format(minified_html_file))
        os.system(cmd)
    except:
        print('Minification failed')


def normalise_align_names(align, outlist):
    """
    The align file always contains names in accession/start-stop format
    because it's generated by esl-reformat.
    However, the outlist file contains sequence names as provided in the SEED,
    which can be in an unexpected format.
    This function makes sure the align dictionary can be accessed with both
    versions of sequence names.
    """
    outlist_names = []
    for outlist_row in outlist:
        if 'seq_name' in outlist_row:
            outlist_names.append(outlist_row['seq_name'])
    for align_name in align.keys():
        if align_name in outlist_names:
            continue
        for outlist_name in outlist_names:
            if outlist_name in align_name:
                align[outlist_name] = align[align_name]
    return align


def verify_species_file_exists(input_path):
    species = os.path.join(input_path, 'species')
    desc = os.path.join(input_path, 'DESC')
    if not os.path.exists(species):
        if os.path.exists(desc):
            suffix = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
            new_desc = "_".join([desc, suffix])
            os.rename(desc, new_desc)
        cmd = 'cd {} && rfsearch.pl -nodesc -cnompi -t 30 -ignoresm && cd -'.format(input_path)
        os.system(cmd)


@click.command()
@click.argument('input_path', type=click.Path(exists=True))
@click.option('--output_path', type=click.Path(exists=True), default='output', help='Path to output folder')
@click.option('--maxhits', default=100000, required=False, type=int, show_default=True, help='Maximum number of hits to output')
@click.option('-t', '--threshold', default=30, required=False, type=int, show_default=True, help='Gathering threshold')
def main(input_path, output_path, maxhits, threshold):
    print('Processing files in {}'.format(input_path))
    basename = os.path.basename(os.path.normpath(input_path))
    verify_species_file_exists(input_path)
    species, ga_threshold, best_reversed = parse_species(os.path.join(input_path, 'species'))

    align, ss_cons, rf_line = parse_align_with_seed(input_path, threshold)
    outlist, num_hits_below_reversed = parse_outlist(os.path.join(input_path, 'outlist'), maxhits)
    align = normalise_align_names(align, outlist)
    outlist_skip = process_large_outlist(outlist, num_hits_below_reversed)
    big_drops = detect_bit_score_drops(outlist)
    seed_nts = get_seed_nts(align, outlist)
    mature_mirna_reference = parse_mature_mirna_file('mature-mirna.tsv')
    mature_mirnas = get_mature_mirna_locations(mature_mirna_reference, outlist, align)
    seed_taxa = process_tax_string(species)
    html_file = write_html(output_path, species, align, ss_cons, rf_line, outlist, basename, ga_threshold, best_reversed, big_drops, outlist_skip, seed_nts, mature_mirnas, seed_taxa)
    minify_html(html_file)


if __name__ == '__main__':
    main()
