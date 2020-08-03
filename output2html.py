# -*- coding: utf-8 -*-

import os
import re
import sys

import emoji
from jinja2 import Environment, FileSystemLoader


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


def parse_align_with_seed(filename):
    align = dict()
    with open(filename, 'r') as f_in:
        for line in f_in:
            if line.startswith('#=GC SS_cons'):
                _, _, ss_cons = re.split(r'\s+', line.strip())
            if line.startswith('#=GC RF'):
                _, _, rf_line = re.split(r'\s+', line.strip())
            if line.startswith('#') or len(line) < 50:
                continue
            name, sequence = re.split(r'\s+', line.strip())
            align[name] = {
                'sequence': sequence,
                'sequence_split': list(sequence),
            }
    return align, ss_cons, rf_line


def parse_outlist(filename):
    """
    # bits  evalue   seqLabel  name                      overlap  start      end        str  qstart  qend  trunc  species                            extra                  description
    #=====  =======  ========  ========================  =======  =========  =========  ===  ======  ====  =====  =================================  =====================  ==========================================================================================================
      98.3  1.1e-17  FULL      CM000306.1                      -   15357359   15357425    +       1    67     no  Macaca_mulatta_(Rhesus_..[9544]    GA:A;RV:A;SO:N[0.000]  Macaca mulatta chromosome 20, whole genome shotgun sequence.
      98.3  1.1e-17  FULL      CM000330.3                      -   15864341   15864407    +       1    67     no  Pan_troglodytes_(chimpa..[9598]    GA:A;RV:A;SO:N[0.000]  Pan troglodytes isolate Yerkes chimp pedigree #C0471 (Clint) chromosome 16, whole genome shotgun sequence.
    """
    outlist = []
    with open(filename, 'r') as f_in:
        for line in f_in:
            if line.startswith('#'):
                continue
            tabbed_line = re.sub(r'\s{2,}', '\t', line.strip())
            bits, evalue, seqLabel, name, overlap, start, end, str, qstart, qend, trunc, species, extra, description = re.split('\t', tabbed_line)
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
                'seq_name': name if name.startswith('URS00') else '{}/{}-{}'.format(name, start, end),
                'urs_taxid': re.sub(r'\/.+', '', name) if name.startswith('URS00') else '',
            })
    return outlist


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
        current = float(entry['bits'])
        if previous - current > 10:
            big_drop[i] = True
        previous = current
    return big_drop


def write_html(data_path, species, align, ss_cons, rf_line, outlist, family, ga_threshold, best_reversed, big_drops):
    env = Environment(
        loader=FileSystemLoader(os.path.dirname(os.path.realpath(__file__)))
    )
    env.globals.update(zip=zip)
    env.globals['get_emoji'] = get_emoji
    template = env.get_template('template.html')
    with open(os.path.join(data_path, 'output.html'), 'w') as f_out:
        output = template.render(species=species, outlist=outlist, align=align, ss_cons=ss_cons, rf_line=rf_line, family=family, big_drops=big_drops)
        f_out.write(output.encode('utf-8'))


def main(data_path):
    basename = os.path.basename(data_path)
    species, ga_threshold, best_reversed = parse_species(os.path.join(data_path, 'species'))
    # align, ss_cons = parse_align('MIPF0000219__mir-484_relabelled/align')

    cmd = 'esl-reformat pfam {} > {}'.format(os.path.join(data_path, 'align-with-seed'), os.path.join(data_path, 'align-with-seed-pfam'))
    os.system(cmd)

    align, ss_cons, rf_line = parse_align_with_seed(os.path.join(data_path, 'align-with-seed-pfam'))
    outlist = parse_outlist(os.path.join(data_path, 'outlist'))
    big_drops = detect_bit_score_drops(outlist)
    write_html(data_path, species, align, ss_cons, rf_line, outlist, basename, ga_threshold, best_reversed, big_drops)


if __name__ == '__main__':
    main(sys.argv[1])
