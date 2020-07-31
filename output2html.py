import os
import re

from jinja2 import Environment, FileSystemLoader


def parse_species(filename):
    species = []
    with open(filename, 'r') as f_in:
        for line in f_in:
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
    return species


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


def write_html(species, align, ss_cons, family):
    env = Environment(
        loader=FileSystemLoader(os.path.dirname(os.path.realpath(__file__)))
    )
    env.globals.update(zip=zip)
    template = env.get_template('template.html')
    with open('output.html', 'w') as f_out:
        f_out.write(template.render(species=species, align=align, ss_cons=ss_cons, family=family))


def main():
    species = parse_species('MIPF0000219__mir-484_relabelled/species')
    align, ss_cons = parse_align('MIPF0000219__mir-484_relabelled/align')
    write_html(species, align, ss_cons, 'MIPF0000219__mir-484')

if __name__ == '__main__':
    main()
