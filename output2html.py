import os
import re

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
            if line.startswith('#') or len(line) < 50:
                continue
            name, sequence = re.split(r'\s+', line.strip())
            align[name] = {
                'sequence': sequence,
                'sequence_split': list(sequence),
            }
    return align, ss_cons


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


def write_html(species, align, ss_cons, outlist, family, ga_threshold, best_reversed):
    env = Environment(
        loader=FileSystemLoader(os.path.dirname(os.path.realpath(__file__)))
    )
    env.globals.update(zip=zip)
    template = env.get_template('template.html')
    with open('output.html', 'w') as f_out:
        f_out.write(template.render(species=species, outlist=outlist, align=align, ss_cons=ss_cons, family=family))


def main():
    species, ga_threshold, best_reversed = parse_species('MIPF0000219__mir-484_relabelled/species')
    # align, ss_cons = parse_align('MIPF0000219__mir-484_relabelled/align')
    align, ss_cons = parse_align_with_seed('MIPF0000219__mir-484_relabelled/align-with-seed')
    outlist = parse_outlist('MIPF0000219__mir-484_relabelled/outlist')
    # import pdb; pdb.set_trace()
    write_html(species, align, ss_cons, outlist, 'MIPF0000219__mir-484', ga_threshold, best_reversed)

if __name__ == '__main__':
    main()
