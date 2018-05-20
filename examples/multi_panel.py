'''
Multi-Panel Example
'''
import glob
from os.path import basename, exists, join, splitext
import genomeweb as gw
import svgutils.transform as sg

root = 'genomes/'
ext = '.fna'
names = [
    '2B2',
    '2B4',
    '2B16',
    'R28427',
    'R5001',
    'R15792',
    'R18932',
    'R30927',
    'R16531',
    '2B17',
    'R26872',
    'R28211',
    'R29976',
    'R30464',
    'R30604',
    'R31249',
    'R32310',
    'R33533',
    'R18528',
    'R24394',
    'R28385',
    'R28400',
    'R32935'
    ]

files = [root + f + ext for f in names]
for f in files:
    assert exists(f),f
working_directory='scratch'

# re-index contigs outside of plot making
reference_genome = [
    'genomes/ref_a.fna', '', '', '',
    'genomes/ref_b.fna', '', '', '', '',
    'genomes/ref_c.fna', '', '', '', '', '', '', '', '',
    'genomes/ref_d.fna', '', '', '', '']

for i, genome in enumerate(files):
    ref = reference_genome[i]
    p_ref = reference_genome[-1]
    if ref:
        ref_to_use = ref
    # only make new db for unique/non-empty/initial
    # genome in referecne list, else use previous entry
    files[i] = gw.reorderctgs.run(
        genome, ref_to_use,
        order_out=join(
            working_directory,
            '%s_reorder.fna' % splitext(basename(genome))[0]),
        make_db=(not i or ref != p_ref) and ref,
        wd=working_directory)


# define panel options
size = 500
ssize = size / 2
kw = dict(
    out_file='multi_panel.svg',
    working_directory=working_directory,
    reorder=False,
    outer_radius=45,
    label_offset=1.2,
    width=size,
    height=size,
    size=ssize,
    font_size=11,
    bezier_max_n=9)

# multi-panel figure 
gw.create_web(files[:4], label_names=names[:4], **kw)
gw.create_web(files[4:9], label_names=names[4:9], x=ssize, append=True, **kw)
gw.create_web(files[9:18], label_names=names[9:18], x=0, y=ssize, append=True, **kw)
fig = gw.create_web(files[18:], label_names=names[18:], x=ssize, y=ssize, append=True, **kw)

# add annotations
fig.append(sg.TextElement(0, 20, 'Species 1', size=12, font='Arial', weight='bold'))
fig.append(sg.TextElement(ssize, 20, 'Species 2 and 3', size=12, font='Arial', weight='bold'))
fig.append(sg.TextElement(0, ssize + 20, 'Species 4', size=12, font='Arial', weight='bold'))
fig.append(sg.TextElement(ssize, ssize + 20, 'Species 5', size=12, font='Arial', weight='bold'))
fig.save('multi_panel.svg')
