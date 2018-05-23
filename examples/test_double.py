'''
Pairwise Example
'''
from os.path import exists
import genomeweb as gw
import svgutils.transform as sg
import numpy as np

root = 'scratch/'
ext = '_reorder.fna'
names = ['R18528', 'R24394']
files = [root + f + ext for f in names]
for f in files:
    assert exists(f),f

# define panel options
size = 500
ssize = size * 2
kw = dict(
    out_file='double.svg',
    working_directory='scratch',
    reorder=False,
    inner_radius=0,
    outer_radius=80,
    label_offset=1.1,
    width=size,
    height=size,
    size=ssize,
    source_angle=np.pi / 8,
    target_angle=-np.pi / 8,
    dummy_axes=2,
    rotation=-90,
    palette_usage=0.8)

gw.create_web(files, x=-size, label_names=names, **kw)

