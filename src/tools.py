from pylab import *
from IPython import nbformat

def run_notebook(nbfile):
    with open(nbfile) as f:
        nb = nbformat.read(f, 4)
    ip = get_ipython()
    for cell in nb.cells:
        if cell.cell_type != 'code':
            continue
        ip.run_cell(cell.source)