'''
This script solves the inability to use pickle exporting
to export a defaultdict. This will allow the defaultdict created by
goodCollapseDictionary.py to be written to a file if needed.
'''
from json import dump
from json import load

def exportDictionary(sequences, pathname):
    dump(sequences, open(pathname, 'w'))

def importDictionary(pathname):
    sequences = load(open(filename))
    return sequences
