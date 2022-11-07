class ANSI:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'
    UP = '\033[1A'
    CLEAR = '\033[K'
    
def complement(seq):
    ret = ''
    for b in seq:
        ret += COMPLEMENT.get(b, 'N')
    return ret

BASEENCODER = {
    'A':0,
    'C':1,
    'G':2,
    'T':3,
    'N':4,
    'Y':5,
    'R':6,
    'S':7,
    'W':8,
    'M':9, 
    'K':10}

BASEDECODER = {
    0:'A',
    1:'C',
    2:'G',
    3:'T',
    4:'N',
    5:'Y',
    6:'R',
    7:'S',
    8:'W',
    9:'M', 
    10:'K'}

STRANDENCODER = {
    '+':0,
    '-':1
}

STRANDDECODER = {
    0:'+',
    1:'-'
}

DATAENCODER = {
    'base':0,
    'omv':1,
    'motif':2,
    'mean':3,
    'std':4,
    'data_density':5,
    'contained_datapoints':6,
    'contained_segments':7,
    'n_datapoints':8,
    'n_segments':9,
    'n_reads':10,
}

COMPLEMENT = {
    'A':'T',
    'C':'G',
    'G':'C',
    'T':'A',
    'N':'N',
    'Y':'R',
    'R':'Y',
    'S':'S',
    'W':'W',
    'M':'K', 
    'K':'M'}

IUPAC = {
    'A':'A',
    'C':'C',
    'G':'G',
    'T':'T',
    'N':'ACGT',
    'Y':'CT',
    'R':'AG',
    'S':'GC',
    'W':'AT',
    'M':'AC',
    'K':'GT'}

OPPOSITION = {
    'A':'CGT',
    'C':'AGT',
    'G':'ACT',
    'T':'ACG',
    'N':'',
    'Y':'AG',
    'R':'CT',
    'S':'AT',
    'W':'GC',
    'M':'GT', 
    'K':'AC'}
