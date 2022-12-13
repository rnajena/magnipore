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

MAGNIPORE_COLUMNS = [
    'strand',                   #1
    'td_score',                 #2
    'kl_divergence',            #3
    'bayesian_p',               #4
    'signal_type',              #5
    'ref_1',                    #6
    'pos_1',                    #7
    'base_1',                   #8
    'motif_1',                  #9
    'signal_mean_1',            #10
    'signal_std_1',             #11
    'n_datapoints_1',           #12
    'contained_datapoints_1',   #13
    'n_segments_1',             #14
    'contained_segments_1',     #15
    'n_reads_1',                #16
    'ref_2',                    #17
    'pos_2',                    #18
    'base_2',                   #19
    'motif_2',                  #20
    'signal_mean_2',            #21
    'signal_std_2',             #22
    'n_datapoints_2',           #23
    'contained_datapoints_2',   #24
    'n_segments_2',             #25
    'contained_segments_2',     #26
    'n_reads_2',                #27
    ]