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
