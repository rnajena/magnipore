import json
import os
import re
from BCBio import GFF
'''
Reading variants json files and parsing keys to be standardized
ATTENTION
---------
    indeces are 1-based
'''
KEY_PARSING = {
    # possible keys in defining mut lists get parsed to GFF3 gene names or product names
    # the regional index points to the region that contains the nucleotide start position
    # if frame shifts present but no new coding region is produced, point to first region
    # (long name, first region, short name)
    'orf1ab':'ORF1ab',
    '1ab':'ORF1ab',
    'orf1a':'ORF1ab polyprotein',
    '1a':'ORF1ab polyprotein',
    'orf1b':'ORF1ab polyprotein',
    '1b':'ORF1ab polyprotein',
    'orf3a':'ORF3a',
    '3a':'ORF3a',
    'orf6':'ORF6',
    '6':'ORF6',
    'orf7a':'ORF7a',
    '7a':'ORF7a',
    'orf7b':'ORF7b',
    '7b':'ORF7b',
    'orf8':'ORF8',
    '8':'ORF8',
    'orf10':'ORF10',
    '10':'ORF10',
    'nsp2':'nsp2',
    'nsp3':'nsp3',
    'nsp4':'nsp4',
    'nsp5':'3C-like proteinase',
    'nsp6':'nsp6',
    'nsp7':'nsp7',
    'nsp8':'nsp8',
    'nsp9':'nsp9',
    'nsp10':'nsp10',
    'nsp11':'nsp11',
    'nsp12':'RNA-dependent RNA polymerase',
    'nsp13':'helicase',
    'nsp14':'3\'-to-5\' exonuclease',
    'nsp15':'endoRNAse',
    's':'surface glycoprotein',
    'spike':'surface glycoprotein', # S
    'e':'envelope protein', # E
    'm':'membrane glycoprotein', # M
    'n':'nucleocapsid phosphoprotein', # N
    'nuc':'SNP',
}

LONG2SHORT = {
    'ORF1ab polyprotein' : 'ORF1ab',
    '3C-like proteinase' : 'proteinase',
    'RNA-dependent RNA polymerase' : 'RdRp',
    '3\'-to-5\' exonuclease' : 'exonuclease',
    'surface glycoprotein' : 'S',
    'envelope protein' : 'E',
    'membrane glycoprotein' : 'M',
    'nucleocapsid phosphoprotein' : 'N',
    '2\'-O-ribose methyltransferase' : 'methyltransferase'
}

class VariantsJSONParser():
    '''
    Reading variants json files and parsing keys to be standardized
    ATTENTION
    ---------
        indices are 1-based
    '''

    def __init__(self, json_path : str, gff3 : str) -> None:
        '''
        Parameters
        ----------
        json_path : str
            path to variant json file
        gff3 : str
            path to annotation file
        '''
        self.json : json = None
        self.muts = []
        self.file = json_path
        self.gff3 = gff3
        self.variant : str = ''

    def read(self) -> None:
        '''
        reads variant json files from given directory
        leave mutations empty if json file does not exist
        '''

        self.anno, self.redundant_anno = readGFF3(self.gff3)

        if self.file is not None:

            if os.path.exists(self.file):
                data = json.load(open(self.file, 'r'))
            else:
                return None

            try:
                self.variant : str = data['variant']['Pango_lineages'][0]
            except:
                self.variant : str = data['variant']['lineage_name']

            self._readMuts(data['sites'])

    def _readMuts(self, muts : dict) -> dict:
        for entry in muts:
            # special case deletion with del:position:length
            if entry.startswith('del'):
                # del:5:5 should look like this
                # Wuhan     AAACCACTTGTT
                # BA.1      AAAC-----GTT
                # 1-based   1234-----567
                # I want to save 4 and 5, the deletion bordering indices -> 5-1 and 5
                gene, mutation, _ = entry.split(':')
                mutation = int(mutation)
                self.muts.extend([int(mutation)-1, int(mutation)])
            
            # special case insertion WHICH IS CALLED NUC (same as SNP), UUFFFF
            elif 'nuc' in entry and '+' in entry:
                # del:5:5 should look like this
                # BA.1      ACCACTTGT
                # Wuhan     AC-----GT
                # 1-based   123456789
                # I want to save 2 and 8, the insertion bordering indices -> 3-1 and 3+len(insertion)
                gene, mutation = entry.split(':')
                position, insertion = mutation.split('+')
                position = int(position)
                self.muts.extend([position-1, position+len(insertion)])

            # substitutions and deletion with substitution to '-' (json list is very consistent ... NOT)
            else:
                gene, mutation = entry.split(':')
                m = re.search('^(\D+)(\d+)', mutation)
                l = len(m.group(1))
                position = int(m.group(2))

                gene = KEY_PARSING.get(gene.lower())
                    
                if gene == 'SNP':
                    self.muts.append(position)
                else:
                    self.muts.extend(
                        amino2NuclPos(position, l, self.redundant_anno[gene][1])
                        )

    def getMuts(self) -> tuple: # [str, list]
        '''
        Returns
        -------
        self.variant : str
            pangolin lineage variant
        self.muts : list
            nucleotide positions of defining mutations according to given annotation and defining mutation list

        ATTENTION
        ---------
            positions are 1-based
        '''
        return self.variant, self.muts

    def getAnno(self) -> dict:
        '''
        Returns
        -------
        self.anno : dict
            0-based [start, end)
            {(start, end) : (depth, geneName)}
        '''
        return self.anno

    def getRedundantAnno(self) -> dict:
        '''
        Returns
        -------
        self.anno : dict
            0-based [start, end)
            {geneName : (depth, start, stop)}
        '''
        return self.redundant_anno

    def getNumMuts(self) -> int:
        return len(self.muts)

def amino2NuclPos(aminoPosition : int, nAminos : int, geneStart : int) -> list:
    '''
    Convert gene amino acid position to nucleotide genome positions 1-based
    '''
    return [geneStart + (aminoPosition-1)*3 + i for i in range(3* nAminos)]

def readGFF3(file : str) -> dict:
    return rec_features_extraction(next(GFF.parse(open(file))).features)

def rec_features_extraction(records, depth : int = 0, features : dict = {}, redundant_features : dict = {}):
    '''
    Returns
    -------
    features : dict
        0-based [start, stop)
        {(start, stop) : (depth, geneName)}
    redundant_features : dict
        0-based [start, stop)
        {geneName : (depth, start, stop)}
    '''

    for rec in records: # list of SeqFeature

        # start is 1-based, transform to 0-based
        # end is included in 1-based, excluded (like python slices) in 0-based
        key = (rec.location.nofuzzy_start - 1, rec.location.nofuzzy_end)

        if key not in features:

            try:
                features[key] = (depth, rec.qualifiers['product'][0])
            except:
                features[key] = (depth, rec.qualifiers['Name'][0])

        try:
            redundant_features[rec.qualifiers['product'][0]] = (depth, rec.location.nofuzzy_start, rec.location.nofuzzy_end)
        except:
            redundant_features[rec.qualifiers['Name'][0]] = (depth, rec.location.nofuzzy_start, rec.location.nofuzzy_end)

        if len(rec.sub_features) > 0:
            rec_features_extraction(rec.sub_features, depth+1, features, redundant_features)

    return features, redundant_features

if __name__ == '__main__':

    gff3 = '/data/fass5/reads/sarscov_kiel/annotation/delta_gisaid.gff3'
    genes = readGFF3(gff3)
    for key in genes:
        print(key, '\t', genes[key])

    # print('==================================')

    # s = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'variants', 'cB.1.1.7.json')
    # parser = VariantsJSONParser('', gff3)
    # parser.read()
    # d = parser.getGeneStarts()
    # for key in d:
    #     print(key, '\t', d[key])
    # variant, muts = parser.getMuts()
    # print(variant, '\n', muts)