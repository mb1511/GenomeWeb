'''
Updated 19DEC2016

@author: Matt Brewer
@organization: University of Bristol
@contact: mb1511@bristol.ac.uk
@summary: Module to handle reading of FASTA files
'''

import re
from string import maketrans
from collections import deque

TRANS_TABLE = maketrans('ATCG', 'TAGC')
RNA_TRANS_TABLE = maketrans('AUCG', 'UAGC')


class Codons(object):
    '''codon list'''
    # TODO: implement different DNA translation tables
    def __init__(self):
        self.codons =  {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
                        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
    def __getitem__(self, key):
        '''catch key errors'''
        try:
            return self.codons[key]
        except KeyError:
            # unidentified codon, assume ANY type
            return 'X'
CODONS = Codons()


class Sequence(object):
    '''sequence object'''

    def __init__(self, name, seq, index=None, seq_type='dna', get_trans=False):
        self.name = name
        self.seq = seq
        self.seq_type = seq_type
        self.index = index
        self.fasta = '\n'.join([self.name, self.seq])
        self._complement = ''
        self._reverse_complement = ''
        self._reverse = ''
        self._ruler = ''
        
        _features = deque(re.findall('(?<=\[)(.*?)(?=\])', name))
        try:
            self.strain = _features.popleft() if '=' not in _features[0] else None
        except IndexError:
            self.strain = None
        if len(_features) > 0:
            self.features = dict(x.split('=') for x in _features if len(x.split('='))==2)
        else:
            self.features = {}
            
        self.t_frames = {1: '',
                         2: '',
                         3: '',
                         4: '',
                         5: '',
                         6: ''}
                
        if self.seq_type == 'dna' and get_trans:
            self.get_all_trans()
            

    def __getitem__(self, key):
        return self.seq[key]

    def __repr__(self):
        return '%s(%s...)' % (self.__class__.__name__, self.name[:20])

    def __len__(self):
        return len(self.seq)
    
    def _trans(self, seq, reverse=False):
        if not reverse:
            return (CODONS[x.group(0)] for x in re.finditer('...', seq))
        else:
            return (CODONS[x.group(0)[::-1]] for x in re.finditer('...', seq))
    
    def get_all_trans(self):
        for f in self.t_frames:
            self.t_frames[f] = self.translate(f, align_to_dna=False)       
    
    @property
    def reverse_complement(self):
        if self._reverse_complement:
            return self._reverse_complement
        if self.seq_type != 'prot':
            self._reverse_complement = self.complement[::-1]
            return self._reverse_complement
        else:
            return 'N/A'         
    
    @property
    def complement(self):
        if self._complement:
            return self._complement
        if self.seq_type == 'dna':
            self._complement = self.seq.translate(TRANS_TABLE)
        elif self.seq_type == 'rna':
            self._complement = self.seq.translate(RNA_TRANS_TABLE)
        else:
            return 'N/A'
        return self._complement
        
    @property
    def reverse(self):
        '''use reversed(seq.seq) for iterator'''
        if self._reverse:
            return self._reverse
        else:
            self._reverse = self.seq[::-1]
            return self._reverse
    
    @property
    def ruler(self):
        if self._ruler:
            return self._ruler
        else:
            mark = ("''''|''''|" for _ in xrange(1, len(self), 10))
            self._ruler = ''.join(mark)
            return self._ruler
        
    def translate(self, frame=4, align_to_dna=True):
        assert frame in range(1,7)
        if align_to_dna:
            space = '  '
            start = ' ' * (frame % 3 if frame % 3 !=0 else 3)
        else:
            space = ''
            start = ''
        if frame < 4:
            return start + space.join(self._trans(self.seq[frame-1:]))
        else:
            return start + space.join(self._trans(self.complement[frame-4:], reverse=True))


def __fasta_read_from_file(file_path):
    '''generator function that reads a file and yields sequence objects'''
    try:
        i = 0
        with open(file_path, 'r') as s:
            try:
                line = next(s)  # is file empty?
            except StopIteration:
                raise IOError() # raising IOError to catch exception
           
            while True:
                if '>' in line:
                    name = line.strip()
                    seq = ''
                    
                    try:
                        line = next(s)
                    except StopIteration:
                        break
                
                    while '>' not in line:
                        seq += line.translate(None, '\n\r \t')
                        try:
                            line = next(s)
                        except StopIteration:
                            break
                    
                    yield Sequence(name, seq, index=i)
                    i += 1
                    
                else:
                    try:
                        line = next(s)
                    except StopIteration:
                        break
    except IOError:
        raise IOError('%s is empty or does not exist.' % file_path)

def __fasta_read_from_text(text):
    '''generator function that reads text and yields sequence objects'''
    
    s = (x.group(0) for x in re.finditer('(.*\n|.+$)', text))
    try:
        i = 0
        line = next(s)                   
    except StopIteration:
        raise IOError() # raising IOError to catch exception
   
    while True:
        if '>' in line:
            name = line.strip()
            seq = ''
            
            try:
                line = next(s)
            except StopIteration:
                break
        
            while '>' not in line:
                seq += line.translate(None, '\n\r \t')
                try:
                    line = next(s)
                except StopIteration:
                    break
            
            yield Sequence(name, seq, index=i)
            i += 1
            
        else:
            try:
                line = next(s)
            except StopIteration:
                break

def format_lines(seq, width=200, frames=[1,], ruler=True):
    r = seq.ruler
    com = seq.complement
    #rev = seq.reverse
    #rcm = seq.reverse_complement
    frm = [seq.translate(i) for i in frames]
    for i in xrange(0, len(seq), width):
        if ruler:
            yield r[i:i + width] + ('    %d' % (i+width)).rjust(14)
        yield seq.seq[i:i + width]
        yield com[i:i + width]
        yield '\n'
        for f in frm:
            yield f[i:i + width]        
        yield '\n' 
    
def fasta_read(file_path=None, from_text=None, generator=True):
    '''returns list of sequence objects from fasta file'''
    
    if generator:
        if file_path is not None:
            return __fasta_read_from_file(file_path)
        else:
            return __fasta_read_from_text(from_text)
    else:
        # evaluate generator
        if file_path is not None:
            return [s for s in __fasta_read_from_file(file_path)]
        else:
            return [s for s in __fasta_read_from_text(from_text)]




if __name__ == '__main__':

    seqs = fasta_read('fnn.faa')
        
    for i in seqs:
        print i.name
        print i.seq
        print i.strain
        break









