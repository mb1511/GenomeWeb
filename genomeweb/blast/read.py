'''
Created on 2 Sep 2016

@author: Matt

Module to quickly read BLAST xml outputs
'''
from __future__ import print_function
from builtins import str, range
import re
import mmap

import time
read_time = 0
i_time = time.time()

class _make_class:
    def __init__(self, dct):
        self.__dict__.update(dct)

class BlastOutput:
    '''
    Reads blast output xml file and formats to a multilayered callable object
    
    Similar to parsing in BioPython, but stripped back into a single module
    '''

    def __init__(self, data):        
        self.iterations = ( Iteration(x) for x in _read('Iteration', data) )
    def __iter__(self):
        return self.iterations
    def __next__(self):
        return next(self.iterations)
    def next(self):
        return self.__next__()

class Iteration:
    def __init__(self, data):
        self.data = data
        self.hits = ( Hit(x) for x in _read('Hit', data) )
        try:
            dct = {
            'num': next(_read('Iteration_iter-num', data)),
            'q_id': next(_read('Iteration_query-ID', data)),
            'q_def': next(_read('Iteration_query-def', data)),
            'q_len': next(_read('Iteration_query-len', data))
            }
            self.__dict__.update(dct)
        except StopIteration:
            raise StopIteration('Iter Error')
    def __iter__(self):
        return self.hits
    def __next__(self):
        return next(self.hits)
    def next(self):
        return self.__next__()
   
class Hit:
    def __init__(self, data):
        self.data = data
        self.hsps = ( Hsp(x) for x in _read('Hsp', data) )
        try:
            dct = {
            'num': next(_read('Hit_num', data)),
            'h_id': next(_read('Hit_id', data)),
            'h_def': next(_read('Hit_def', data)),
            'h_acc': next(_read('Hit_accession', data)),
            'h_len': next(_read('Hit_len', data))
            }
            self.__dict__.update(dct)
        except StopIteration:
            raise StopIteration('Hit Error')
    def __iter__(self):
        return self.hsps
    def __next__(self):
        return next(self.hsps)
    def next(self):
        return self.__next__()

class Hsp:
    def __init__(self, data):
        self.data = data
        try:
            dct = {
            'num': next(_read('Hsp_num', data)),
            'bit': next(_read('Hsp_bit-score', data)),
            'score': next(_read('Hsp_score', data)),
            'expect': next(_read('Hsp_evalue', data)),
            'q_from': next(_read('Hsp_query-from', data)),
            'q_to': next(_read('Hsp_query-to', data)),
            'h_from': next(_read('Hsp_hit-from', data)),
            'h_to': next(_read('Hsp_hit-to', data)),
            'q_frame': next(_read('Hsp_query-frame', data)),
            'h_frame': next(_read('Hsp_hit-frame', data)),
            'identity': next(_read('Hsp_identity', data)),
            'positive': next(_read('Hsp_positive', data)),
            'gaps': next(_read('Hsp_gaps', data)),
            'align_len': next(_read('Hsp_align-len', data)),
            'q_seq': next(_read('Hsp_qseq', data)),
            'h_seq': next(_read('Hsp_hseq', data)),
            'midline': next(_read('Hsp_midline', data))
            }
            self.__dict__.update(dct)
        except StopIteration:
            raise StopIteration('HSP Error')

def _read(tag, search_space):
    '''
    Gets data between tags - returns list for all tag instances
    '''
    for i, d in enumerate(re.split( ('<%s>|</%s>' % (tag, tag)), search_space)):
        x = d.replace(' ', '').strip()
        if d == ' ' or x == '\n' or i == 0 or x == '\r' or x == '':
            continue
        else:
            yield d

def _get(xml, num_ret=0, join_hsps=False, expect=1, r_a=False, **kw):
    '''
    xml - file path to xml output
    num_ret - number of hits to return for each iteration, default = 0 -> returns all
    expect - max e-value threashold for hsp, default = 1
    '''
    
    def weight_avg(v, w):
        return sum([ i*j for i, j in zip(v, w) ])/float(sum(w))
    
    def hits(iteration):
        for i, hit in enumerate(iteration):
            if i < num_ret or num_ret == 0:
                h0 = next(hit)
                if float(h0.expect) <= expect:
                    r = {}
                    
                    r['q_len'] = iteration.q_len
                    r['q_def'] = iteration.q_def
                    r['h_def'] = '>' + hit.h_acc + ' ' + hit.h_def
                    r['h_acc'] = hit.h_acc
                    r['h_len'] = hit.h_len
                    
                    if not join_hsps:
                        hsp = h0    # get best hsp
                        # add hsp data
                        if not r_a:
                            r['score'] = hsp.score
                            r['expect'] = hsp.expect
                            r['psim'] = 100 * float(hsp.positive)/float(hsp.align_len)
                            r['pid'] = 100 * float(hsp.identity)/float(hsp.align_len)
                            r['gaps'] = hsp.gaps
                            r['cover'] = 100 * float(hsp.align_len)/float(iteration.q_len)
                            
                            r['q_from'] = hsp.q_from
                            r['q_to'] = hsp.q_to
                            r['h_from'] = hsp.h_from
                            r['h_to'] = hsp.h_to
                            r['align_len'] = hsp.align_len
                            
                            r['idnt'] = hsp.identity
                            r['pos'] = hsp.positive
                            
                            r['q_seq'] = hsp.q_seq
                            r['h_seq'] = hsp.h_seq
                            r['mid'] = hsp.midline
                        else:
                            
                            r['score'] = []
                            r['expect'] = []
                            r['psim'] = []
                            r['pid'] = []
                            r['gaps'] = []
                            r['cover'] = []
                            r['q_from'] = []
                            r['q_to'] = []
                            r['h_from'] = []
                            r['h_to'] = []
                            r['align_len'] = []
                            r['q_seq'] = []
                            r['h_seq'] = []
                            r['mid'] = []
                            
                            r['idnt'] = []
                            r['pos'] = []
                            
                            while True:
                                # iterate over hsps
                                if float(hsp.expect) <= expect:
                                    r['score'].append(float(hsp.score))
                                    r['expect'].append(float(hsp.expect))
                                    r['psim'].append( 100 * float(hsp.positive)/float(hsp.align_len) )
                                    r['pid'].append( 100 * float(hsp.identity)/float(hsp.align_len ) )
                                    r['gaps'].append(int(hsp.gaps))
                                    r['cover'].append( 100 * float(hsp.align_len)/float(iteration.q_len ) )
                                    
                                    # keep alignments as lists
                                    r['q_from'].append( hsp.q_from )
                                    r['q_to'].append( hsp.q_to )
                                    r['h_from'].append( hsp.h_from )
                                    r['h_to'].append( hsp.h_to )
                                    r['align_len'].append( hsp.align_len )
                                    r['idnt'].append( hsp.identity )
                                    r['pos'].append( hsp.positive )
                                    r['q_seq'].append( hsp.q_seq )
                                    r['h_seq'].append( hsp.h_seq )
                                    r['mid'].append( hsp.midline )
                                try:
                                    hsp = next(hit)     # set to next hsp
                                except StopIteration:
                                    break   
                            
                    else:                        
                        r['score'] = []     # sum
                        r['expect'] = []    # weighted avg
                        r['psim'] = []      # weighted avg
                        r['pid'] = []       # weighted avg
                        r['gaps'] = []      # sum
                        r['cover'] = []     # sum
                        r['q_from'] = []
                        r['q_to'] = []
                        r['h_from'] = []
                        r['h_to'] = []
                        r['align_len'] = []
                        r['q_seq'] = []
                        r['h_seq'] = []
                        r['mid'] = []
                        
                        r['idnt'] = []
                        r['pos'] = []
                        
                        hsp = h0
                        while True:
                            # iterate over hsps
                            if float(hsp.expect) <= expect:
                                r['score'].append(float(hsp.score))
                                r['expect'].append(float(hsp.expect))
                                r['psim'].append( float(hsp.positive)/float(hsp.align_len) )
                                r['pid'].append( float(hsp.identity)/float(hsp.align_len ) )
                                r['gaps'].append(int(hsp.gaps))
                                r['cover'].append( float(hsp.align_len)/float(iteration.q_len ) )
                                
                                # keep alignments as lists
                                r['q_from'].append( hsp.q_from )
                                r['q_to'].append( hsp.q_to )
                                r['h_from'].append( hsp.h_from )
                                r['h_to'].append( hsp.h_to )
                                r['align_len'].append( hsp.align_len )
                                r['idnt'].append( hsp.identity )
                                r['pos'].append( hsp.positive )
                                r['q_seq'].append( hsp.q_seq )
                                r['h_seq'].append( hsp.h_seq )
                                r['mid'].append( hsp.midline )
                            try:
                                hsp = next(hit)     # set to next hsp
                            except StopIteration:
                                break
                        
                        # convert lists to sums and averages    
                        r['score'] = sum(r['score'])
                        r['psim'] = 100 * weight_avg(r['psim'], r['cover']) 
                        r['pid'] = 100 * weight_avg(r['pid'], r['cover'])
                        r['expect'] = 100 * weight_avg(r['expect'], r['cover'])                                                          
                        r['gaps'] = sum(r['gaps'])
                        r['cover'] = 100 * sum(r['cover'])      # should not exceed 100%
                    
                            
                    yield _make_class(r)
    
    record = BlastOutput(xml)
    for iteration in record:
        yield hits(iteration)
        

def _format(s1, s2=None, m=None, num_pl = 60, pad=0, det=False, 
            qf=0, qt=0, hf=0, ht=0, **kw):
    
    # chop seqs into lines of len 60 with 2 white spaces either side
    l1 = []
    l2 = []
    l3 = []
    
    num_pad = len(str(max((qf, qt, hf, ht)))) + 1
    
    qx, qy, hx, hy = 0, int(qf)-1, 0, int(hf)-1
    
    for i in range(0, len(s1), num_pl):
        
        l1_s1 = s1[ i : i + num_pl ]
        l1_line =  l1_s1.ljust(num_pl + pad).rjust(num_pl + 2*pad) 
        
        if det:
            qx = qy + 1
            qy = qy + len(l1_s1) - l1_s1.count('-')
            l1_line = str(qx).ljust(num_pad) + l1_line + str(qy).ljust(num_pad)
        
        l1.append(l1_line)
        
        if s2 is not None and m is not None:
            
            l2_m  =  m[ i : i + num_pl ]
            l3_s2 = s2[ i : i + num_pl ]
            l2_line =  l2_m.ljust(num_pl + pad).rjust(num_pl + 2*pad) 
            l3_line = l3_s2.ljust(num_pl + pad).rjust(num_pl + 2*pad)
            
            if det:
                hx = hy + 1
                hy = hy + len(l3_s2) - l3_s2.count('-')
                l2_line = ''.ljust(num_pad) + l2_line + ''.ljust(num_pad)
                l3_line = str(hx).ljust(num_pad) + l3_line + str(hy).ljust(num_pad)
            
            l2.append(l2_line)
            l3.append(l3_line)
            
        elif s2 is not None:
            l2_s2 = s2[ i : i + num_pl ]
            l2_line = l2_s2.ljust(num_pl + pad).rjust(num_pl + 2*pad)
            
            if det:
                hx = hy + 1
                hy = hy + len(l2_s2) - l2_s2.count('-')
                l2_line = str(hx).ljust(num_pad) + l2_line + str(hy).ljust(num_pad)
            
            l2.append(l2_line)
     
            
    if not l2:
        return '\n'.join(l1)
    elif not l3:
        out = []
        for i, j in zip(l1, l2):
            out.append('\n'.join((i, j)) + '\n')
        return '\n'.join(out)
    else:
        out = []
        for i, j, k in zip(l1, l2, l3):
            out.append('\n'.join((i, j, k)) + '\n')
        return '\n'.join(out)
            

def fasta(xml, raw_data=False, *args, **kw):
    '''
    Returns hits in fasta format
    '''
    
    if raw_data:
        data = xml
    else:
        f = open(xml, 'r+')
        data = mmap.mmap(f.fileno(), 0)
    
    
    records = _get(data, *args, **kw)

    out = ''
    
    for i, r in enumerate(records):
        if i > 0:
            out += '\n'
        h = next(r)
        out += 'Query: >' + h.q_def + '\n\n'
        while True:
            out += '>' + h.h_def + '\n'
            out += _format(h.h_seq, **kw) + '\n'
            try:
                h = next(r)
            except StopIteration:
                break
    
    if not raw_data:
        data.close()
                
    return out

def clustal(xml, raw_data=False, *args, **kw):
    '''
    Returns alignments in clustal style
    '''
    
    if raw_data:
        data = xml
    else:
        f = open(xml, 'r+')
        data = mmap.mmap(f.fileno(), 0)
    
    records = _get(data, *args, **kw)
    
    out = ''
    
    if 'details' in kw.keys():
        det = kw['details']
    elif 'det' in kw.keys():
        det = kw['det']
        del kw['det']
    else:
        det = True
    if 'pad' in kw.keys():
        pad = kw['pad']
        del kw['pad']
    else:
        pad = 1
    
    for i, r in enumerate(records):
        h = next(r)
        if i > 0:
            out += '\n'
        out += 'Query: >' + h.q_def + '\n\n'
        while True:
            out += '>' + h.h_def + '\n'
            # alignment info
            out += 'Score: %s\nCoverage: %.1f%%\nIdentity: %.1f%%\nSimilarity: %.1f%%\nGaps: %s\nE-value: %s\n\n' % (h.score, 
                                                                                                                     h.cover, 
                                                                                                                     h.pid, 
                                                                                                                     h.psim, 
                                                                                                                     h.gaps, 
                                                                                                                     h.expect)
            out += _format(h.h_seq, s2=h.q_seq, m=h.mid, pad=pad, qf=h.q_from, qt=h.q_to,
                           hf=h.h_from, ht=h.h_to, det=det, **kw) + '\n'
            try:
                h = next(r)
            except StopIteration:
                break
    
    if not raw_data:
        data.close()
        
    return out

def details(xml, raw_data=False, *args, **kw):
    '''
    Returns alignment details as list of classes
    '''
    
    if raw_data:
        data = xml
    else:
        f = open(xml, 'r+')
        data = f.read()
        f.close()
        # for big files:
        #data = mmap.mmap(f.fileno(), 0)
    
    records = _get(data, *args, **kw) # returns deeply-nested generator

    out = []
    
    for r in records:
        out.append([])
        for h in r:
            out[-1].append(h)
    
    if not raw_data:
        pass#f.close()
        
    return out

if __name__ == '__main__':
    print(clustal('temp_blast.xml', num_ret=0))
