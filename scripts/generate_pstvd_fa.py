import os
import sys



def read_fasta(fpath):
    vdict = {}
    sid  = ''
    sseq = ''
    with open(fpath, 'r') as infh:
        for fline in infh:
            fline = fline.rstrip()
            if fline[0] == '>':
                if sid != '':
                    vdict[sid] = sseq
                sid = fline[1:]
                sseq = ''
            else:
                sseq += fline
    vdict[sid] = sseq
    return vdict



def filter_PSTVd_seq(vdict):
    rm_list = []
    for sid in vdict.keys():
        if not (('potato spindle tuber viroid' in sid.lower()) and ('complete' in sid.lower())):
            rm_list.append(sid)
    for sid in rm_list:
        del vdict[sid]
    return vdict
    


def filter_abnormal_seq(vdict):
    rm_list = []
    for sid in vdict.keys():
        if not (340 <= len(vdict[sid]) and len(vdict[sid]) <= 365):
            rm_list.append(sid)
    for sid in rm_list:
        del vdict[sid]
    return vdict


   
def filter_uniq_seq(vdict):
    seq2id = {}
    for sid, sseq in vdict.items():
        if sseq not in seq2id:
            seq2id[sseq] = []
        seq2id[sseq].append(sid.split(' ')[0])
    uniq_vdict = {}
    for sseq, sids in seq2id.items():
        sid = ' '.join(sorted(sids))
        uniq_vdict[sid] = sseq
    return uniq_vdict



def write_fasta(vdict, fpath):
    with open(fpath, 'w') as outfh:
        for sid in sorted(vdict.keys()):
            outfh.write('>' + sid + '\n')
            outfh.write(vdict[sid]+ '\n')
    
    
    

if __name__ == '__main__':
    infa_fpath = sys.argv[1]
    outfa_fpath = sys.argv[2]
    
    vdict = read_fasta(infa_fpath)
    print('#sequence           ', len(vdict))

    vdict = filter_PSTVd_seq(vdict)
    print('#sequence (PSTVd)   ', len(vdict))

    vdict  = filter_abnormal_seq(vdict)
    print('#sequence (abnormal)', len(vdict))
    
    vdict = filter_uniq_seq(vdict)
    print('#sequence (unique)  ', len(vdict))

    write_fasta(vdict, outfa_fpath)


