import os
import sys

def main(fa_fpath, output_fpath):
    with open(fa_fpath, 'r') as infh, open(output_fpath, 'w') as outfh:
        for line in infh:
            if line[0] == '>':
                vids = line[1:-1].split(' ')
                vids.insert(0, vids[0])
                outfh.write('{}\t{}\n'.format(vids[0], ','.join(vids[1:])))
    



if __name__ == '__main__':
    
    fa_fpath = sys.argv[1]
    output_fpath = sys.argv[2]
    
    main(fa_fpath, output_fpath)
    

