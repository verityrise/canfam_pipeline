#!/usr/bin/python
'''
This programs is to integrate dog reference genome from chr to a single one.
Author: Hongzhi Luo
'''
import gzip
import glob
import shutil

path='/vlsci/LSC0007/shared/canine_alport_syndrome/ref_files/'
#path=''
prefix='cfa_ref_CanFam3.1'

def integrate_genome():
    '''
    @param: num: chr1...chrMT in the list.
    '''
    files=glob.glob(path+'*_CanFam3.1_*.fa.gz')
    print files
    #cat_together=[]
    for f in files:
        #cat_together.append(f)
    #for files in cat_together:
        outfile=gzip.open(path+prefix+".fa.gz",'wb')
        for f in files:
            gfile=gzip.open(f)
            outfile.write(gfile.read())
            gfile.close()
        outfile.close()
    
if __name__ == '__main__':
    integrate_genome()