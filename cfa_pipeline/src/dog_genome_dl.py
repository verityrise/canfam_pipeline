#!/usr/bin/python
'''
This program uses urlretrive function under urllib to download canine genome data.
Author: Hongzhi Luo
'''

import urllib

#defined address of ftp server, file prefix, and target path
canine_genome_url='ftp://ftp.ncbi.nih.gov/genomes/Canis_lupus_familiaris/Assembled_chromosomes/seq/'
prefix='cfa_ref_CanFam3.1_'
path='/vlsci/LSC0007/shared/canine_alport_syndrome/ref_files/'

def dl_dog_genome():
    '''
    @param: num: chr1...chrMT in the list.
    '''
    chr_num=['chr' + str(n) for n in range(1, 39)]
    chr_num.append('chrX')
    chr_num.append('chrMT')
    chr_num.append('unplaced')
    for num in chr_num:
        print 'downloading...',num
        #download the fasta file with .fa extension
        urllib.urlretrieve(canine_genome_url+prefix+"{0}.fa.gz".format(num), path+prefix+"{0}.fa.gz".format(num))
    print 'download finished'

if __name__ == '__main__':
    dl_dog_genome()