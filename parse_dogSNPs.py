
import pandas as pd
import os
import urllib

##########  hard coded parameters
LOCATION = {
    "HZL-PC": os.path.normpath("D:/EDUCATION/BINF90007/ResearchProj/dogs/ref_snp/"),
}

path="E:/EDUCATION/BINF90007/ResearchProj/dogs/ref_snp/"
snp_url = "ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/dog_9615/rs_fasta/"
chroms = ['ch' + str(n) for n in range(1, 39)]
chroms.append('chX')
chroms.append('chMulti')
chroms.append('chNotOn')
chroms.append('chUn')

def dl_snp():
    for ch in chroms:
        print 'downloading...    ', ch
        urllib.urlretrieve(
            snp_url+"rs_{0}.fas.gz".format(ch),
            os.path.join(path,"rs_{0}.fas.gz".format(ch)))
    print 'download finished~~~~~~~~~~~~~~~'

if __name__ == '__main__':
    print chroms
    print path
    dl_snp()