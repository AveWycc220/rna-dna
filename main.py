import eel
import os
import skbio.io.format
import skbio.io
import skbio.io.format.blast6
import skbio.io.format.blast7
import skbio.io.format.clustal
import skbio.io.format.embl
import skbio.io.format.emptyfile
import skbio.io.format.fasta
import skbio.io.format.fastq
import skbio.io.format.genbank
import skbio.io.format.gff3
import skbio.io.format.lsmat
import skbio.io.format.newick
import skbio.io.format.ordination
import skbio.io.format.phylip
import skbio.io.format.qseq
import skbio.io.format.stockholm
import skbio.io.format.tests
import skbio.io.format._base
import skbio.io.format._blast
import skbio.io.format._sequence_feature_vocabulary
from middle.controler import *

""" CONST """
PATH = os.path.dirname(os.path.abspath(__file__)) + '\\'

if __name__ == '__main__':
    eel.init(PATH)
    eel.start(PATH + 'front\\index.html', mode='chrome', port=0, size=(1000, 720))
