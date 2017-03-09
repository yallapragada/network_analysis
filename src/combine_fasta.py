from Bio import SeqIO
import sys, os

#combine per-protein fasta files from folder1 and folder2
#write output to outfolder
def combine_fasta(infolder1, infolder2, outfolder):
   
    afile1 = infolder1 + os.sep + 'ha.fasta'
    afile2 = infolder1 + os.sep + 'na.fasta'
    afile3 = infolder1 + os.sep + 'm1.fasta'
    afile4 = infolder1 + os.sep + 'm2.fasta'
    afile5 = infolder1 + os.sep + 'np.fasta'
    afile6 = infolder1 + os.sep + 'pb1.fasta'
    afile7 = infolder1 + os.sep + 'pb2.fasta'
    afile8 = infolder1 + os.sep + 'pa.fasta'
    afile9 = infolder1 + os.sep + 'ns1.fasta'
    afile10 = infolder1 + os.sep + 'ns2.fasta'

    bfile1 = infolder2 + os.sep + 'ha.fasta'
    bfile2 = infolder2 + os.sep + 'na.fasta'
    bfile3 = infolder2 + os.sep + 'm1.fasta'
    bfile4 = infolder2 + os.sep + 'm2.fasta'
    bfile5 = infolder2 + os.sep + 'np.fasta'
    bfile6 = infolder2 + os.sep + 'pb1.fasta'
    bfile7 = infolder2 + os.sep + 'pb2.fasta'
    bfile8 = infolder2 + os.sep + 'pa.fasta'
    bfile9 = infolder2 + os.sep + 'ns1.fasta'
    bfile10 = infolder2 + os.sep + 'ns2.fasta'

    cfile1 = outfolder + os.sep + 'ha.fasta'
    cfile2 = outfolder + os.sep + 'na.fasta'
    cfile3 = outfolder + os.sep + 'm1.fasta'
    cfile4 = outfolder + os.sep + 'm2.fasta'
    cfile5 = outfolder + os.sep + 'np.fasta'
    cfile6 = outfolder + os.sep + 'pb1.fasta'
    cfile7 = outfolder + os.sep + 'pb2.fasta'
    cfile8 = outfolder + os.sep + 'pa.fasta'
    cfile9 = outfolder + os.sep + 'ns1.fasta'
    cfile10 = outfolder + os.sep + 'ns2.fasta'
    
    outfile = open(cfile1, "w")
    infile1 = open(afile1, "r")
    infile2 = open(bfile1, "r")
    outfile.write(infile1.read())
    outfile.write(infile2.read())
    outfile.close()
    infile1.close()
    infile2.close()

    outfile = open(cfile2, "w")
    infile1 = open(afile2, "r")
    infile2 = open(bfile2, "r")
    outfile.write(infile1.read())
    outfile.write(infile2.read())
    outfile.close()
    infile1.close()
    infile2.close()

    outfile = open(cfile3, "w")
    infile1 = open(afile3, "r")
    infile2 = open(bfile3, "r")
    outfile.write(infile1.read())
    outfile.write(infile2.read())
    outfile.close()
    infile1.close()
    infile2.close()

    outfile = open(cfile4, "w")
    infile1 = open(afile4, "r")
    infile2 = open(bfile4, "r")
    outfile.write(infile1.read())
    outfile.write(infile2.read())
    outfile.close()
    infile1.close()
    infile2.close()

    outfile = open(cfile5, "w")
    infile1 = open(afile5, "r")
    infile2 = open(bfile5, "r")
    outfile.write(infile1.read())
    outfile.write(infile2.read())
    outfile.close()
    infile1.close()
    infile2.close()

    outfile = open(cfile6, "w")
    infile1 = open(afile6, "r")
    infile2 = open(bfile6, "r")
    outfile.write(infile1.read())
    outfile.write(infile2.read())
    outfile.close()
    infile1.close()
    infile2.close()

    outfile = open(cfile7, "w")
    infile1 = open(afile7, "r")
    infile2 = open(bfile7, "r")
    outfile.write(infile1.read())
    outfile.write(infile2.read())
    outfile.close()
    infile1.close()
    infile2.close()

    outfile = open(cfile8, "w")
    infile1 = open(afile8, "r")
    infile2 = open(bfile8, "r")
    outfile.write(infile1.read())
    outfile.write(infile2.read())
    outfile.close()
    infile1.close()
    infile2.close()

    outfile = open(cfile9, "w")
    infile1 = open(afile9, "r")
    infile2 = open(bfile9, "r")
    outfile.write(infile1.read())
    outfile.write(infile2.read())
    outfile.close()
    infile1.close()
    infile2.close()

    outfile = open(cfile10, "w")
    infile1 = open(afile10, "r")
    infile2 = open(bfile10, "r")
    outfile.write(infile1.read())
    outfile.write(infile2.read())
    outfile.close()
    infile1.close()
    infile2.close()

if __name__ == '__main__':
    infolder1     = sys.argv[1]
    infolder2     = sys.argv[2]
    outfolder     = sys.argv[3]
    combine_fasta(infolder1=infolder1, infolder2=infolder2, outfolder=outfolder)