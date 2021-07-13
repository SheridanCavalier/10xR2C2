#!/usr/bin/env python3
# Roger Volden

'''
Makes a sorted histogram of how many reads there are per
barcode. I'll determine where to have the cutoff for the
barcodes, and then save those separately.

Usage:
    python3 detBarcodes.py R2C2_10x_postprocessed.fasta 737K-august-2016.txt >1500_most_frequent_bcs.fasta

R2C2_10x_postprocessed.fasta contains the cell barcode information
737K-august-2016.txt contains the valid barcodes from 10x
'''

import sys

def fastaReader(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            readDict[line[1:]] = ''
            lastHead = line[1:]
        else:
            readDict[lastHead] += line
    return readDict

def readBarcodes(inFile):
    '''Reads in the canonical 10x barcodes'''
    barcodes = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        barcodes[line] = 0
    return barcodes

def tallyReads(reads, barcodes):
    '''Counts the number of occurrences of each barcode'''
    totalmatched=0
    total=0
    for read in reads:
        total += 1
        cellBC = read[:16]
        if cellBC in barcodes:
            barcodes[cellBC] += 1
            totalmatched +=1
    #print('%s out of %s total reads had exact matches to canonical barcodes'%(totalmatched, total))
    return barcodes,totalmatched

def makeHist(counts):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mplpatches
    # plt.style.use('BME163')

    plt.figure(figsize = (6, 2))
    hist = plt.axes([0.1, 0.1, 8/10, 4/5], frameon = True)

    for i in range(2000):
        width, height = 1, counts[i][0]
        bar = mplpatches.Rectangle((i, 0), width, height,
                                   lw=0, fc='black')
        hist.add_patch(bar)

    hist.set_xlim(0, 1800)
    hist.set_ylim(0, 1500)
    hist.set_title('Nanopore barcode counts')
    hist.set_xlabel('Barcode')
    hist.set_ylabel('Count')
    hist.tick_params(axis='both',which='both',\
                     bottom='on', labelbottom='on',\
                     left='on', labelleft='on',\
                     right='off', labelright='off',\
                     top='off', labeltop='off')
    plt.savefig('combined_ro_hist.png', bbox_inches="tight")

def main():
    readDict = fastaReader(sys.argv[1])
    reads = list(readDict.values())
    barcodes = readBarcodes(sys.argv[2])
    barCounts, totalmatched = tallyReads(reads, barcodes)
    counts = []
    for bc, count in barCounts.items():
        counts.append((count, bc))
    counts = sorted(counts, reverse=True)

    n = 0
    toptotal=0
    for count in counts[:2500]: #manually change N most frequent barcodes here
        print('>barcode_' + str(n) + '_' + str(count[0]))
        print(count[1])
        toptotal += count[0]
        n += 1
    #makeHist(counts)
    #percentage=(toptotal/totalmatched)*100
    #print('Of reads with exact matches to canonical barcodes, reads with the top x most frequent barcodes account for %s percent' %(percentage)) 
if __name__ == '__main__':
    main()
