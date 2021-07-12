#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_fasta', type=str)
parser.add_argument('-o', '--output_path', type=str)
parser.add_argument('-x', '--splint_number', type=str)
parser.add_argument('-i5', '--input5_fasta', type=str)
parser.add_argument('-i3', '--input3_fasta', type=str)
args = parser.parse_args()
output_path = args.output_path + '/'
input_reads = args.input_fasta
splint = args.splint_number
input5 = args.input5_fasta
input3 = args.input3_fasta

def read_fasta(infile):
    reads = {}
    sequence = ''
    for line in open(infile):
      if line:
        a = line.strip()
        if len(a) > 0:
          if a[0] == '>':
            if sequence != '':
                reads[name] = sequence
            name = a[1:].split()[0]
            sequence = ''
          else:
            sequence += a
    if sequence != '':
        reads[name] = sequence
    return reads

def getFlanking(splint):
    flankingDict = {}
    flankingDict['1'] = ('CCATA', 'ATCAC', 'ATCGC', 'TTAGT')
    flankingDict['2'] = ('AATAA', 'ATTGC', 'TCACG', 'CTGCC')
    flankingDict['3.5'] = ('TGGGT', 'TAAAA', 'CAGCT', 'ATATT')
    flankingDict['4.5'] = ('TCCGT', 'TACGA', 'AGGCG', 'ACCTG')
    flankingDict['5'] = ('GATAG', 'TTCTG', 'AAGAG', 'GAACC')
    flankingDict['6'] = ('CTGGT', 'ACGTC', 'ATTAG', 'TAGCC')
    flankingDict['7'] = ('TTATA', 'TGGAC', 'AGAGG', 'CAGAT')
    flankingDict['8'] = ('AGAAT', 'TTCTC', 'AGGCT', 'AGGGA')
    flankingDict['9'] = ('CGGGG', 'AAGAT', 'GTAAC', 'ACCTA')
    flankingDict['10'] = ('TTCAG', 'ATTAA', 'CTGGT', 'AGCCT')
    flankingDict['11'] = ('CGATA', 'ATTCT', 'TGGTG', 'TCATA')
    flankingDict['12'] = ('CAAGT', 'CATAC', 'AAGTC', 'CACAA')
    flankingDict['13'] = ('ATCTG', 'ATGCA', 'ATGAC', 'CGTAG')
    flankingDict['14'] = ('TGGAC', 'ACTTC', 'TGAGG', 'TTAAT')
    flankingDict['15'] = ('CCTGT', 'CATGG', 'TCGTA', 'AATTT')
    flankingDict['16'] = ('GTACA', 'TCTAT', 'AGAAT', 'TTTCT')

    return flankingDict[splint]

def reverse_complement(sequence):
  Seq = ''
  complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '-':'-'}
  for item in sequence[::-1]:
    Seq = Seq + complement[item]
  return Seq

def find_UMI(sequence, frame1, frame2):
    UMI = ''
    for i in range(0, max(len(sequence)-25, 1), 1):
        match1 = sequence[i:i+5]
        if match1 == frame1:

            for j in range(i+12, i+25, 1):
                match2 = sequence[j:j+5]
                if match2 == frame2:
                    UMI = sequence[i+5:j]
    return UMI




flanking1, flanking2, flanking3, flanking4 = getFlanking(splint)

reads5 = read_fasta(input5)
reads3 = read_fasta(input3)
reads = read_fasta(input_reads)

out = open(output_path + '/R2C2_Consensus.UMI', 'w')

for read in reads:
    name=read
    total=total+[1]
    UMI5 = ''
    UMI3 = ''
    if reads5['5Prime_adapter'] in reads[read]:
        sequence5=reads[read]
        UMI5=find_UMI(sequence5, reverse_complement(flanking4), reverse_complement(flanking3))
    elif reverse_complement(reads5['5Prime_adapter']) in reads[read]:
        sequence5=reads[read]
        UMI5=find_UMI(sequence5, flanking3, flanking4)
    else:
        sequence5 = ''
        UMI5=''
    if reads3['3Prime_adapter'] in reads[read]:
        sequence3=reads[read]
        UMI3 = find_UMI(sequence3, reverse_complement(flanking2), reverse_complement(flanking1))
    elif reverse_complement(reads3['3Prime_adapter']) in reads[read]:
        sequence3 = reads[read]
        UMI3 = find_UMI(sequence3, flanking1, flanking2)
    else:
        sequence3 = ''
        UMI3=''
    #if UMI3==''or UMI5=='':
        #excluded=excluded+[1]
    out.write('%s\t%s\t%s\n' %(name, UMI5, UMI3))

