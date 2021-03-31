#!/usr/bin/env python3
import argparse
import editdistance
import numpy as np
import os
import sys
import mappy as mm #mappy is an interface for minimap2, importing as mm
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta_file', type=str)
parser.add_argument('-s', '--subreads_file', type=str)
parser.add_argument('-o', '--output_path', type=str)
parser.add_argument('-u', '--umi_file', type=str)
parser.add_argument('-c', '--config_file', type=str)
parser.add_argument('-m', '--score_matrix', type=str)

args = parser.parse_args()
path = args.output_path + '/'
fasta_file = args.fasta_file
subreads_file = args.subreads_file

umi_file = args.umi_file
config_file= args.config_file
score_matrix = args.score_matrix
subsample = 200

def configReader(configIn): #configIn is our config file where each line has the name of the program, a tab separation, and the path to the executable
    '''Parses the config file.'''
    progs = {} #defines progs as an empty dictionary
    for line in open(configIn): #for each line in the config file
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t') #splits each line by tab, creating a new list called "line" that has the name at index 0 and path at index 1
        progs[line[0]] = line[1] #appends new value to empty progs dictionary where the key is the program name, and the value stored under the key is the path like {'poa':'/usr/bin/poaV2/poa'}
        #Now we have a progs dictionary where the path to each program is stored under the program name
    possible = set(['poa', 'minimap2', 'gonk', 'consensus', 'racon', 'blat','emtrey', 'psl2pslx']) #creates set of possible program names; do I need emtrey and psl2pslx?
    inConfig = set() #creates an empty set called inConfig
    for key in progs.keys(): #for each program name in the config file, add the program name to inConfig set
        inConfig.add(key)
        if key not in possible: #if the program name from the program file is not in the listof possible programs, raise exception
            raise Exception('Check config file')
    # check for missing programs
    # if missing, default to path
    for missing in possible-inConfig: #create a list that is the subtraction of programs in the config file from possible programs, for each of the elements in this list...
        if missing == 'consensus': #if the program that is missing from inConfig is the consensus program...
            path = 'consensus.py' #the path is just the local consensus.py file
        else:
            path = missing #otherwise path is equal to the name of the missing program 
        progs[missing] = path 
        sys.stderr.write('Using ' + str(missing)
                         + ' from your path, not the config file.\n')
    return progs

progs = configReader(config_file)
poa = progs['poa'] # for an empty config file, the key 'poa' doesn't store the path, it just stores the name 'poa'. Which is fine because it's already an executable, you can just type 'poa' into the terminal and it works.
minimap2 = progs['minimap2']
racon = progs['racon']
consensus = progs['consensus']
#looks like the only programs needed are the consensus.py script, poa-biopipeline, minimap2 and racon (I got emtrey but can't find psl2pslx anywhere)

def determine_consensus(name, fasta, fastq, temp_folder): #defining function 'determine_consensus' with four arguments 
    '''Aligns and returns the consensus''' #so if alignment occurs does it need minimap2? Presumably this is the same alignment and consensus strategy that C3POa uses but the C3POa aligner is conk...
    corrected_consensus = '' #defines corrected_consensus variable as an empty string
    out_F = fasta #out_F is the fasta argument input  
    fastq_reads = read_fastq_file(fastq, False) #read_fastq_file function determined futher down; it's inputs are (seq file, and check).  
    out_Fq = temp_folder + '/subsampled.fastq' #out_Fq is the temp_folder input/subsampled.fastq (this creates a directory to store the out fq)
    out = open(out_Fq, 'w') #write out to the /subsampled.fastq file defined above
    indexes = np.random.choice(np.arange(0, len(fastq_reads), 1), min(len(fastq_reads), subsample), replace=False) #uses NumPy function np.random.choice 
    #np.arange takes arguments in (start, stop, step) form and returns an array with evenly-spaced values. Here the start value is 0, the stop value is the length of the fastq_reads file (so maybe how many reads?), and 1 is the interval
    #then np.random.choice takes the following inputs (a, size, replace, p) where a is an array (in this case it is the np.arange array), size is the size of the returned random sample
    #here that is the smallest of the two values: fastq_reads length or subsample (defined as 200 above), and replace=False (still don't really know what this means)
    
    #Basically this is saying give me 200 random values between 0 and the total number of reads in the fastq, these values are then the positions of 200 random reads in the fastq input file that will be used for the rest of this function

    subsample_fastq_reads = [] #define subasmple_fastq_reads as an empty list
    for index in indexes: #for each value in indexes (that list of 200 random numbers)...
        subsample_fastq_reads.append(fastq_reads[index]) #add the random read at that position in the fastq file to this empty subsample_fastq_reads list
    #At the end of this 'for' loop we now have a populated list of 200 randomly-selected reads from our fastq file
    for read in subsample_fastq_reads: #now for each of these randomly-selected reads, write out to our fastq_out directory...
        out.write('@' + read[0] + '\n' + read[1] + '\n+\n' + read[2] + '\n') #write out @ (the start of the first line for a fastq) plus the read at position [0] which would be the name??
    out.close()                                                              #then the read sequence [1] then the scores [2]??? Need to investigate what this actually makes...
#at the end of all this you have a temporary fastq file named 'subsampled.fastq' that should have 200 randomly selected reads from your input fastq
    
    poa_cons = temp_folder + '/consensus.fasta' #defines poa_cons as a consensus.fasta file in the output path
    final = temp_folder + '/corrected_consensus.fasta' #defines final as a corrected_consensus.fasta file in the output path
    overlap = temp_folder + '/overlaps.sam' #defines overlap as a overlaps.sam file in the output path
    pairwise = temp_folder + '/prelim_consensus.fasta' #defines pairwise as a prelim_consensus.fasta file in the output path

    max_coverage, repeats = 0, 0 #assigns two variables at once, max_coverage=0 and repeats=0
    reads = read_fasta(out_F) #defines variable 'reads' as the PRODUCT of the read_fasta function operated on out_F (the input fasta to the determine_consensus function)
    qual, raw, before, after = [], [], [], [] #defines four empty list variables qual, raw, before, after

    header_log = path + 'header_associations.tsv' #creates a new log file called header_associations.tsv
    header_fh = open(header_log, 'a+') #open header log file and APPEND (as opposed to 'write' which will add to the beginning of the file, append will add to the end)
    headers = [] #create empty list called 'headers'
    for read in reads: #the out_F file is defined as the fasta that is input into the determine_consensus function, when the determine_consensus function is actually called in
                       #the make_consensus function, the given fasta is something called 'temp_consensus_reads.fasta'; if I look at the structure of that file, the names that are
                       #being referred to by read are structured: 'f11ac981-9408-4643-9b65-a250b0232be2_21.57_30141_25_1284|10'
        info = read.split('_') #splits a string into a list using '_' character as the break, so for the aforementioned name it would be readName_averageQuality_originalReadLength_numberOfRepeats_subreadLength
        coverage = int(info[3]) #coverage is info at index 3, which is 'Number of Repeats'
        headers.append(info[0]) #headers is info at index 0, which is 'Read Name'
        qual.append(float(info[1])) #qual is the decimal value of info at index 1, which is "Average Quality"
        raw.append(int(info[2])) #raw is the integer value of info at index 2, which is "Original Read Length"
        repeats += int(info[3]) #repeats is the integer value of info at index 3, which is "Number of Repeats" now ADDED to the variable 'repeats' defined outside the loop as 0
        before.append(str(info[4])) #python is having trouble here "invalid literal for int() with base 10: '1284|10'"--this is where the error is coming from, the last '_' delimited
                                    #item in the read name '1284|10' so either they messed up or this is not the right fasta file...
                                    #the R2C2_Consensus and R2C2_Subreads files generated from C3POa don't have the |10 on the end...
                                    #'before' is the empty list created earlier in the function; but info[4] is int|10, python thinks we are trying to do a bitwise operation
        after.append(int(info[5].split('|')[0])) #'after' is the empty list created earlier; info[5] shouldn't exist but it is clear that this was written with the int|10 value in mind
                                                 #am I putting the wrong fasta in here? 

        if coverage >= max_coverage: #max_coverage was defined as 0, coverage is the number of repeats 
            best = read #if coverage is greater than 0 then 'best' is the name of the read
            max_coverage = coverage #max_coverage becomes the number of repeats

    print('\t'.join(headers), file=header_fh) 
    header_fh.close()

    out_cons_file = open(poa_cons, 'w')
    out_cons_file.write('>' + best + '\n' + reads[best].replace('-', '') + '\n')
    out_cons_file.close()

    final = poa_cons
    input_cons = poa_cons
    output_cons = poa_cons.replace('.fasta', '_1.fasta')

    os.system('%s --secondary=no -ax map-ont \
              %s %s > %s 2> ./minimap2_messages.txt'
              % (minimap2, input_cons, out_Fq, overlap))
    os.system('%s -q 5 -t 1 \
               %s %s %s >%s 2>./racon_messages.txt'
               %(racon, out_Fq, overlap, input_cons, output_cons))
    final = output_cons

    # print(final)
    reads = read_fasta(final)
    for read in reads:
        corrected_consensus = reads[read]

    return corrected_consensus, repeats, headers[0], round(np.average(qual), 2), int(np.average(raw)), int(np.average(before)), int(np.average(after))

def read_subreads(seq_file, chrom_reads): #define read_subreads function (referenced in processing with the inputs subreads_file, chrom_reads)--this is the subreads input fastq
    for read in mm.fastx_read(seq_file, read_comment=False): #use mappy to read this file line by line
        root_name = read[0].split('_')[0] #this pulls the rootname from the subreads '@ed20f548-6dac-45ac-acf1-076b26267bf9'--make note! This name has @ in front, but chrom_reads does NOT
        if root_name in chrom_reads: #chrom_reads has NO group information...it's literally just ALL of the read names; also this won't work because the '@' in front
            # root_name : [(header, seq, qual), ...]
            chrom_reads[root_name].append(read) # read = (header, seq, qual); append read info in this empty chrom_reads dictionary
    return chrom_reads #so now we have a dictionary that contains the root read name and all of the fastq subread info for that root read; but we have NO group info...

def read_fasta(infile): #defining read_fasta (this function is called in the define_consensus function) 
    reads = {} #generates empty dictionary called 'reads'
    for read in mm.fastx_read(infile, read_comment=False): #mm. (similar to np. syntax) denotes a mappy function; here mm.fastx_read with read_comment=False generates a (name, seq, qual) tuple for each sequence entry
        #tuples are used to store multiple items in one variable and they are indexed
        reads[read[0]] = read[1] #creates a dictionary calleds 'reads' where the sequence is stored under the name of the sequence 
    return reads

def read_fastq_file(seq_file, check): #defining the read_fastq_file with two arguments (seq_file, check)
    '''
    Takes a FASTQ file and returns a list of tuples
    In each tuple:
        name : str, read ID
        seed : int, first occurrence of the splint
        seq : str, sequence
        qual : str, quality line
        average_quals : float, average quality of that line
        seq_length : int, length of the sequence
    Has a check mode where if it sees one read, it'll return True
    '''
    read_list = []
    for read in mm.fastx_read(seq_file, read_comment=False):
        split_name = read[0].split('_')
        name, seed = split_name[0], 0
        seq, qual = read[1], read[2]
        if check:
            return True
        avg_q = np.average([ord(x)-33 for x in qual])
        s_len = len(seq)
        read_list.append((read[0], seq, qual, avg_q, s_len))
    return read_list

def make_consensus(Molecule, UMI_number, subreads): #define make_consensus function with input arguments (Molecule, UMI_number, subreads)
    subread_file = path + '/temp_subreads.fastq' #creates subread_file 'temp_subreads.fastq'
    fastaread_file = path + '/temp_consensus_reads.fasta' #creates fastaread_file 'temp_consensus_reads.fasta'
    subs = open(subread_file, 'w') #open subread_file (temp_subreads.fastq) and write to the beginning of the file
    fasta = open(fastaread_file, 'w') #open fastaread_file (temp_consensus_reads.fasta) and write to the beginning of the file
    for read in Molecule: #make_consensus function is called later in the group_reads function (Molecule variable is defined there); 
        fasta.write(read) 
        # print(read)
        root_name = read[1:].split('_')[0]
        raw = subreads[root_name] # [(header, seq, qual), ...], includes rootname_subread_n
        for entry in raw:
            subs.write('@' + entry[0] + '\n' + entry[1] + '\n+\n' + entry[2] + '\n')
    subs.close()
    fasta.close()
    if read_fastq_file(subread_file, True):
        corrected_consensus, repeats, name, qual, raw, before, after = determine_consensus(str(UMI_number), fastaread_file, subread_file, path) #it's the determine_consensus function from earlier! This function takes the following inputs (name, fasta, fastq, temp_folder)
        #so here the determine_consensus function is being operated on the following inputs, name=UMI_number (even tho name is never actually called in the original function), fastaread_file (which is defined above), subread_file (defined above), and the path
        return '>%s_%s_%s_%s_%s_%s|%s\n%s\n' %(name, str(qual), str(raw), str(repeats), str(before), str(after), str(UMI_number), corrected_consensus)
    else:
        return ''

def parse_reads(reads, sub_reads, UMIs): #defining function parse_reads with inputs (reads, sub_reads, UMIs); from the main function:reads is the read:sequence dictionary from the input R2C2_Consensus.fasta
    #sub_reads is the readname:group dictionary created in the main function (this encompasses all the reads)
    #UMIs is the readname:UMI5,UMI3 dictionary generated from the read_UMIs(umi_file) function 
    group_dict, chrom_reads= {}, {} #defining variables group_dict and chrom_reads as empty dictionaries
    groups = [] #defining groups as an empty list
    UMI_group, previous_start, previous_end = 0, 0, 0 #setting the variables UMI_group, previous_start, and previous_end as 0
    for name, group_number in sub_reads.items(): #the .items() function in python returns the key:value pairs in a dictionary; so it's the name from R2C2_Consensus.fasta (readName_averageQuality_originalReadLength_numberOfRepeats_subreadLength) and then the group number
        root_name = name.split('_')[0] #root_name equals the read name split by '_' and the 0 index position--the name of the read (not qual scores or anything else)
        UMI5 = UMIs[name][0] #UMI dictionary again is readname: UMI5, UMI3 so UMIs[name][0] will be UMI5, and UMIs[name][1] will be UMI3
        UMI3 = UMIs[name][1] #as above
        if root_name not in chrom_reads: #if root_name is not in chrom_reads (which it shouldn't be because these are all unique consensus read names, not subread names...)
            chrom_reads[root_name] = [] #chrom_reads list at this name is empty
            if group_number not in group_dict: #if the group number is not in the group_dict (group number IS something that can be repeated); but for the first time it appears...
                group_dict[group_number] = [] #group_dict at [group number] is empty
            group_dict[group_number].append((name, UMI5, UMI3, reads[name])) #group_dict at the group_number, append the read name, UMI5, UMI3, and read sequence (which we didn't have before associated with the UMI or group#)
    for group_number in sorted(group_dict): #sorts group dictionary by ascending number
        group = group_dict[group_number] 
        groups.append(list(set(group))) #this basically creates a list version of the group_dict, just without the 'key' the group number; this list can be indexed
    return groups, chrom_reads #creates a dictionary root_name:[]; this dictionary is empty because there is no root_name in chrom_reads as these are consensus reads so the root read name will never be the same

def group_reads(groups, reads, subreads, UMIs, final, final_UMI_only, matched_reads): #define function group_reads with inputs (groups, reads, subreads UMIs, final, final_UMI_only, matched_reads)
    #in the main function these inputs are (annotated_groups, reads, subreads, UMIs, final, final_UMI_only, matched_reads)
    #annotated_groups is a list version of group_dict; reads are still grouped together with read_name, UMI5, UMI3, sequence but there is no group number
    #reads are the reads in the R2C2_Consensus.fasta file
    #subreads is the chrom_reads dictionary filled with all of the subread names, sequences, and quality scores
    #UMIs are the UMIs from read_UMIs function
    #final, final_UMI_only, and matched_reads are out files 
   UMI_group = 0 #set variable UMI_group equal to 0 
    # print(len(groups))
    for group in groups: # 
        group = list(set(group))
        UMI_dict = {}
        set_dict = {}
        group_counter = 0 #now each group is it's own set, so for each set...let's carry through with an example where there are multiple reads per group
        if len(group) > 1: #so for each group, if there are more than two reads in this group...
            group_counter += 1 #group_counter goes up by one...but this variable is within the for loop so it will be reset to 0 after each group...
            # print('test', group_counter, len(group))
            UMI_counter = 0 #UMI_counter equals 0 
            for i in range(0, len(group), 1): #let's say our group has two reads, so for i=0 and 1
                UMI_dict[group[i][0]] = set() #creates a UMI_dict that is empty if there is only one read in a group, but for groups with two or more reads it creates a temporary UMI_dict tuple
                                              #where each read name in the group is the key to an empty set; so for a group with read_1 and read_2, UMI_dict is {read_1:set(), read_2:set()}
                                              #I'm assuming this set gets populated 
            # if len(group) == 2:
            #     print(group[0][1], group[0][2])
            #     print(group[1][1], group[1][2])

            for i in range(len(group)): #comparing each read in the group with the next read (last read is compared ot itself??); violations give the read pair the 'single' attribute rather than 'both'
                UMI_counter += 1 
                UMI_dict[group[i][0]].add(UMI_counter) 
                for j in range(i+1, len(group)): 
                    if np.abs(len(group[i][3])-len(group[j][3]))/len(group[i][3]) < 0.1: #if the read pair is about the same length...
                        status = 'both' #status = 'both'
                        if len(group[i][1]) > 0 and len(group[j][1]) > 0: #if the UMI5 for each of the reads in the pair exists...
                            dist5 = editdistance.eval(group[i][1], group[j][1]) #dist5 is the Levenshtein distance between the UMI5s

                        else: #if the UMI5 doesn't exist for one or both reads in the compared pair...
                            dist5 = 15 #then L=15 because UMI5 is 15bp long
                            status = 'single' #status is changed from 'both' to 'single' 
                        if len(group[i][2]) > 0 and len(group[j][2]) > 0: #NOW is the UMI3 for both reads exists...
                            dist3 = editdistance.eval(group[i][2], group[j][2]) #probe the Levenshtein distance between them 
                        else: #if the read pair is missing one or both UMI3s...
                            dist3 = 15 #dist3 is 15 because UMI3 is 15bp long
                            status = 'single' #status is changed from 'both' to 'single'

                        match = 0 #set match variable to 0
                        if status == 'both': #if the read pair passed the length and UMI requirements...
                            if dist5 + dist3 <= 2: #ifthe cumulative Levenshtein distance between the two UMIs is 2 or less...
                                match = 1 #match variable = 1
                        if match == 1: #if match is 1
                            UMI_dict[group[j][0]] = UMI_dict[group[j][0]]|UMI_dict[group[i][0]] #this creates a list {1,2} that has the matching values
                            UMI_dict[group[i][0]] = UMI_dict[group[j][0]]|UMI_dict[group[i][0]]

            for i in range(0, len(group), 1):
                for j in range(i+1, len(group), 1):
                    if np.abs(len(group[i][3])-len(group[j][3]))/len(group[i][3]) < 0.1:
                        status = 'both'
                        if len(group[i][1]) > 0 and len(group[j][1]) > 0:
                            dist5 = editdistance.eval(group[i][1], group[j][1])
                        else:
                            dist5 = 15
                            status = 'single'
                        if len(group[i][2]) > 0 and len(group[j][2]) > 0:
                            dist3 = editdistance.eval(group[i][2], group[j][2])
                        else:
                            dist3 = 15
                            status = 'single'

                        match = 0
                        if status == 'both':
                            if dist5 + dist3 <= 2:
                                match = 1
                        if match == 1:
                            UMI_dict[group[j][0]] = UMI_dict[group[j][0]]|UMI_dict[group[i][0]]
                            UMI_dict[group[i][0]] = UMI_dict[group[j][0]]|UMI_dict[group[i][0]]
     #UMI_dict is now an ordered dictionary of each read in the group, with matching reads having the same dictionary values, i.e. {2,3}

            for entry in UMI_dict: #now we have a UMI dictionary where the reads are still grouped but read pairs that have passed the L distance test have their set populated with their read match
                counter_set = UMI_dict[entry] #counter_set is the UMI_dict value for that entry (so the tuple values); i.e {1,2}, {1,2}, and {3}
                if not set_dict.get(tuple(counter_set)): #set_dict is an empty dictionary that is re-set for each new group
                    UMI_group += 1
                    set_dict[tuple(counter_set)] = UMI_group

            read_list = [] 
            for i in range(0, len(group), 1): 
                UMI_number = set_dict[tuple(UMI_dict[group[i][0]])] 
                read_list.append(('>%s|%s\n%s\n' % (group[i][0], str(UMI_number), group[i][3]), UMI_number, group[i][1], group[i][2]))
            
            previous_UMI = ''
            Molecule = set() #defines Molecule as an empty set 
            for read, UMI_number, umi5, umi3 in sorted(read_list, key=lambda x:int(x[1])):
                matched_reads.write(str(UMI_number) + '\t' + read.split('|')[0] + '\t' + umi5 + '\t' + umi3 + '\n')
                if UMI_number != previous_UMI: #previous_UMI is empty, so UMI_number 
                    if len(Molecule) == 1:
                        final.write(list(Molecule)[0])
                    elif len(Molecule) > 1:
                        new_read = make_consensus(list(Molecule), previous_UMI, subreads) #make_consensus function defined earlier has the input arguments
                                                                                          #(Molecule, UMI_number, subreads)
                        if not new_read:
                            continue
                        final.write(new_read)
                        final_UMI_only.write(new_read)
                    Molecule = set()#if the current_UMI does not match previous UMI and Molecule is empty
                    Molecule.add(read) #if the current_UMI does not match 'previous_UMI' then add the read name to 'Molecule'
                    previous_UMI = UMI_number #if the previous_UMI is empty, then 'previous_UMI' becomes the current UMI_number 

                elif UMI_number == previous_UMI:
                    Molecule.add(read)

            if len(Molecule) == 1:
                final.write(list(Molecule)[0])
            elif len(Molecule) > 1:
                new_read = make_consensus(list(Molecule), previous_UMI, subreads)
                if not new_read:
                    continue
                # print('new_read', new_read)
                # print('written')
                final.write(new_read)
                # print('wrote')
                final_UMI_only.write(new_read)
                # print('wrote')
        elif len(group) > 0:
            UMI_group += 1
            final.write('>%s|%s\n%s\n' % (group[0][0], str(UMI_group), group[0][3]))
    # print(group_counter)

def read_UMIs(UMI_file): #defines function read_UMIs, that is used in 'main' script with the umi_file input 
    UMI_dict, group_dict, kmer_dict = {}, {}, {} #setting variables UMI_dict, group_dict, kmer_dict as empty dictionaries
    group_number = 0 #sets variable group_number as 0
    for line in open(UMI_file): #for each line in the UMI file, again organized as 'name'\t'UMI5'\t'UMI3'...
        a = line.strip('\n').split('\t') #remove spaces at beginning and end of each line and split line by tabs
        name, UMI5, UMI3 = a[0], a[1], a[2] #name and UMI5 and UMI3 are now defined correctly by their position on each line
        kmer_list = ['', '', '', ''] #defining variable kmer_list as an empty list with four indix positions

        UMI_dict[name] = (UMI5, UMI3) #populate the UMI dictionary; the name of the read is the key and the UMI5 and the UMI3 values are stored under it
        if UMI5[5:10] == 'TATAT': #if the middle of UMI5 is TATAT...
            kmer_list[0] = UMI5[:5] #kmer_list at position 0 is the UMI5 sequence from 0 to 5 (NNNNN) first random pentamer barcode sequence
            kmer_list[1] = UMI5[10:] #kmer_list at position 1 is the UMI5 sequence from 10 onward (NNNNN) second random pentamer barcode sequence
        if UMI3[5:10] == 'ATATA': #if the middle of UMI3 is ATATA (which we have confirmed will NOT happen with a UMI5 that is TATAT, because they are both on the same strand...)
            kmer_list[2] = UMI3[:5] #kmer_list at position 2 is the first random pentamer barcode sequence of UMI3
            kmer_list[3] = UMI3[10:] #kmer_list at position 3 is the second random pentamer barcode sequence of UMI3

        combination_list = []
        for x in range(4):
            for y in range(x+1, 4):
                first_kmer = kmer_list[x]
                second_kmer = kmer_list[y]
                if first_kmer != '' and second_kmer != '':
                    combination = '_' * x+first_kmer + '_' * (y-x)+second_kmer + '_' * (3-y)
                    combination_list.append(combination) #combination list is a coded list that represents the relationship between the four kmers in kmer list 
        match = ''
        for combination in combination_list: #for each of the combinations in the combination list...
            if combination in kmer_dict: #if the combination is in the kmer_dict...kmer_dict is an empty dictionary so this will not be true for the first combination
                match = kmer_dict[combination]
                for combination1 in combination_list: 
                    kmer_dict[combination1] = match
                break
        if match == '': #match will be empty for the first pass of the 'for' loop
            group_number += 1 #group number for the first pass will be 1 (it is defined outside of the loop as 0
            group_dict[group_number] = [] #the group dictionary at index 1 will be empty
            match = group_number #match is equal to 1
            for combination in combination_list: #for each combination in combination_list
                kmer_dict[combination] = match #kmer_dict for each of these key combinations will be 1
        group_dict[match].append(name) #group dictionary at position 'match' will also store the name of the read
    return group_dict, UMI_dict
    #this function creates two outputs 1) group_dict which is a dictionary under which reads are grouped together the key is the group_number and the values are the grouped reads, if we think of the four UMI kmers as UMI5.1,5.2, 3.1, 3.2 then 
    #reads have to have at least two matching kmers (same sequence and same positions--which is virtually impossible for any random pair of reads) and 2) UMI_dict which stores UMI5 UMI3
    #under the read name
def processing(reads, sub_reads, UMIs, groups, final, final_UMI_only, matched_reads): #defines the processing function with the given arguments
    annotated_groups, chrom_reads = parse_reads(reads, sub_reads, UMIs) #annotated_groups, chrom_reads are the products of parse_reads which is the groups list, and the rootname:[] dictionary
    subreads = read_subreads(subreads_file, chrom_reads) #subreads is the product of read_subreads(subreads_file, chrom_reads) function; it is just a dictionary where root read names are the key and all the fastq subreads and info is stored under that key
                                                         #but no group info is stored here!!!
    # print('grouping and merging consensus reads')
    group_reads(annotated_groups, reads, subreads, UMIs, final, final_UMI_only, matched_reads) #call to group_reads function--how does it work?

def main():
    final = open(path + '/R2C2_full_length_consensus_reads_UMI_merged.fasta', 'w') #final is this R2C2_full_length_consensus_reads_UMI_merged.fasta file
    final_UMI_only = open(path + '/R2C2_full_length_consensus_reads_UMI_only.fasta', 'w') #final_UMI_only is this R2C2_full_length_consensus_reads_UMI_only.fasta file
    matched_reads = open(path + '/matched_reads.txt', 'w') #matched_reads is this matched_reads.txt file
    print('kmer-matching UMIs')
    groups, UMIs = read_UMIs(umi_file) #'groups' variable and 'UMIs' variable are defined at the two outputs from the read_UMIs file: 1) group_dict and 2)UMI_dict 
    print('reading consensus reads')
    reads = read_fasta(fasta_file) #'reads' is the result of the read_fasta function performed on the fasta_file (the R2C2_Consensus.fasta input file), read_fasta produces a dictionary where sequence is stored under read name key
    count = 0 #set variable count=0
    sub_reads = {} #set variable 'sub_reads' as empty tuple
    for group in tqdm(sorted(groups)): #tqdm creates a progress bar, sorted(groups) returns the group numbers
        count += len(groups[group]) #this counts all the reads, essentially right? Mine are around 300k
        for name in groups[group]: #for each name in group 1 for example
            sub_reads[name] = group #sub_reads at position 'name' contains the group number. This basically creates a dictionary of reads names that store which matching group they are in
        if count > 500000: #my count is 300k---I don't know what this loop does, does it break up processing?
            # print('processing')
            processing(reads, sub_reads, UMIs, groups, final, final_UMI_only, matched_reads)
            count = 0
            sub_reads = {}
    processing(reads, sub_reads, UMIs, groups, final, final_UMI_only, matched_reads) #this calls the processing function with all of the listed arguments--go to processing function to see what it does
    #reads is the read:sequence dictionary from the input R2C2_Consensus.fasta
    #sub_reads is the readname:group dictionary created in this function
    #UMIs is the readname:UMI5,UMI3 dictionary generated from the read_UMIs(umi_file) function
    #groups is the group:readname dictionary created by read_UMIs(umi_file)
    #final, final_UMI_only and matched_reads are the files created in this function

if __name__ == '__main__':
    main()
