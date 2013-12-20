#!/usr/bin/python
DIRECTORY='/home/carles/QCB301/promoters'
PROMOTERS='/home/carles/QCB301/ReferenceGenomesandAnnotations/promoter_sequences.fasta'
MOTIFS='/home/carles/QCB301/ReferenceGenomesandAnnotations/1.02/ALIGNED_ENOLOGO_FORMAT_PWMS/'
# Some default list, to be given.
TARGETS='/home/carles/QCB301/Lists/AllGenes'

import os
import sys

lsdir = os.listdir(MOTIFS)
pwm = []
for file in lsdir:
    if (file != 'AllMotifs.pwm'):
        pwm.append(file.split("_")[0])

pwm.sort()
pwm = set(pwm)
zeros = [0] *len(pwm)

d = dict(zip(pwm,zeros))

lsdir = os.listdir(DIRECTORY)

for file in lsdir:
    if (file.split(".")[1] == 'pcounts'):
        with open(DIRECTORY + "/" + file,'r') as readin:
            header=readin.readline()
            for line in readin:
                a = line.split("\n")[0].split(" ")
                d[a[-1]] = d[a[-1]] + int(a[-2])

with open('QCB301/AllMotifs.pcounts','w') as out:
    for key in d:
        out.write(key + " " + str(d[key]) + "\n")
    
            
