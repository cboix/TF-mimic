#!/usr/bin/python
# Scan for motifs that correspond to tfs expressed above a certain level.
REFRESH=False#canonically false, but can be overwritten. THIS whole section should be parsed as args.
DIRECTORYtwo='/home/carles/QCB301/threeprom'
DIRECTORY='/home/carles/QCB301/promoters'
PROMOTERS='/home/carles/QCB301/ReferenceGenomesandAnnotations/promoter_sequences.fasta'
# Some default list, to be given.
TARGETS='/home/carles/QCB301/Lists/AllGenes'
THRESHOLD=0.025

## CONDITIONS FOR ADDING THE CUT THREE TFS! comment if not needed...
REFRESH=True
MOTIFS='/home/carles/QCB301/MotifsT/PWMs/AllThree.pwm'

import os.path
import sys
from math import exp,log

# Must run script on its own / w/out "python"
if (len(sys.argv) > 1):
    DIRECTORY = sys.argv[1]
    TARGETS = sys.argv[2]

def read_motif_table(lines):
    first=lines[0].split("\n")[0]
    header=first.split("_")
    gene=header[0].split(">")
    num=header[-1]
    d = dict(gene=gene[1],number=num,iden=first)
    for line in lines[1:5]:
        proc=line.split("\n")[0]
        proc=proc.split("\t")
        name=proc[0]
        info=proc[1:]
        d.update(zip(name,[info]))
    d['length'] = len(d['T'])
    d["consensus"] = 0
    for i in range(0,len(d['T'])):
       d["consensus"] = d["consensus"] + max(float(d['A'][i]),float(d['G'][i]),float(d['C'][i]),float(d['T'][i]))
    return(d)

# TODO if I want to do my own transformation from count data:
  #  for i in range(0,d['length']):
  #      # TODO check order.
  #      total = float(d['A'][i]) + float(d['G'][i]) + float(d['T'][i]) + float(d['C'][i]) + 4
  #      d['A'][i] = log((float(d['A'][i])+1)/total/.25)
  #      d['C'][i] = log((float(d['C'][i])+1)/total/.25)
  #      d['G'][i] = log((float(d['G'][i])+1)/total/.25)
  #      d['T'][i] = log((float(d['T'][i])+1)/total/.25)

#CREATION of the motif dictionaries: (NOTE to look at Labwork for restricting motifs)
def read_PWMs(filename,threshold):
    with open(filename,'r') as infile:
        PWM=[]
        lines = []
        # Append lines until we have 5 (one motif)
        for line in infile:
            lines.append(line)
            if len(lines) >= 5:
                dic = read_motif_table(lines)
                PWM.append(dic)
                lines = []
        if len(lines) > 0:
            print('something left!')
        return(PWM)

# Could make this more efficient by cutting off a threshold earlier
# To make this more efficient would be to have a matrix (length *4) which you can multiply against the pwm
def compare_motif(lines,pwm,consensus,threshold,strand):
    value=0
    breaker=0
    i=0
    header=lines[0]
    last = header[1]
    for line in lines:
        this = line[1]
        #if not continuous, break!
        if abs(this - last) > 1:
            breaker=1
            break
        #Compare the dmelanogaster(first one) or yeast sequence only.
        value = value + float(pwm[line[0].upper()][i])
        # This can and should be adjusted.
        if value < -30:
            breaker=1
            break
        last = this
        i=i+1
    # compute the occupancy and see if it fullfills the threshold:
    # DUE TO OVERFLOW:
    if consensus - value < 100:
        o = 1/(1+exp(consensus - value))
    else:
        o = 0
    if (o >= threshold):
        d = dict(start=header[1],end=last,v=value,occ=o,pwm=pwm['gene'],iden=pwm['iden'],strand=strand,consensus=consensus)
        # The old definition of dictionary:
        # d = dict(start=header[2],end=end,position=header[1],peaksep=header[3],v=value,occ=o,pwm=pwm['gene'],iden=pwm['iden'])
        # All of these 0s will be cut out later (in find_motifs):
        if breaker==0: 
            return(d)
        else:
            return(0)
    else:
        return(0)

def revline(line):
    if(line[0].upper() == 'A'):
        return(['T',line[1]])
    if(line[0].upper() == 'G'):
        return(['C',line[1]])
    if(line[0].upper() == 'C'):
        return(['G',line[1]])
    if(line[0].upper() == 'T'):
        return(['A',line[1]])
    else:
        return(line)

def find_motifs(filename,pwm,threshold):
    motif = []
    consensus = pwm['consensus']
    with open(filename,'r') as promoter:
        place = -1000
        lines = []
        backlines = []
        for line in promoter:
            # define position, make it a simple array
            line = [line[0].split("\n")[0] , place]
            lines.append(line)
            backlines.insert(0,revline(line))
            if len(lines) == pwm['length']:
                motif.append(compare_motif(lines,pwm,consensus,threshold,"+"))
                # Here get rid of the earliest line:
                lines = lines[1:]
                place = place + 1
        # REPEAT exact same process but backwards:
        for line in backlines:
            lines.append(line)
            if len(lines) == pwm['length']:
                motif.append(compare_motif(lines,pwm,consensus,threshold,"-"))
                # Here get rid of the earliest line:
                lines = lines[1:]
    matches=[]
    for m in motif:
        if m != 0:
            matches.append(m)
    return(matches)

def compute_file(filename,PWM,threshold):
    motifs = []
    for pwm in PWM:
        out = find_motifs(filename,pwm,threshold)
        motifs.append(out)
    d = []
    for m in motifs:
        if len(m) > 0:
            for dic in m:
                d.append(dic)
    #processing to be done here to select motif threshold. 
    return(d)

def gene_compute(gene,directory,refresh):
    namein = directory + "/" + gene + '.promoter'
    nameout = DIRECTORYtwo + "/" + gene + '.tfbs'
    headerline = 'pwm gene start end occ iden consensus value strand\n'
    if (os.path.exists(namein)):
        # If out file exists, don't do this unless refresh.
        if (not os.path.exists(nameout) or refresh):
            motifs=compute_file(namein,PWM,THRESHOLD)
            with open(nameout,'a') as out:
                out.write(headerline);
                for m in motifs:
                    outline = m['pwm'] + " " + gene + " " + str(m['start']) + " " + str(m['end']) + " " + str(m['occ']) + " " + m['iden'] + " " + " " + str(m['consensus']) + " " + str(m['v']) + " " + m['strand']  + "\n"
                    out.write(outline)

PWM = read_PWMs(MOTIFS,3)

# TODO Loop over all gene targets given and compute the enhancer searches for all:
with open(TARGETS,'r') as tar:
    for geneline in tar:
        gene = geneline.split("\n")[0]
        gene_compute(gene,DIRECTORY,REFRESH)

# NOTE 
# AFTER THIS, run the conversion from tfbs to pcounts in the 'createlists' file

#To test:
# pwm=PWM[1]
# motif=find_motifs('/home/carles/Labwork/StarrSeq/cbt.yak.div',pwm,0.25)
#start = time.time()
#motifs=compute_file('/home/carles/Labwork/StarrSeq/cbt.yak.div',PWM,0.25)
#elapsed = time.time() - start
