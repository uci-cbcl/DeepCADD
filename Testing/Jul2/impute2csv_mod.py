#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *18.06.2014
"""

import sys, os
import random
from optparse import OptionParser
import math
import gzip

random.seed(1234567890)

transform = [("y",None),
             ("Chrom","REMOVE"),
             ("Pos","REMOVE"),
             ("Ref","OK"),
             ("Alt","OK"),
             ("Type","OK"),
             ("Length","OK"),
             ("isTv","OK"),
             ("Consequence","OK"),
             ("GC","OK"),
             ("CpG","OK"),
             ("priPhCons","OK"),
             ("mamPhCons","OK"),
             ("verPhCons","OK"),
             ("priPhyloP","OK"),
             ("mamPhyloP","OK"),
             ("verPhyloP","OK"),
             ("GerpN","OK"),
             ("GerpS","OK"),
             ("GerpRS","OK"),
             ("GerpRSpval","OK"),
             ("bStatistic","OK"),
             ("EncExp","OK"),
             ("EncH3K27Ac","OK"),
             ("EncH3K4Me1","OK"),
             ("EncH3K4Me3","OK"),
             ("EncNucleo","OK"),
             ("EncOCC","OK"),
             ("EncOCCombPVal","OK"),
             ("EncOCDNasePVal","OK"),
             ("EncOCFairePVal","OK"),
             ("EncOCpolIIPVal","OK"),
             ("EncOCctcfPVal","OK"),
             ("EncOCmycPVal","OK"),
             ("EncOCDNaseSig","OK"),
             ("EncOCFaireSig","OK"),
             ("EncOCpolIISig","OK"),
             ("EncOCctcfSig","OK"),
             ("EncOCmycSig","OK"),
             ("Segway","OK"),
             ("tOverlapMotifs","OK"),
             ("motifDist","OK"),
             ("motifECount","OK"),
             ("motifEHIPos","OK"),
             ("motifEScoreChng","OK"),
             ("TFBS","OK"),
             ("TFBSPeaks","OK"),
             ("TFBSPeaksMax","OK"),
             ("DAF","REMOVE"),
             ("minDistTSS","OK"),
             ("minDistTSE","OK"),
             ("cDNApos","IND"),
             ("relcDNApos","IND"),
             ("CDSpos","IND"),
             ("relCDSpos","IND"),
             ("protpos","IND"),
             ("relProtpos","IND"),
             ("Dst2Splice","OK"),
             ("Dst2SplType","OK"),
             ("oAA","OK"),
             ("nAA","OK"),
             ("Grantham","IND"),
             ("PolyPhenCat","OK"),
             ("PolyPhenVal","IND"),
             ("SIFTcat","OK"),
             ("SIFTval","IND"),
             ("IND_cDNApos",None),
             ("IND_relcDNApos",None),
             ("IND_CDSpos",None),
             ("IND_relCDSpos",None),
             ("IND_protpos",None),
             ("IND_relProtpos",None),
             ("IND_Dst2Splice","REMOVE"),
             ("IND_Grantham",None),
             ("IND_PolyPhenVal",None),
             ("IND_SIFTval",None)]

categorical = { "Type" : ["SNV","INS","DEL"],
                "Consequence" : ["U3","U5","DN","IG","I","NC","IF","FS","NS","R","CS","S","SG","SL","SN","O","UP"],
                #["3PRIME_UTR","5PRIME_UTR","DOWNSTREAM","INTERGENIC","INTRONIC","NONCODING_CHANGE","INFRAME","FRAME_SHIFT","NON_SYNONYMOUS","REGULATORY","CANONICAL_SPLICE","SPLICE_SITE","STOP_GAINED","STOP_LOST","SYNONYMOUS","UNKNOWN","UPSTREAM"],
                "Segway":["C0","C1","D","E/GM","F0","F1","GE0","GE1","GE2","GM0","GM1","GS","H3K9me1","L0","L1","R0","R1","R2","R3","R4","R5","TF0","TF1","TF2","TSS","UD"],
                "Ref" : ["A","C","G","T","N"],
                "Alt" : ["A","C","G","T","N"],
                "Dst2SplType" : ["DONOR","ACCEPTOR","UD"],
                "PolyPhenCat" : ["benign","possibly_damaging","probably_damaging","unknown","UD"],
                "SIFTcat" : ["deleterious","tolerated","unknown","UD"],
                "oAA" : ["*","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","UD"],
                "nAA" : ["*","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","UD"] }

consequence_set = set(categorical["Consequence"])

parser = OptionParser("%prog [options]")
parser.add_option("--noGC",dest="noGC", help="Remove GC/CpG values from output lines",default=False,action="store_true")
parser.add_option("--noInteractions",dest="noInteractions", help="Remove additional interaction terms with Consequence",default=False,action="store_true")
parser.add_option("--noPolyPhen",dest="noPolyPhen", help="Remove PolyPhen values",default=False,action="store_true")
parser.add_option("-p","--portion", dest="portion", help="Portion of data to be used for training (def None)",default=None,type='float')
parser.add_option("--training", dest="training", help="Training output file (def train.csv.gz)",default="train.csv.gz")
parser.add_option("--testing", dest="testing", help="Testing output file (def test.csv.gz)",default="test.csv.gz")
parser.add_option("--noexp", dest="noexp", help="Replace exponential numbers",default=False,action="store_true")
parser.add_option("--assignment", dest="assignment", help="Print assignment of features",default=False,action="store_true")
(options, args) = parser.parse_args()

if options.noGC:
  for ind,(name,value) in enumerate(transform):
    if name == "GC" or name == "CpG":
      transform[ind] = (name,"REMOVE")
  sys.stderr.write('Removing GC/CpG...\n')

if options.noPolyPhen:
  for ind,(name,value) in enumerate(transform):
    if name.startswith("PolyPhen") or name.startswith("IND_PolyPhenVal"):
      transform[ind] = (name,"REMOVE")
  sys.stderr.write('Removing PolyPhen...\n')

catvals = {}
intvals = {}
indicators = {}
products = {}
AApairs = {}
BasePairs = {}
SplicePairs = {}

ind = 1
label_assignment = []
for key,operation in transform:
  if operation == None: continue
  elif operation == "REMOVE": continue
  else:
    if key in categorical:
      for category in categorical[key]:
        if key not in catvals: catvals[key] = [(category,ind)]
        else: catvals[key].append((category,ind))
        label_assignment.append((str(ind),(key,category)))
        ind += 1
    else:
      intvals[key] = ind
      label_assignment.append((str(ind),key))
      ind += 1
    if operation == "IND":
      indicators[key] = ind
      label_assignment.append((str(ind),"IND_"+key))
      ind += 1

# for key,value in intvals.iteritems():
#   if key.startswith("Enc"): continue
#   print key,value
# sys.exit()

prodlst = []
products = {}
if not options.noInteractions:
  prodlst = ["bStatistic","cDNApos","CDSpos","Dst2Splice","GerpN","GerpS","mamPhCons","mamPhyloP","minDistTSE","minDistTSS","priPhCons","priPhyloP","protpos","relcDNApos","relCDSpos","relProtpos","verPhCons","verPhyloP"]

  for key in categorical["Consequence"]:
    products[key] = {}
    for key2 in prodlst:
      products[key][key2] = ind
      label_assignment.append((str(ind),(key,key2)))
      ind += 1

for elem in categorical["Dst2SplType"][:-1]:
  SplicePairs[elem] = ind
  label_assignment.append((str(ind),(elem,"Dst2Splice")))
  ind += 1
  indicators[elem] = ind
  label_assignment.append((str(ind),"IND_"+elem))
  ind += 1

for oB in categorical["Ref"]:
  for nB in categorical["Alt"]:
    BasePairs[(oB,nB)] = ind
    label_assignment.append((str(ind),(oB,nB)))
    ind += 1

for oAA in categorical["oAA"][:-1]:
  for nAA in categorical["nAA"][:-1]:
    AApairs[(oAA,nAA)] = ind
    label_assignment.append((str(ind),(oAA,nAA)))
    ind += 1

if options.assignment:
  for key,value in label_assignment:
    print key,value
  sys.exit()

training = None
testing = None
if options.portion != None:
  training = gzip.open(options.training,'w')
  testing = gzip.open(options.testing,'w')

countlines = 0
transformNames = map(lambda (x,y):x,transform)
len_lassign = len(label_assignment)

# Print CSV header
if options.portion == None:
  sys.stdout.write("y,%s\n"%(",".join(map(lambda (key,value): "x".join(value) if type(value) == type(("A","B")) else value, label_assignment))))
else:
  training.write("y,%s\n"%(",".join(map(lambda (key,value): "x".join(value) if type(value) == type(("A","B")) else value, label_assignment))))
  testing.write("y,%s\n"%(",".join(map(lambda (key,value): "x".join(value) if type(value) == type(("A","B")) else value, label_assignment))))

for line in sys.stdin:
  countlines += 1
  fields = line.rstrip().split('\t')
  if line.startswith('#') or line.startswith('y'):
    if len(fields) != len(transform):
      sys.stderr.write('Unexpected number of fields in line %d: T:%d F:%d\n'%(countlines,len(transform),len(fields)))
      sys.stderr.write(str(transform)+"\n")
      sys.stderr.write(str(fields)+"\n")
      for ind in range(min(len(transform),len(fields))):
        sys.stderr.write("%s\t%s\n"%(transform[ind],fields[ind]))
      sys.exit()
  elif len(fields) == len(transform):
    fieldsDict = dict(zip(transformNames,fields))
    new_fields=[]
    cvalues = {}
    oAA,nAA = fieldsDict['oAA'],fieldsDict['nAA']
    Ref,Alt = fieldsDict['Ref'],fieldsDict['Alt']

    splice,spliceval = fieldsDict["Dst2SplType"],fieldsDict["Dst2Splice"]
    state = "-1" if fieldsDict['y'] == "0" else "1"

    for ind,(name,operation) in enumerate(transform):
      if (operation == "REMOVE"): continue
      if name == "y": continue

      found = False
      if name.startswith("IND_"):
        found = True
        if fields[ind] != "0":
          new_fields.append((indicators[name[4:]],"1"))
      else:
        cval = fields[ind]
        if name in intvals:
          found = True
          cvalues[name]=cval
          # SURPPRESS ZEROS FOR SPARSE MATRIX FORMAT
          if float(cval) != 0:
            new_fields.append((intvals[name],cval))
        elif name in catvals:
          cval = fields[ind]
          for cat,sind in catvals[name]:
            cvalues[(name,cat)]="0"
            if cat == cval:
              # ONLY REPORT ACTUAL CATEGORY FOR SPARSE MATRIX FORMAT
              new_fields.append((sind,"1"))
              found = True
              cvalues[cat]="1"

        if not found:
          sys.stderr.write("Unexpected field: %d %s %s %s!\n"%(ind,name,operation,fields[ind]))

    if fieldsDict["Consequence"] in products:
      for key2,pind in products[fieldsDict["Consequence"]].iteritems():
        hval2 = cvalues.get(key2,None)
        if hval2 != None and float(hval2) != 0: new_fields.append((pind,hval2))

    if splice != None and spliceval != "0":
      new_fields.append((SplicePairs[splice],spliceval))
    for splType in ["DONOR","ACCEPTOR"]:
      if splice != splType: new_fields.append((indicators[splType],"1"))

    new_fields.append((BasePairs[Ref,Alt],"1"))

    if len(oAA) == len(nAA) and len(oAA) == 1:
      new_fields.append((AApairs[oAA,nAA],"1"))

    # MAKE SURE LABELS ARE SORTED!
    new_fields.sort()
    if options.noexp:
      new_fields = map(lambda (x,y): (x,y) if ("e" not in y) else "%.8f"%(float(y)),new_fields)

    all_fields = []
    lind = 1
    for (key,value) in new_fields:
      while (lind < key):
        all_fields.append("0")
        lind+=1
      all_fields.append(value)
      lind+=1
    while (lind <= len_lassign):
      all_fields.append("0")
      lind+=1

    #sys.stderr.write("%d\t%d\t%d\n"%(len(all_fields),len_lassign,len(new_fields)))
    if options.portion == None:
      sys.stdout.write("%s,%s\n"%(state,",".join(all_fields)))
    else:
      if random.random() < options.portion: training.write("%s,%s\n"%(state,",".join(all_fields)))
      else: testing.write("%s,%s\n"%(state,",".join(all_fields)))

  else:
    sys.stderr.write('Unexpected number of fields in line %d: T:%d F:%d\n'%(countlines,len(transform),len(fields)))
    sys.stderr.write(str(transform)+"\n")
    sys.stderr.write(str(fields)+"\n")
    sys.exit()
