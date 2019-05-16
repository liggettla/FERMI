#!/usr/bin/env/ python

def identifyVars(indiv, samples=False, cutoff=False, confLevel=0.9999999999, indir='sampleData', ref=False, suffix=False):
    from util import defineProbes, populatePandasDataframe, assessBackground
    from util import findSigVars

    # set directories
    if not ref:
        ref = '../ReferenceGenomes/hg19.fa'

    # specify default samples to use as background
    if not samples:
        samples = ['15','7','80', '81','82','83','84','85']

    # exclude indiv of interest from background calculation
    if indiv in samples:
        samples.remove(indiv)
    indiv = [indiv]

    # compute background variant rates
    rates, allRates, bgRange = assessBackground(indir, samples, confLevel)

    # set probe regions
    probes = defineProbes()

    # build some sample lists to be passed to
    # populatePandasDataframe()
    if not suffix:
        suffix = '.fastq/onlyProbedRegions.vcf'
    fullSamples = []
    #probes = defineProbes()

    for i in indiv:
        fullSamples.append('%s%s' % (i, suffix))

    # number of flanking bases
    # up, down = 1, 1
    up, down = False, False

    # populate DataFrame with individual of focus
    CR = populatePandasDataframe(indir, fullSamples, probes, ref, up, down)
    if cutoff:
        CR = CR[CR.AO > cutoff]

    # calculate significantly elevated variants
    df = findSigVars(rates, CR)
    #print('Number of vars: %s' %(len(df)))
    return df, bgRange, allRates

# this is a method to understand the background rates of variants
def assessBackground(indir, samples, conf = 0.95):
    import math
    from scipy.stats import t
    import numpy as np
    from collections import defaultdict

    #print('Computing background variant rates...')

    allRates = defaultdict(list)
    rates = {}
    bgRange = {}

    # compile dict of all vafs
    for vcf in samples:
        vcf = open(('%s/%s.fastq/onlyProbedRegions.vcf' % (indir, vcf)), 'r')
        for line in vcf:
            if '#' not in line and 'chr' in line: # skip the info
                lineobj = VCFObj(line)

                # exclude germline and clonally expanded mutations from
                # the distribution calculations
                if lineobj.af < 0.3:
                    allRates[lineobj.index].append(lineobj.af)

    #print('Computing confidence intervals...')

    for var in allRates:
        confidence = (
        t.interval(conf,len(allRates[var])-1,loc=np.mean(allRates[var]),scale=np.std(allRates[var])/math.sqrt(len(allRates[var])))[1])
        rates[var] = [confidence]

        backgroundRange = (
        t.interval(conf,len(allRates[var])-1,loc=np.mean(allRates[var]),scale=np.std(allRates[var])/math.sqrt(len(allRates[var]))))
        bgRange[var] = [backgroundRange]

    return rates, allRates, bgRange

def populatePandasDataframe(dirinput, fileList, probes, ref, upstream=10, downstream=10, ignore_indels=True):
    import pandas as pd
    from Bio.Seq import Seq
    #print('Building data structure...')

    allSamples = []
    columns = ['Chrom','Loc','WT','Var','Change','ConvChange','AO','DP','VAF','IntEx','Gene','Upstream','Downstream','Individual', 'ID']
    dat = []

    tempAllVariants = []
    sampleCount = 0
    for sample in fileList:
        inFile = open(dirinput + '/' + sample, 'r')
        sampleCount += 1

        for line in inFile:
            if '#' not in line and 'chr' in line: # skip the info
                lineobj = VCFObj(line)
                # convert to six changes
                if lineobj.wt == 'G' or lineobj.wt == 'A':
                    wt = str(Seq(lineobj.wt).complement())
                    var = str(Seq(lineobj.var).complement())
                else:
                    wt = str(lineobj.wt)
                    var = str(lineobj.var)

                if upstream and downstream:
                    surrounding = getRefSequence(lineobj, upstream, downstream, ref)
                    up = str(surrounding[:upstream])
                    down = str(surrounding[-downstream:])
                else:
                    up = 'NaN'
                    down = 'NaN'

                probeRegion = ''
                for probe in probes:
                    if len(probeRegion) < 1:
                        for loc in probes[probe]:
                            if str(loc) in str(lineobj.location):
                                if probe[0] == 'T':
                                    probeRegion = 'TIII'
                                    specificProbe = probe
                                else:
                                    probeRegion = 'Exon'
                                    specificProbe = probe

                # catch other probes not in list that might
                # be part of a new probeset
                if probeRegion == '':
                    probeRegion = 'Other'
                    specificProbe = 'Other'


                if ignore_indels:
                    if len(lineobj.wt) == 1 and len(lineobj.var) == 1:
                        dat = [lineobj.chrom, lineobj.location, str(lineobj.wt), str(lineobj.var), str(lineobj.wt) + '>' + str(lineobj.var), wt + '>' + var, lineobj.ao, lineobj.dp, lineobj.af, probeRegion, specificProbe, up, down, sampleCount, sample]
                        tempdat = pd.DataFrame(dat, index=columns)
                        tempAllVariants.append(tempdat.T)
                else:
                    dat = [lineobj.chrom, lineobj.location, str(lineobj.wt), str(lineobj.var), str(lineobj.wt) + '>' + str(lineobj.var), wt + '>' + var, lineobj.ao, lineobj.dp, lineobj.af, probeRegion, specificProbe, up, down, sampleCount, sample]
                    tempdat = pd.DataFrame(dat, index=columns)
                    tempAllVariants.append(tempdat.T)


        inFile.close()
    allVariants = pd.concat(tempAllVariants, ignore_index=True)

    return allVariants

