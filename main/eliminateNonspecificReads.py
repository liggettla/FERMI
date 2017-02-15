#!/usr/bin/env python
# The purpose of this script is to eliminate all variants that did not fall within
# the specifically probed regions. It is strange that these can occur with such
# large frequencies, but ignoring them may be a good choice until a better
# solution is found.

def runArgParse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--inFile', '-i', type=str, required=True, help='The input vcf to be analyzed.')
    parser.add_argument('--outFile', '-o', type=str, required=True, help='The output file containing only variants that fall within probed regions.')

    args = parser.parse_args()

    inFile = args.inFile
    outFile = args.outFile

    return inFile, outFile

def elimBadAligns(inFile, outFile):
    from parseLine import seqRead

    inTarget = open(inFile, 'r')
    outTarget = open(outFile, 'w')

    targetLocs = {'chr1':['1152278','1152279','1152564','1152565','1152566','1152587','1152588'],'chr2':['254572','254573','2091130','2091131','2091132','1982668','1982669','2231906','2231907','2231908','2290411','2290412'],'chr4':['1061972','1061973','1061974','1061551','1061552','1105411','1105412','1105413','1129972','1129973','1211677','1211678','1235477','1235478','1235479','1244286','1244287'],'chr9':['50737','50738'],'chr11':['21262','21263','21264','23899','2390','2593','2594','114865','114866','114867','5342','5343'],'chr12':['253982','253983','253984','253802','253803'],'chr15':['925270','925271','906318','906319'],'chr16':['733796','733797','733798','824550','824551','859491','859492'],'chr17':['75775','75776','75783','75784','75785','75770','75771','75772'],'chrX':['486496','486497','486498']}

    for line in inTarget:
        if '#' not in line and 'chr' in line: # skip the damn info
            lineObj = seqRead(line)
            if lineObj.chrom() in targetLocs: # is this chrom probed?
                goodLoc = False
                for i in targetLocs[lineObj.chrom()]:
                    if i in lineObj.loc(): # is location probed?
                        goodLoc = True
                if goodLoc == True:
                    outTarget.write(line)
        elif '#' in line:
            outTarget.write(line)

    inTarget.close()
    outTarget.close()

if __name__ == "__main__":
    inFile, outFile = runArgParse()
    elimBadAligns(inFile, outFile)
