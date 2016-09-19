#This calculates the average coverage of each individual capture
def outputCov(input_file, final_output_file, distance_stringency, coverage_file, averageErrorRate):
    target = open(coverage_file, 'w')

    inputLines = sum(1 for line in open(input_file))
    outputLines = sum(1 for line in open(final_output_file))
    target.write("Eliminating only exact and close UMI matches:\n")
    target.write("# of Allowed UMI Mismatches: %d\n" %(distance_stringency))
    if not inputLines == 0:
        target.write("Total # of Original UMIs: %d\n" % (inputLines/4))
        target.write("# of Unique UMIs: %d\n" % (outputLines/4))
        avgCov = float(inputLines)/outputLines
        target.write("Avg UMI Coverage: %r\n" % (avgCov))
    target.write("Average Error Rate: %f\n" % (averageErrorRate))

    target.close()
