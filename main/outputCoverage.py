#This calculates the average coverage of each individual capture
from numpy import asscalar

def outputCov(input_file, final_output_file, distance_stringency, coverage_file, averageErrorRate, averageCoverage):
    target = open(coverage_file, 'w')

    inputLines = sum(1 for line in open(input_file))
    outputLines = sum(1 for line in open(final_output_file))
    target.write("Eliminating only exact and close UMI matches:\n")
    target.write("# of Allowed UMI Mismatches: %d\n" %(distance_stringency))
    if not inputLines == 0:
        target.write("Total # of Original UMIs: %d\n" % (inputLines/4))
        target.write("# of Unique UMIs: %d\n" % (outputLines/4))
        avgCov = float(inputLines)/outputLines
        target.write("Avg UMI Coverage (Incl Unused Reads): %r\n" % (avgCov))

        averageCoverage = asscalar(averageCoverage)
        target.write("Avg UMI Coverage (Only Used Reads): %f\n" % (averageCoverage))

    averageErrorRate = asscalar(averageErrorRate)
    target.write("Average Error Rate: %f\n" % (averageErrorRate))

    target.close()
