# ChIPSeqPeakCalling

## Peak calling with hiddenDomains3.0

http://hiddendomains.sourceforge.net/

## Sliding window domain merging

domain.merger.function.R

Merging of called peaks in close proximity using a sliding window approach. The input is the *.analysis.bed file output from hiddenDomains. The function uses three paramter bins, windows and treshold which have different values for broad an narrow histone modifications.

bins - bin size for the binned genome

windows - window size of the sliding windows

treshold -  if more then treshold positions in a window are called as peak, the whole window is set to peak

for narrow histone modifications the settings are: bins = 200, window = 1000 and treshold = 800

for broad histone modifications the settings are: bins = 800, window = 5600, treshold = 4000

source file in R and call the function with

domainMergerSeq(<analysis.bed file>, <broad or narrow>)


## Calulate intersect peak file from biological replicates

getPeakIntersect.R

input are the window merged peak bed files

input is a vector with 2 or 3 bed files and the name of the output file


source file in R and call the function

getIntersectVectNoID(<vector with bed files>, < save name>)

