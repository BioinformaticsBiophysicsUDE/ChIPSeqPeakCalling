# ChIPSeqPeakCalling

## Peak calling with hiddenDomains3.0

http://hiddendomains.sourceforge.net/

## sliding window domain merging

Merging of called peaks in close proximity using a sliding window approach. The input is the *.analysis.bed file output from hiddenDomains. The function uses three paramter bins, windows and treshold which have different values for broad an narrow histone modifications.

bins - bin size for the binned genome

windows - window size of the sliding windows

treshold -  if more then treshold positions in a window are called as peak, the whole window is set to peak

for narrow histone modifications the settings are: bins = 200, window = 1000 and treshold = 800

for broad histone modifications the settings are: bins = 800, window = 5600, treshold = 4000

call the function in R 

domainMergerSeq(<analysis.bed file>, <broad or narrow>)
