import sys
import pandas as pd
import numpy as np

#print sys.argv[1]

infile = sys.argv[1]
bedfile = pd.read_csv(filepath_or_buffer=infile, sep="\t",
                      header=None,
                      names=["chrom", "start", "end", "strand", "peak"],
                      dtype={'start': np.int64})

#print(bedfile.head())
for index, row in bedfile.iterrows():
    peaks = row['peak'].split(',')
    peaks = np.array(peaks)
    peaks = peaks.astype(np.float)
    peak_max = peaks.max()

    if len(peaks) > 1:
        start = np.array(range(row['start'], row['start']+len(peaks), 1))
        dompeak = start[peaks == peak_max]
        #print(list(row))
        print(row['chrom'], dompeak[0], dompeak[0]+1,
              'dominantPeak(%d)' % len(peaks), peak_max, row['strand'],
              sep="\t")
    else:
        print(row['chrom'], row['start'], row['end'],
              'dominantPeak', peak_max, row['strand'],
              sep="\t")
