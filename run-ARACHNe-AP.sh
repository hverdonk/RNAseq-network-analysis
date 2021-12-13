#!/bin/bash
# more information can be found at https://github.com/califano-lab/ARACNe-AP

output_dir='your/path/here'
transcript_counts='your/path/here/'mini-expression-data.tsv   # enter the full path to the count matrix
transcription_factors='your/path/here/'TFs.txt           # enter the full path to the transcription factor list
ARACHNe_dir='your/path/here/'ARACNe-AP    # enter the path to the ARACHNe program root directory

# move to ARACHNe root directory to run the analysis
curr_dir=$PWD
cd $ARACHNe_dir

# later, replace '/usr/bin/javac' with 'java' to avoid user-specific pathing issues

# Calculate a threshold for Mutual Information
java -Xmx5G -jar dist/aracne.jar -e $transcript_counts  -o $output_dir --tfs $transcription_factors --pvalue 1E-8 --seed 1 --calculateThreshold

# Run ARACNe on 10 reproducible bootstraps of the input matrix (your expression count data)
# For my true data, I ran ARACHNe-AP on 300 reproducible bootstraps
for i in {1..10}
do
java -Xmx5G -jar dist/aracne.jar -e $transcript_counts  -o $output_dir --tfs $transcription_factors --pvalue 1E-8 --seed $i --threads 16
done

# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar dist/aracne.jar -o $output_dir --consolidate

# return to the original directory
cd $curr_dir

