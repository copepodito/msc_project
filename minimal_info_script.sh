#!/usr/bin/env bash

# This bash script creates a file named "minimal_info.csv" for each subject.
# Each subject folder must contain the consensus files (18) from https://github.com/briney/grp_paper.

for d in ~/Desktop/Project/*;
do
  if [ -d "$d" ]; then
    cd $d
    # Create a .csv file and add headers
    echo "v_gene,j_gene,cdr3_length,cdr3_aa,var_mut_count_aa,isotype" > minimal_info.csv

    # Loop trough all files within a directory and populate the file
    for f in ./*.txt
    do
      # Avoid header (first line), avoid non-productive abs and light chains, and append only selected columns
      tail -n +2 $f | grep "yes" | grep "heavy" | cut -f6,10,11,13,22,25 -d "," >> minimal_info.csv
    done
  fi
done

