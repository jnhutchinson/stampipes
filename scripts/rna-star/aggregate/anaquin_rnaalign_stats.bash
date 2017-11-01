# !/bin/bash

inputfile=$1
outputfile=$2

# Dilution
echo -ne "sequins-dilution\t" > $2
cat $inputfile | grep 'Dilution' | sed -e 's/Dilution://g' | tr -d " \t\r" >> $2

# Alignments nonspliced
echo -ne "sequins-alignments-non-spliced\t" >> $2
cat $inputfile | grep -A 3 "Alignments (Synthetic)" | grep "Non-spliced" | sed -e 's/Non-spliced://g' | tr -d " \t\r" >> $2

# Alignments spliced
echo -ne "sequins-alignments-spliced\t" >> $2
cat $inputfile | grep -A 3 "Alignments (Synthetic)" | grep "Spliced" | sed -e 's/Spliced://g' | tr -d " \t\r" >> $2

# Intron level Sensitivity
echo -ne "sequins-intron-level-sensitivity\t" >> $2
cat $inputfile | grep -A 8 "Comparison of alignments to reference annotation" | grep -A 3 'Intron level' | grep 'Sensitivity' | sed -e 's/Sensitivity://g' | tr -d ' \t\r' >> $2
 
# Intron level Precision
echo -ne "sequins-intron-level-precision\t" >> $2
cat $inputfile | grep -A 8 "Comparison of alignments to reference annotation" | grep -A 3 'Intron level' | grep 'Precision' | sed -e 's/Precision://g' | tr -d ' \t\r' >> $2

# Base level Sensitivity
echo -ne "sequins-base-level-sensitivity\t" >> $2
cat $inputfile | grep -A 8 "Comparison of alignments to reference annotation" | grep -A 3 'Base level' | grep 'Sensitivity' | sed -e 's/Sensitivity://g' | tr -d ' \t\r' >> $2
 
# Base level Precision
echo -ne "sequins-base-level-precision\t" >> $2
cat $inputfile | grep -A 8 "Comparison of alignments to reference annotation" | grep -A 3 'Base level' | grep 'Precision' | sed -e 's/Precision://g' | tr -d ' \t\r' >> $2
