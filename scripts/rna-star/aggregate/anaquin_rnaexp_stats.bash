# !/bin/bash

inputfile=$1
outputfile=$2

# Isoform Slope
echo -ne "sequins-isoforms-slope\t" > $2
cat $inputfile | grep -A 3 'Linear regression (isoform expression)' | grep 'Slope' | sed -e 's/Slope://g' | tr -d " \t\r" >> $2

# Isoform Pearson Correlation, Log2 Transformed
echo -ne "sequins-isoforms-log2-pearson-cor\t" >> $2
cat $inputfile | grep -A 3 'Linear regression (isoform expression)' | grep 'Correlation' | sed -e 's/Correlation://g' | tr -d " \t\r" >> $2

# Gene Slope
echo -ne "sequins-genes-slope\t" >> $2
cat $inputfile | grep -A 3 'Linear regression (gene expression)' | grep 'Slope' | sed -e 's/Slope://g' | tr -d " \t\r" >> $2

# Gene Pearson Correlation, Log2 Transformed
echo -ne "sequins-genes-log2-pearson-cor\t" >> $2
cat $inputfile | grep -A 3 'Linear regression (gene expression)' | grep 'Correlation' | sed -e 's/Correlation://g' | tr -d " \t\r" >> $2

# Percent Isoforms Found
echo -ne "sequins-percent-isoforms-found\t" >> $2
num1=$(cat anaquin_kallisto/RnaExpression_summary.stats | grep -A 3 "Detected Isoforms" | grep "Sequin" | sed -e 's/Sequin://g' | tr -d " \t\n\r")
denom1=$(cat anaquin_kallisto/RnaExpression_summary.stats | grep -A 3 "Reference Transcript Annotations" | grep 'isoforms' | sed -e 's/Sequin://g' | sed -e 's/isoforms//g' | tr -d " \t\n\r")
echo "scale=5; $num1 / $denom1 " | bc >> $2

# Percent Genes Found
echo -ne "sequins-percent-genes-found\t" >> $2
num2=$(cat $inputfile | grep -A 3 "Detected Genes" | grep "Sequin" | sed -e 's/Sequin://g' | tr -d " \t\n\r")
denom2=$(cat $inputfile | grep -A 3 "Reference Transcript Annotations" | grep 'genes' | sed -e 's/Sequin://g' | sed -e 's/genes//g' | tr -d " \t\n\r")
echo "scale=5; $num2 / $denom2 " | bc >> $2
 
# Detection Sensitivity of Isoforms
echo -ne "sequins-detection-sensitivity-isoforms\t" >> $2
cat $inputfile | grep -A 3 'Detected Isoforms' | grep 'Detection Sensitivity' | sed -e 's/Detection Sensitivity://g' | sed -e 's/(.*//g' | tr -d " \t\r" >> $2

# Detection Sensitivity of Gene
echo -ne "sequins-detection-sensitivity-genes\t" >> $2
cat $inputfile | grep -A 3 'Detected Genes' | grep 'Detection Sensitivity' | sed -e 's/Detection Sensitivity://g' | sed -e 's/(.*//g' | tr -d " \t\r" >> $2



