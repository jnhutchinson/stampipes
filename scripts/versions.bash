echo "---"
echo "STAMPIPES:" 
echo "VERSION: " `jq '.version' $STAMPIPES/version.json`
echo "SHA: " `git --git-dir $STAMPIPES/.git rev-parse HEAD`
echo "---"
echo "BEDOPS:"
bedops --help | grep "version" | head -n1
echo "---"
echo "BWA:"
bwa  2>&1 | grep "Version"
echo "---"
echo "SAMTOOLS:"
samtools  2>&1 | grep "Version"
echo "---"
echo "FASTQC:"
fastqc --version
echo "---"
echo "PICARD APPS: "
echo "CollectInsertSizeMetrics"
java -jar `which CollectInsertSizeMetrics.jar`  --version
echo "MarkDuplicates"
java -jar `which MarkDuplicates.jar` --version
echo "---"
echo "AWK:"
awk --version | head -n 1
echo "---"
echo "JQ:"
jq --version
echo "---"
echo "PYTHON:" 
python -V 2>&1
python $STAMPIPES/scripts/versions.py
echo "---"
echo "JAVA:"
java -version 2>&1
echo "---"
echo "BASH:"
bash --version
echo "---"
echo "R:"
R --version | grep "R version"
