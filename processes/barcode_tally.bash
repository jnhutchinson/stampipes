# Dependencies
source $MODULELOAD
module load python/2.7.3

MAKEFILE="$STAMPIPES/makefiles/barcodes.mk"
bash "$STAMPIPES/scripts/versions.bash" &> barcode.versions.txt

qsub -cwd -N .tb$FLOWCELL -V -q all.q -pe threads 8 <<< make -f "$MAKEFILE" -j 16
