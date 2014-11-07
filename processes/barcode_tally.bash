# Dependencies
source $MODULELOAD

MAKEFILE="$STAMPIPES/makefiles/barcodes.mk"
bash "$STAMPIPES/scripts/versions.bash" &> barcode.versions.txt

qmake -cwd -N .tb$FLOWCELL -V -q all.q -- -f "$MAKEFILE" -j 40
