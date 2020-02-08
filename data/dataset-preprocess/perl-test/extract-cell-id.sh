#!/bin/sh
INPUT_DIR=$1
echo "INPUT_DIR:$1"
LINE_START=
for filename in $INPUT_DIR/*; do
    # echo "$filename"
    OUTPUT_FILE="$(basename $filename)"
    echo "${OUTPUT_FILE}"
    sed -n '12,1011 p' $filename | cut -f 1 > ./cell-names/$OUTPUT_FILE.cells.txt
done
