# !/bin/bash/

preseq=$1
targets=$2

FRAGMENT_TARGETS=(
      20000000
      30000000
      40000000
     100000000
     150000000
     250000000
)

for target in "${FRAGMENT_TARGETS[@]}" ; do
    fragments_needed=$(awk -v "target=$target" 'NR>1 && $2 > target {printf "%d", $1; exit}' "$preseq")
    if [ -n "$fragments_needed" ]; then
       echo -e "preseq-est-for-$target\t$fragments_needed" >> "$targets"
    fi
done
maximum=$(awk 'END{print int($2)}' "$preseq")
echo -e "preseq-est-max\t$maximum" >> "$targets"