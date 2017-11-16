# !/bin/bash

# output from anaquin RnaExpress
item=$1
outfile=$2
reference=$3

# set list of sequin isoforms
array=(R1_101_1 R1_101_2 R1_102_1 R1_102_2 R1_103_1 R1_103_2 R1_11_1 R1_11_2 R1_12_1 R1_12_2 R1_13_1 R1_13_2 R1_14_1 R1_21_1 R1_21_2 R1_22_1 R1_22_2 R1_23_1 R1_23_2 R1_24_1 R1_31_1 R1_31_2 R1_32_1 R1_32_2 R1_33_1 R1_33_2 R1_41_1 R1_41_2 R1_42_1 R1_42_2 R1_43_1 R1_43_2 R1_51_1 R1_51_2 R1_52_1 R1_52_2 R1_53_1 R1_61_1 R1_61_2 R1_62_1 R1_62_2 R1_63_1 R1_63_2 R1_71_1 R1_71_2 R1_72_1 R1_72_2 R1_73_1 R1_73_2 R1_81_1 R1_81_2 R1_82_1 R1_82_2 R1_83_1 R1_83_2 R1_91_1 R1_91_2 R1_92_1 R1_92_2 R1_93_1 R1_93_2 R2_1_1 R2_105_1 R2_115_1 R2_115_2 R2_116_1 R2_116_2 R2_116_3 R2_117_1 R2_117_3 R2_14_1 R2_14_2 R2_14_3 R2_140_1 R2_143_1 R2_150_1 R2_150_2 R2_151_1 R2_151_2 R2_151_3 R2_152_1 R2_152_2 R2_153_1 R2_153_2 R2_153_3 R2_154_1 R2_154_2 R2_18_1 R2_18_2 R2_19_1 R2_19_2 R2_20_1 R2_20_2 R2_24_1 R2_24_2 R2_26_1 R2_26_2 R2_27_1 R2_27_2 R2_28_1 R2_28_2 R2_28_3 R2_32_1 R2_32_2 R2_32_3 R2_33_1 R2_37_1 R2_37_2 R2_37_3 R2_38_1 R2_38_2 R2_38_3 R2_38_4 R2_41_1 R2_41_2 R2_42_1 R2_42_2 R2_42_3 R2_45_1 R2_45_2 R2_45_3 R2_45_4 R2_46_1 R2_46_2 R2_46_3 R2_47_1 R2_47_2 R2_53_1 R2_53_2 R2_53_3 R2_54_1 R2_54_2 R2_55_2 R2_55_3 R2_57_1 R2_57_2 R2_59_1 R2_59_2 R2_59_3 R2_6_1 R2_6_2 R2_60_1 R2_60_2 R2_63_1 R2_63_3 R2_65_1 R2_66_1 R2_66_2 R2_67_1 R2_68_1 R2_68_2 R2_7_1 R2_7_2 R2_71_1 R2_71_2 R2_72_1 R2_72_2 R2_72_3 R2_72_4 R2_73_1 R2_73_2 R2_76_1 R2_76_2 R2_76_3)

#
for term in ${array[*]}
do
    saveval=$(cat *$item* | grep $term | awk '{print $4}')
    if [ -z $saveval ]
    then
        echo "0" >> $outfile
    else
        echo "$saveval" >> $outfile
    fi
done

paste $reference $outfile > tmp.txt
mv tmp.txt $outfile

Rscript $STAMPIPES/scripts/rna-star/aggregate/anaquin_neat_comparison.Rscript $outfile $outfile.info