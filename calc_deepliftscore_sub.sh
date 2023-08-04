#!/bin/bash

datadir="$2"
datadir2="$3"
fr1="$4"
fr2="$5"
FR1=${fr1^^}
FR2=${fr2^^}
str1="${datadir}Apromoter${fr1}/"
str2="${datadir}${datadir2}"

if [ ! -d $str2 ]; then
	mkdir $str2
fi

tfid="$1"

declare -a array1=("${tfid}_${FR1}" "${tfid}_${FR2}")
declare -a array2=("${datadir}Agenedata${fr1}/test_${tfid}_${FR1}.txt.gz" "${datadir}Agenedata${fr1}/test_${tfid}_${FR2}.txt.gz")
declare -a array3=("${str2}trainout_${tfid}_${FR1}/" "${str2}trainout_${tfid}_${FR2}/")

rm -rf ${str2}trainout_${tfid}_${FR1}/
rm -rf ${str2}trainout_${tfid}_${FR2}/
mkdir -p ${str2}trainout_${tfid}_${FR1}/
mkdir -p ${str2}trainout_${tfid}_${FR2}/
cp -rp ./train_out/* ${str2}trainout_${tfid}_${FR1}/
cp -rp ./train_out/* ${str2}trainout_${tfid}_${FR2}/
read line < ./train_out/summary/best_model.txt
line=`echo ${line} | sed -e "s/[\r\n]\+//g"`
rm -rf ${str2}trainout_${tfid}_${FR1}/${line}
rm -rf ${str2}trainout_${tfid}_${FR2}/${line}
outfile="${str2}calc_deepliftscore_out_${tfid}.txt"

echo ${#array1[@]}
for ((i = 0; i < ${#array1[@]}; i++)) {
	filename1=${array2[i]}
	filename1=${filename1/test/train}
	filename2=${array2[i]}
	filename2=${filename2/test/validate}
	filename3=${array2[i]}
	filename3=${filename3/test/test2}
	echo "zcat ${array2[i]} ${filename1} ${filename2} | sort -u | gzip -9c > ${filename3}"
	zcat ${array2[i]} ${filename1} ${filename2} | sort -u | gzip -9c > ${filename3}
        echo "python3 calc_deepliftscore_sub.py ${str1} ${array1[i]} ${filename3} ${array3[i]} ${outfile}"
        python3 calc_deepliftscore_sub.py ${str1} ${array1[i]} ${filename3} ${array3[i]} ${outfile}
}
