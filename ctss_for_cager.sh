#!/bin/bash

for file in ./*.R1.ctss.bed; do
	name=$(basename $file .R1.ctss.bed)
	cat $file | awk 'BEGIN{OFS="\t"}{if($5 != 0) print $1,$2,$6,$5}' > $name.cager.ctss.bed
done

