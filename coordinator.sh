#!/bin/bash
title="title"
search_specs="search_specs.txt"
t="$(grep $title $search_specs)"
t=${t#*\"}
t=${t%\"*}
mkdir $t$3
_d="_data.txt"
touch "./$t$_d"
for ((i=0;i<$1;i++))
do
	sbatch nwbh.sh $i 0 $2 $1 $3
done
