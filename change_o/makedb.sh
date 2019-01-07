#!/bin/bash
dir="$(cd "$(dirname "$0")" && pwd)"

input=$1
noparse=$2
scores=$3
regions=$4
output=$5

if [ "true" == "$noparse" ] ; then
	noparse="--noparse"
else
	noparse=""
fi

if [ "true" == "$scores" ] ; then
	scores="--scores"
else
	scores=""
fi

if [ "true" == "$regions" ] ; then
	regions="--regions"
else
	regions=""
fi

mkdir $PWD/outdir

echo "makedb: $PWD/outdir"

python3 $dir/MakeDb.py imgt -i $input --outdir $PWD/outdir --outname output $noparse $scores $regions
#/data/users/david/anaconda3/bin/python $dir/MakeDb.py imgt -i $input --outdir $PWD/outdir --outname output $noparse $scores $regions
#/home/galaxy/anaconda3/bin/python $dir/MakeDb.py imgt -i $input --outdir $PWD/outdir --outname output $noparse $scores $regions

mv $PWD/outdir/output_db-pass.tab $output

rm -rf $PWD/outdir/
