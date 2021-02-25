#!/bin/bash
# something to make it commit, with unix line ends
dir="$(cd "$(dirname "$0")" && pwd)"

#define_clones.sh $input $noparse $scores $regions $out_file

type=$1
input=$2

mkdir -p $PWD/outdir

cp $input $PWD/input.tab #file has to have a ".tab" extension

if [ "bygroup" == "$type" ] ; then	
	mode=$3
	act=$4
	model=$5
	norm=$6
	sym=$7
	link=$8
	dist=$9
	output=${10}
	output2=${11}
	
	DefineClones.py -d $PWD/input.tab --nproc 4 --outdir $PWD/outdir --outname output --mode $mode --act $act --model $model --dist $dist --norm $norm --sym $sym --link $link
	
	Rscript $dir/define_clones.r $PWD/outdir/output_clone-pass.tab $output2 2>&1
else
	method=$3
	output=$4
	output2=$5
	
	DefineClones.py hclust -d $PWD/input.tab --nproc 4 --outdir $PWD/outdir --outname output --method $method
	
	Rscript $dir/define_clones.r $PWD/outdir/output_clone-pass.tab $output2 2>&1
fi

cp $PWD/outdir/output_clone-pass.tab $output

rm -rf $PWD/outdir/
