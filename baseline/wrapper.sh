#!/bin/bash
dir="$(cd "$(dirname "$0")" && pwd)"

testID=$1
species=$2
substitutionModel=$3
mutabilityModel=$4
clonal=$5
fixIndels=$6
region=$7
inputs=$8
inputs=($inputs)
IDs=$9
IDs=($IDs)
ref=${10}
output=${11}
selection=${12}
output_table=${13}
outID="result"

echo "$PWD"

echo "testID = $testID"
echo "species = $species"
echo "substitutionModel = $substitutionModel"
echo "mutabilityModel = $mutabilityModel"
echo "clonal = $clonal"
echo "fixIndels = $fixIndels"
echo "region = $region"
echo "inputs = ${inputs[@]}"
echo "IDs = ${IDs[@]}"
echo "ref = $ref"
echo "output = $output"
echo "outID = $outID"

fasta="$PWD/baseline.fasta"


count=0
for current in ${inputs[@]}
do
	f=$(file $current)
	zipType="Zip archive"
	if [[ "$f" == *"Zip archive"* ]] || [[ "$f" == *"XZ compressed data"* ]]
	then
		id=${IDs[$count]}
		echo "id=$id"
		if [[ "$f" == *"Zip archive"* ]] ; then
			echo "Zip archive"
			echo "unzip $input -d $PWD/files/"
			unzip $current -d "$PWD/$id/"
		elif [[ "$f" == *"XZ compressed data"* ]] ; then
			echo "ZX archive"
			echo "tar -xJf $input -C $PWD/files/"
			mkdir -p "$PWD/$id/files"
			tar -xJf $current -C "$PWD/$id/files/"
		fi
		filtered="$PWD/filtered_${id}.txt"
		imgt_1_file="`find $PWD/$id -name '1_*.txt'`"
		imgt_2_file="`find $PWD/$id -name '2_*.txt'`"
		echo "1_Summary file: ${imgt_1_file}"
		echo "2_IMGT-gapped file: ${imgt_2_file}"
		echo "filter.r for $id"
		Rscript $dir/filter.r ${imgt_1_file} ${imgt_2_file} "$selection" $filtered 2>&1
		
		final="$PWD/final_${id}.txt"
		cat $filtered | cut -f2,4,7 > $final
		python $dir/script_imgt.py --input $final --ref $ref --output $fasta --id $id
	else
		python $dir/script_xlsx.py --input $current --ref $ref --output $fasta
	fi
	count=$((count+1))
done
workdir="$PWD"
cd $dir
echo "file: ${inputs[0]}"
#Rscript --verbose $dir/Baseline_Main.r $testID $species $substitutionModel $mutabilityModel $clonal $fixIndels $region ${inputs[0]} $workdir/ $outID 2>&1
Rscript --verbose $dir/Baseline_Main.r $testID $species $substitutionModel $mutabilityModel $clonal $fixIndels $region $fasta $workdir/ $outID 2>&1

echo "$workdir/${outID}.txt"

rows=`tail -n +2 $workdir/${outID}.txt | grep -v "All sequences combined" | grep -n 'Group' | grep -Eoh '^[0-9]+' | tr '\n' ' '`
rows=($rows)
#unset rows[${#rows[@]}-1]

cd $dir
Rscript --verbose $dir/comparePDFs.r $workdir/${outID}.RData $output ${rows[@]} 2>&1
cp $workdir/result.txt ${output_table}




