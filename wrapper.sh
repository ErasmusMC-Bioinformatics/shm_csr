#!/usr/bin/env bash
set -e -o pipefail
dir="$(cd "$(dirname "$0")" && pwd)"
input=$1
method=$2
log=$3 #becomes the main html page at the end
outdir=$4
output="$outdir/index.html" #copied to $log location at the end
title="$5"
include_fr1=$6
functionality=$7
unique=$8
naive_output=$9
naive_output_ca=${10}
naive_output_cg=${11}
naive_output_cm=${12}
naive_output_ce=${13}
naive_output_all=${14}
filter_unique=${15}
filter_unique_count=${16}
class_filter=${17}
empty_region_filter=${18}
fast=${19}

#exec 5> debug_output.txt
#BASH_XTRACEFD="5"
## Busybox date does not support '+%s.%N'. So use the slower python instead.
## Using -S python does not do 'import site' which shortens the command
## to 10 milliseconds.
#PS4='$(python -Sc "import time; print(time.time())") $LINENO: '
#set -x

mkdir -p $outdir

tar -xzf $dir/style.tar.gz -C $outdir

echo "---------------- read parameters ----------------"
echo "---------------- read parameters ----------------<br />" > $log

echo "unpacking IMGT file"

type="`file $input`"
if [[ "$type" == *"Zip archive"* ]] ; then
	echo "Zip archive"
	echo "unzip $input -d $PWD/files/"
	unzip $input -d $PWD/files/
elif [[ "$type" == *"XZ compressed data"* ]] ; then
	echo "ZX archive"
	echo "tar -xJf $input -C $PWD/files/"
	mkdir -p "$PWD/files/$title"
	tar -xJf $input -C "$PWD/files/$title"
else
	echo "Unrecognized format $type"
	echo "Unrecognized format $type" > $log
	exit 1
fi

cat "`find $PWD/files/ -name "1_*"`" > $PWD/summary.txt
cat "`find $PWD/files/ -name "2_*"`" > $PWD/gapped_nt.txt
cat "`find $PWD/files/ -name "3_*"`" > $PWD/sequences.txt
cat "`find $PWD/files/ -name "4_*"`" > $PWD/gapped_aa.txt
cat "`find $PWD/files/ -name "5_*"`" > $PWD/aa.txt
cat "`find $PWD/files/ -name "6_*"`" > $PWD/junction.txt
cat "`find $PWD/files/ -name "7_*"`" > $PWD/mutationanalysis.txt
cat "`find $PWD/files/ -name "8_*"`" > $PWD/mutationstats.txt
cat "`find $PWD/files/ -name "9_*"`" > $PWD/aa_change_stats.txt
cat "`find $PWD/files/ -name "10_*"`" > $PWD/hotspots.txt

echo "---------------- unique id check ----------------"

Rscript $dir/check_unique_id.r $PWD/summary.txt $PWD/gapped_nt.txt $PWD/sequences.txt $PWD/gapped_aa.txt $PWD/aa.txt $PWD/junction.txt $PWD/mutationanalysis.txt $PWD/mutationstats.txt $PWD/aa_change_stats.txt $PWD/hotspots.txt

if [[ ${#BLASTN_DIR} -ge 5 ]] ; then
	echo "On server, using BLASTN_DIR env: ${BLASTN_DIR}"
else
	BLASTN_DIR="/home/galaxy/Downloads/ncbi-blast-2.4.0+/bin"
	echo "Dev Galaxy set BLASTN_DIR to: ${BLASTN_DIR}"
fi

echo "---------------- class identification ----------------"
echo "---------------- class identification ----------------<br />" >> $log

python $dir/gene_identification.py --input $PWD/summary.txt --output $outdir/identified_genes.txt

echo "---------------- merge_and_filter.r ----------------"
echo "---------------- merge_and_filter.r ----------------<br />" >> $log

Rscript $dir/merge_and_filter.r $PWD/summary.txt $PWD/sequences.txt $PWD/mutationanalysis.txt $PWD/mutationstats.txt $PWD/hotspots.txt "$PWD/gapped_aa.txt" $outdir/identified_genes.txt $outdir/merged.txt $outdir/before_unique_filter.txt $outdir/unmatched.txt $method $functionality $unique ${filter_unique} ${filter_unique_count} ${class_filter} ${empty_region_filter} 2>&1

if [[ "${naive_output}" == "yes" ]] || [[ "$fast" == "no" ]] ; then

	echo "---------------- creating new IMGT zips ----------------"
	echo "---------------- creating new IMGT zips ----------------<br />" >> $log

	python $dir/split_imgt_file.py --outdir $outdir $input $outdir/merged.txt \
	  - IGA IGA1 IGA2 IGG IGG1 IGG2 IGG3 IGG4 IGM IGE

fi

echo "---------------- shm_csr.r ----------------"
echo "---------------- shm_csr.r ----------------<br />" >> $log

classes="IGA,IGA1,IGA2,IGG,IGG1,IGG2,IGG3,IGG4,IGM,IGE,unmatched"
echo "R mutation analysis"
Rscript $dir/shm_csr.r $outdir/merged.txt $classes $outdir ${empty_region_filter} 2>&1

echo "---------------- plot_pdfs.r ----------------"
echo "---------------- plot_pdfs.r ----------------<br />" >> $log

echo "Rscript $dir/shm_csr.r $outdir/pdfplots.RData $outdir 2>&1"

Rscript $dir/plot_pdf.r "$outdir/pdfplots.RData" "$outdir" 2>&1

echo "---------------- shm_csr.py ----------------"
echo "---------------- shm_csr.py ----------------<br />" >> $log

python $dir/shm_csr.py --input $outdir/merged.txt --genes $classes --empty_region_filter "${empty_region_filter}" --output $outdir/hotspot_analysis.txt

echo "---------------- aa_histogram.r ----------------"
echo "---------------- aa_histogram.r ----------------<br />" >> $log

Rscript $dir/aa_histogram.r $outdir/aa_id_mutations.txt $outdir/absent_aa_id.txt "IGA,IGG,IGM,IGE" $outdir/ 2>&1
if [ -e "$outdir/aa_histogram_.png" ]; then
        mv $outdir/aa_histogram_.png $outdir/aa_histogram.png
        mv $outdir/aa_histogram_.pdf $outdir/aa_histogram.pdf
        mv $outdir/aa_histogram_.txt $outdir/aa_histogram.txt
        mv $outdir/aa_histogram_absent_.txt $outdir/aa_histogram_absent.txt
        mv $outdir/aa_histogram_count_.txt $outdir/aa_histogram_count.txt
        mv $outdir/aa_histogram_sum_.txt $outdir/aa_histogram_sum.txt
fi

genes=(IGA IGA1 IGA2 IGG IGG1 IGG2 IGG3 IGG4 IGM IGE)

funcs=(sum mean median)
funcs=(sum)

echo "---------------- sequence_overview.r ----------------"
echo "---------------- sequence_overview.r ----------------<br />" >> $log

mkdir $outdir/sequence_overview

Rscript $dir/sequence_overview.r $outdir/before_unique_filter.txt $outdir/merged.txt $outdir/sequence_overview $classes $outdir/hotspot_analysis_sum.txt ${empty_region_filter} 2>&1

echo "<table border='1'>" > $outdir/base_overview.html

while IFS=$'\t' read ID class seq A C G T
do
	echo "<tr><td>$ID</td><td>$seq</td><td>$class</td><td>$A</td><td>$C</td><td>$G</td><td>$T</td></tr>" >> $outdir/base_overview.html
done < $outdir/sequence_overview/ntoverview.txt

echo "<html><center><h1>$title</h1></center>" > $output
echo "<meta name='viewport' content='width=device-width, initial-scale=1'>" >> $output
echo "<script type='text/javascript' src='jquery-1.11.0.min.js'></script>" >> $output
echo "<script type='text/javascript' src='tabber.js'></script>" >> $output
echo "<script type='text/javascript' src='script.js'></script>" >> $output
echo "<link rel='stylesheet' type='text/css' href='style.css'>" >> $output
echo "<link rel='stylesheet' type='text/css' href='pure-min.css'>" >> $output

matched_count="`cat $outdir/merged.txt | grep -v 'unmatched' | tail -n +2 | wc -l`"
unmatched_count="`cat $outdir/unmatched.txt | tail -n +2 | wc -l`"
total_count=$((matched_count + unmatched_count))
perc_count=$((unmatched_count / total_count * 100))
perc_count=`bc -l <<< "scale=2; ${unmatched_count} / ${total_count} * 100"`
perc_count=`bc -l <<< "scale=2; (${unmatched_count} / ${total_count} * 100 ) / 1"`

echo "<center><h2>Total: ${total_count}</h2></center>" >> $output
echo "<center><h2>Matched: ${matched_count} Unmatched: ${unmatched_count}</h2></center>" >> $output
echo "<center><h2>Percentage unmatched: ${perc_count}</h2></center>" >> $output

echo "---------------- main tables ----------------"
echo "---------------- main tables ----------------<br />" >> $log

echo "<div class='tabber'>" >> $output
echo "<div class='tabbertab' title='SHM Overview' style='width: 3000px;'>" >> $output

for func in ${funcs[@]}
do
	
	echo "---------------- $func table ----------------"
	echo "---------------- $func table ----------------<br />" >> $log
	
	cat $outdir/mutations_${func}.txt $outdir/shm_overview_tandem_row.txt $outdir/hotspot_analysis_${func}.txt > $outdir/data_${func}.txt
	
	echo "---------------- pattern_plots.r ----------------"
	echo "---------------- pattern_plots.r ----------------<br />" >> $log

	Rscript $dir/pattern_plots.r $outdir/data_${func}.txt $outdir/aid_motives $outdir/relative_mutations $outdir/absolute_mutations $outdir/shm_overview.txt 2>&1
	
	echo "<table class='pure-table pure-table-striped'>" >> $output
	echo "<thead><tr><th>info</th>" >> $output
	
	if [ "${class_filter}" != "101_101" ] ; then
	
		for gene in ${genes[@]}
		do
			tmp=`cat $outdir/${gene}_${func}_n.txt`
			echo "<th><a href='matched_${gene}_${func}.txt'>${gene} (N = $tmp)</a></th>" >> $output
		done
		
		tmp=`cat $outdir/all_${func}_n.txt`
		echo "<th><a href='matched_all_${func}.txt'>all (N = $tmp)</a></th>" >> $output
		tmp=`cat $outdir/unmatched_${func}_n.txt`
		echo "<th><a href='unmatched.txt'>unmatched (N = ${unmatched_count})</a></th><tr></thead>" >> $output

		while IFS=, read name cax cay caz ca1x ca1y ca1z ca2x ca2y ca2z cgx cgy cgz cg1x cg1y cg1z cg2x cg2y cg2z cg3x cg3y cg3z cg4x cg4y cg4z cmx cmy cmz cex cey cez unx uny unz allx ally allz 
		do
			if [ "$name" == "FR R/S (ratio)" ] || [ "$name" == "CDR R/S (ratio)" ] || [ "$name" == "Tandems/Expected (ratio)" ] ; then #meh
				echo "<tr><td>$name</td><td>${cax}/${cay} (${caz})</td><td>${ca1x}/${ca1y} (${ca1z})</td><td>${ca2x}/${ca2y} (${ca2z})</td><td>${cgx}/${cgy} (${cgz})</td><td>${cg1x}/${cg1y} (${cg1z})</td><td>${cg2x}/${cg2y} (${cg2z})</td><td>${cg3x}/${cg3y} (${cg3z})</td><td>${cg4x}/${cg4y} (${cg4z})</td><td>${cmx}/${cmy} (${cmz})</td><td>${cex}/${cey} (${cez})</td><td>${allx}/${ally} (${allz})</td><td>${unx}/${uny} (${unz})</td></tr>" >> $output
			elif [ "$name" == "Median of Number of Mutations (%)" ] ; then
				echo "<tr><td>$name</td><td>${caz}%</td><td>${ca1z}%</td><td>${ca2z}%</td><td>${cgz}%</td><td>${cg1z}%</td><td>${cg2z}%</td><td>${cg3z}%</td><td>${cg4z}%</td><td>${cmz}%</td><td>${cez}%</td><td>${allz}%</td><td>${unz}%</td></tr>" >> $output
			else
				echo "<tr><td>$name</td><td>${cax}/${cay} (${caz}%)</td><td>${ca1x}/${ca1y} (${ca1z}%)</td><td>${ca2x}/${ca2y} (${ca2z}%)</td><td>${cgx}/${cgy} (${cgz}%)</td><td>${cg1x}/${cg1y} (${cg1z}%)</td><td>${cg2x}/${cg2y} (${cg2z}%)</td><td>${cg3x}/${cg3y} (${cg3z}%)</td><td>${cg4x}/${cg4y} (${cg4z}%)</td><td>${cmx}/${cmy} (${cmz}%)</td><td>${cex}/${cey} (${cez}%)</td><td>${allx}/${ally} (${allz}%)</td><td>${unx}/${uny} (${unz}%)</td></tr>" >> $output
			fi
		done < $outdir/data_${func}.txt
		
	else
		tmp=`cat $outdir/all_${func}_n.txt`
		echo "<th><a href='matched_all_${func}.txt'>all (N = $tmp)</a></th>" >> $output
		
		while IFS=, read name cax cay caz ca1x ca1y ca1z ca2x ca2y ca2z cgx cgy cgz cg1x cg1y cg1z cg2x cg2y cg2z cg3x cg3y cg3z cg4x cg4y cg4z cmx cmy cmz cex cey cez unx uny unz allx ally allz
		do
			if [ "$name" == "FR R/S (ratio)" ] || [ "$name" == "CDR R/S (ratio)" ] ; then #meh
				echo "<tr><td>$name</td><td>${allx}/${ally}</td></tr>" >> $output
			elif [ "$name" == "Median of Number of Mutations (%)" ] ; then
				echo "<tr><td>$name</td><td>${allz}%</td></tr>" >> $output
			else
				echo "<tr><td>$name</td><td>${allx}/${ally} (${allz}%)</td></tr>" >> $output
			fi
		done < $outdir/data_${func}.txt
		
	fi
	echo "</table>" >> $output
	#echo "<a href='data_${func}.txt'>Download data</a>" >> $output
done

echo "<a href='aid_motives.pdf'><img src='aid_motives.png' /></a><br />" >> $output
echo "<a href='relative_mutations.pdf'><img src='relative_mutations.png' /></a><br />" >> $output
echo "<a href='absolute_mutations.pdf'><img src='absolute_mutations.png' /></a><br />" >> $output
echo "<br />" >> $output
cat $dir/shm_overview.htm >> $output
echo "</div>" >> $output #SHM overview tab end

echo "---------------- images ----------------"
echo "---------------- images ----------------<br />" >> $log

echo "<div class='tabbertab' title='SHM Frequency' style='width: 3000px;'></a>" >> $output

if [ -a $outdir/scatter.png ]
then
	echo "<a href='scatter.pdf'><img src='scatter.png'/><br />" >> $output
fi
if [ -a $outdir/frequency_ranges.png ]
then
	echo "<a href='frequency_ranges.pdf'><img src='frequency_ranges.png'/></a><br />" >> $output
fi

echo "<br />" >> $output
cat $dir/shm_frequency.htm >> $output

echo "</div>" >> $output #SHM frequency tab end

echo "<div class='tabbertab' title='Transition tables' style='width: 3000px;'>" >> $output

echo "<table border='0'>" >> $output

for gene in ${genes[@]}
do
	echo "<tr>" >> $output
	echo "<td><h1>${gene}</h1></td>" >> $output
	
	if [ -e $outdir/transitions_heatmap_${gene}.png ]
	then
		echo "<td><a href='transitions_heatmap_${gene}.pdf'><img src='transitions_heatmap_${gene}.png' /></a></td>" >> $output
	else
		echo "<td></td>" >> $output
	fi
	
	if [ -e $outdir/transitions_stacked_${gene}.png ]
	then
		echo "<td><a href='transitions_stacked_${gene}.pdf'><img src='transitions_stacked_${gene}.png' /></a></td>" >> $output
	else
		echo "<td></td>" >> $output
	fi
	
	echo "<td><table style='border-left-width: 1;' class='pure-table transition-table pure-table-bordered'>" >> $output
	echo "<tr><td></td><td colspan="5"><center>To</center></td></tr>" >> $output
	first="true"
	while IFS=, read from a c g t
		do
			if [ "$first" == "true" ] ; then
				echo "<tr><td rowspan='5'>From</td><td>$from</td><td>$a</td><td>$c</td><td>$g</td><td>$t</td></tr>" >> $output
				first="false"
			else
				echo "<tr><td>$from</td><td>$a</td><td>$c</td><td>$g</td><td>$t</td></tr>" >> $output
			fi
	done < $outdir/transitions_${gene}_sum.txt
	echo "</table></td>" >> $output
	
	echo "</tr>" >> $output
done

echo "<tr>" >> $output
echo "<td><h1>All</h1></td>" >> $output
echo "<td><a href='transitions_heatmap_all.pdf'><img src='transitions_heatmap_all.png' /></a></td>" >> $output
echo "<td><a href='transitions_stacked_all.pdf'><img src='transitions_stacked_all.png' /></a></td>" >> $output
echo "<td><table style='border-left-width: 1;' class='pure-table transition-table pure-table-bordered'>" >> $output
echo "<tr><td></td><td colspan="5"><center>To</center></td></tr>" >> $output
first="true"
while IFS=, read from a c g t
	do
		if [ "$first" == "true" ] ; then
			echo "<tr><td rowspan='5'>From</td><td>$from</td><td>$a</td><td>$c</td><td>$g</td><td>$t</td></tr>" >> $output
			first="false"
		else
			echo "<tr><td>$from</td><td>$a</td><td>$c</td><td>$g</td><td>$t</td></tr>" >> $output
		fi
done < $outdir/transitions_all_sum.txt
echo "</table></td>" >> $output

echo "</tr>" >> $output

echo "</table>" >> $output

echo "<br />" >> $output
cat $dir/shm_transition.htm >> $output

echo "</div>" >> $output #transition tables tab end

echo "<div class='tabbertab' title='Antigen Selection'>" >> $output

if [ -e $outdir/aa_histogram.png ]
then
	echo "<a href='aa_histogram.pdf'><img src='aa_histogram.png'/></a><br />" >> $output
fi

if [ -e $outdir/aa_histogram_IGA.png ]
then
	echo "<a href='aa_histogram_IGA.pdf'><img src='aa_histogram_IGA.png'/></a><br />" >> $output
fi

if [ -e $outdir/aa_histogram_IGG.png ]
then
	echo "<a href='aa_histogram_IGG.pdf'><img src='aa_histogram_IGG.png'/></a><br />" >> $output
fi

if [ -e $outdir/aa_histogram_IGM.png ]
then
	echo "<a href='aa_histogram_IGM.pdf'><img src='aa_histogram_IGM.png'/></a><br />" >> $output
fi

if [ -e $outdir/aa_histogram_IGE.png ]
then
	echo "<a href='aa_histogram_IGE.pdf'><img src='aa_histogram_IGE.png'/></a><br />" >> $output
fi



if [[ "$fast" == "no" ]] ; then

    

	echo "---------------- baseline ----------------"
	echo "---------------- baseline ----------------<br />" >> $log
	tmp="$PWD"

	mkdir -p $outdir/baseline
	
	echo "<center><h1>BASELINe</h1>" >> $output
	header_substring="Based on CDR1, FR2, CDR2, FR3 (27:27:38:55:65:104:-)"
	
	baseline_boundaries="27:27:38:55:65:104:-"
	
	if [[ "${empty_region_filter}" == "leader" ]] ; then
		baseline_boundaries="1:26:38:55:65:104:-"
		header_substring="Based on FR1, CDR1, FR2, CDR2, FR3 (1:26:38:55:65:104,-)"
	fi
	
	echo "<p>${header_substring}</p></center>" >> $output

	mkdir $outdir/baseline/IGA_IGG_IGM
	if [[ $(wc -l < $outdir/new_IMGT/1_Summary.txt) -gt "1" ]]; then
		cd $outdir/baseline/IGA_IGG_IGM
		bash $dir/baseline/wrapper.sh 1 1 1 1 0 0 "${baseline_boundaries}" $outdir/new_IMGT.txz "IGA_IGG_IGM_IGE" "$dir/baseline/IMGTVHreferencedataset20161215.fa" "$outdir/baseline.pdf" "Sequence.ID" "$outdir/baseline.txt"
	else
		echo "No sequences" > "$outdir/baseline.txt"
	fi

	mkdir $outdir/baseline/IGA
	if [[ $(wc -l < $outdir/new_IMGT_IGA/1_Summary.txt) -gt "1" ]]; then
		cd $outdir/baseline/IGA
		bash $dir/baseline/wrapper.sh 1 1 1 1 0 0 "${baseline_boundaries}" $outdir/new_IMGT_IGA.txz "IGA" "$dir/baseline/IMGTVHreferencedataset20161215.fa" "$outdir/baseline_IGA.pdf" "Sequence.ID" "$outdir/baseline_IGA.txt"
	else
		echo "No IGA sequences" > "$outdir/baseline_IGA.txt"
	fi

	mkdir $outdir/baseline/IGG
	if [[ $(wc -l < $outdir/new_IMGT_IGG/1_Summary.txt) -gt "1" ]]; then
		cd $outdir/baseline/IGG
		bash $dir/baseline/wrapper.sh 1 1 1 1 0 0 "${baseline_boundaries}" $outdir/new_IMGT_IGG.txz "IGG" "$dir/baseline/IMGTVHreferencedataset20161215.fa" "$outdir/baseline_IGG.pdf" "Sequence.ID" "$outdir/baseline_IGG.txt"
	else
		echo "No IGG sequences" > "$outdir/baseline_IGG.txt"
	fi

	mkdir $outdir/baseline/IGM
	if [[ $(wc -l < $outdir/new_IMGT_IGM/1_Summary.txt) -gt "1" ]]; then
		cd $outdir/baseline/IGM
		bash $dir/baseline/wrapper.sh 1 1 1 1 0 0 "${baseline_boundaries}" $outdir/new_IMGT_IGM.txz "IGM" "$dir/baseline/IMGTVHreferencedataset20161215.fa" "$outdir/baseline_IGM.pdf" "Sequence.ID" "$outdir/baseline_IGM.txt"
	else
		echo "No IGM sequences" > "$outdir/baseline_IGM.txt"
	fi

	mkdir $outdir/baseline/IGE
	if [[ $(wc -l < $outdir/new_IMGT_IGE/1_Summary.txt) -gt "1" ]]; then
		cd $outdir/baseline/IGE
		bash $dir/baseline/wrapper.sh 1 1 1 1 0 0 "${baseline_boundaries}" $outdir/new_IMGT_IGE.txz "IGE" "$dir/baseline/IMGTVHreferencedataset20161215.fa" "$outdir/baseline_IGE.pdf" "Sequence.ID" "$outdir/baseline_IGE.txt"
	else
		echo "No IGE sequences" > "$outdir/baseline_IGE.txt"
	fi

	cd $tmp

	echo "Cleaning up *.RData files"
	find $outdir/baseline -name "*.RData" -type f -delete
	
	if [ -e $outdir/baseline.pdf ]
	then
		echo "<embed src='baseline.pdf' width='700px' height='1000px'>" >> $output
	fi

	if [ -e $outdir/baseline_IGA.pdf ]
	then
		echo "<embed src='baseline_IGA.pdf' width='700px' height='1000px'>" >> $output
	fi

	if [ -e $outdir/baseline_IGG.pdf ]
	then
		echo "<embed src='baseline_IGG.pdf' width='700px' height='1000px'>" >> $output
	fi

	if [ -e $outdir/baseline_IGM.pdf ]
	then
		echo "<embed src='baseline_IGM.pdf' width='700px' height='1000px'>" >> $output
	fi

	if [ -e $outdir/baseline_IGE.pdf ]
	then
		echo "<embed src='baseline_IGE.pdf' width='700px' height='1000px'>" >> $output
	fi
fi

echo "<br />" >> $output
cat $dir/shm_selection.htm >> $output

echo "</div>" >> $output #antigen selection tab end

echo "<div class='tabbertab' title='CSR'>" >> $output #CSR tab

if [ -e $outdir/IGA.png ] 
then
	echo "<a href='IGA.pdf'><img src='IGA.png'/></a><br />" >> $output
fi
if [ -e $outdir/IGG.png ]
then
	echo "<a href='IGG.pdf'><img src='IGG.png'/></a><br />" >> $output
fi

echo "<br />" >> $output
cat $dir/shm_csr.htm >> $output

echo "</div>" >> $output #CSR tab end

if [[ "$fast" == "no" ]] ; then

	echo "---------------- change-o MakeDB ----------------"

	mkdir -p $outdir/change_o

	tmp="$PWD"

	cd $outdir/change_o

	bash $dir/change_o/makedb.sh $outdir/new_IMGT.txz false false false $outdir/change_o/change-o-db.txt
	bash $dir/change_o/define_clones.sh bygroup $outdir/change_o/change-o-db.txt gene first ham none min complete 3.0 $outdir/change_o/change-o-db-defined_clones.txt $outdir/change_o/change-o-defined_clones-summary.txt
	Rscript $dir/change_o/select_first_in_clone.r $outdir/change_o/change-o-db-defined_clones.txt $outdir/change_o/change-o-db-defined_first_clones.txt 2>&1
	
	mkdir $outdir/new_IMGT_changeo
	cp $outdir/new_IMGT/* $outdir/new_IMGT_changeo
	
	Rscript $dir/new_imgt.r $outdir/new_IMGT_changeo $outdir/change_o/change-o-db-defined_first_clones.txt "-" 2>&1
	
	cd $outdir/new_IMGT_changeo
	tar -cJf ../new_IMGT_first_seq_of_clone.txz *
	cd $outdir/change_o
	
	rm -rf $outdir/new_IMGT_changeo
	
	Rscript $dir/merge.r $outdir/change_o/change-o-db-defined_clones.txt $outdir/merged.txt "all" "Sequence.ID,best_match" "SEQUENCE_ID" "Sequence.ID" $outdir/change_o/change-o-db-defined_clones.txt 2>&1
	echo "Rscript $dir/merge.r $outdir/change_o/change-o-db-defined_clones.txt $outdir/$outdir/merged.txt 'all' 'Sequence.ID,best_match' 'Sequence.ID' 'Sequence.ID' '\t' $outdir/change_o/change-o-db-defined_clones.txt 2>&1"
	
	if [[ $(wc -l < $outdir/new_IMGT_IGA/1_Summary.txt) -gt "1" ]]; then
		bash $dir/change_o/makedb.sh $outdir/new_IMGT_IGA.txz false false false $outdir/change_o/change-o-db-IGA.txt
		bash $dir/change_o/define_clones.sh bygroup $outdir/change_o/change-o-db-IGA.txt gene first ham none min complete 3.0 $outdir/change_o/change-o-db-defined_clones-IGA.txt $outdir/change_o/change-o-defined_clones-summary-IGA.txt
		Rscript $dir/change_o/select_first_in_clone.r $outdir/change_o/change-o-db-defined_clones-IGA.txt $outdir/change_o/change-o-db-defined_first_clones-IGA.txt 2>&1
		
		mkdir $outdir/new_IMGT_IGA_changeo
		cp $outdir/new_IMGT/* $outdir/new_IMGT_IGA_changeo
		
		Rscript $dir/new_imgt.r $outdir/new_IMGT_IGA_changeo $outdir/change_o/change-o-db-defined_first_clones-IGA.txt "-" 2>&1
		
		cd $outdir/new_IMGT_IGA_changeo
		tar -cJf ../new_IMGT_IGA_first_seq_of_clone.txz *
		
		rm -rf $outdir/new_IMGT_IGA_changeo
		
		cd $outdir/change_o
	else
		echo "No IGA sequences" > "$outdir/change_o/change-o-db-defined_clones-IGA.txt"
		echo "No IGA sequences" > "$outdir/change_o/change-o-defined_clones-summary-IGA.txt"
	fi
	
	if [[ $(wc -l < $outdir/new_IMGT_IGG/1_Summary.txt) -gt "1" ]]; then
		bash $dir/change_o/makedb.sh $outdir/new_IMGT_IGG.txz false false false $outdir/change_o/change-o-db-IGG.txt
		bash $dir/change_o/define_clones.sh bygroup $outdir/change_o/change-o-db-IGG.txt gene first ham none min complete 3.0 $outdir/change_o/change-o-db-defined_clones-IGG.txt $outdir/change_o/change-o-defined_clones-summary-IGG.txt
		Rscript $dir/change_o/select_first_in_clone.r $outdir/change_o/change-o-db-defined_clones-IGG.txt $outdir/change_o/change-o-db-defined_first_clones-IGG.txt 2>&1
		
		mkdir $outdir/new_IMGT_IGG_changeo
		cp $outdir/new_IMGT/* $outdir/new_IMGT_IGG_changeo
		
		Rscript $dir/new_imgt.r $outdir/new_IMGT_IGG_changeo $outdir/change_o/change-o-db-defined_first_clones-IGG.txt "-" 2>&1
		
		cd $outdir/new_IMGT_IGG_changeo
		tar -cJf ../new_IMGT_IGG_first_seq_of_clone.txz *
		rm -rf $outdir/new_IMGT_IGG_changeo
		
		cd $outdir/change_o
	else
		echo "No IGG sequences" > "$outdir/change_o/change-o-db-defined_clones-IGG.txt"
		echo "No IGG sequences" > "$outdir/change_o/change-o-defined_clones-summary-IGG.txt"
	fi

	if [[ $(wc -l < $outdir/new_IMGT_IGM/1_Summary.txt) -gt "1" ]]; then
		bash $dir/change_o/makedb.sh $outdir/new_IMGT_IGM.txz false false false $outdir/change_o/change-o-db-IGM.txt
		bash $dir/change_o/define_clones.sh bygroup $outdir/change_o/change-o-db-IGM.txt gene first ham none min complete 3.0 $outdir/change_o/change-o-db-defined_clones-IGM.txt $outdir/change_o/change-o-defined_clones-summary-IGM.txt
		Rscript $dir/change_o/select_first_in_clone.r $outdir/change_o/change-o-db-defined_clones-IGM.txt $outdir/change_o/change-o-db-defined_first_clones-IGM.txt 2>&1
		
		mkdir $outdir/new_IMGT_IGM_changeo
		cp $outdir/new_IMGT/* $outdir/new_IMGT_IGM_changeo
		
		Rscript $dir/new_imgt.r $outdir/new_IMGT_IGM_changeo $outdir/change_o/change-o-db-defined_first_clones-IGM.txt "-" 2>&1
		
		cd $outdir/new_IMGT_IGM_changeo
		tar -cJf ../new_IMGT_IGM_first_seq_of_clone.txz *
		
		rm -rf $outdir/new_IMGT_IGM_changeo
		
		cd $outdir/change_o
	else
		echo "No IGM sequences" > "$outdir/change_o/change-o-db-defined_clones-IGM.txt"
		echo "No IGM sequences" > "$outdir/change_o/change-o-defined_clones-summary-IGM.txt"
	fi

	if [[ $(wc -l < $outdir/new_IMGT_IGE/1_Summary.txt) -gt "1" ]]; then
		bash $dir/change_o/makedb.sh $outdir/new_IMGT_IGE.txz false false false $outdir/change_o/change-o-db-IGE.txt
		bash $dir/change_o/define_clones.sh bygroup $outdir/change_o/change-o-db-IGE.txt gene first ham none min complete 3.0 $outdir/change_o/change-o-db-defined_clones-IGE.txt $outdir/change_o/change-o-defined_clones-summary-IGE.txt
		Rscript $dir/change_o/select_first_in_clone.r $outdir/change_o/change-o-db-defined_clones-IGE.txt $outdir/change_o/change-o-db-defined_first_clones-IGE.txt 2>&1
		
		mkdir $outdir/new_IMGT_IGE_changeo
		cp $outdir/new_IMGT/* $outdir/new_IMGT_IGE_changeo
		
		Rscript $dir/new_imgt.r $outdir/new_IMGT_IGE_changeo $outdir/change_o/change-o-db-defined_first_clones-IGE.txt "-" 2>&1
		
		cd $outdir/new_IMGT_IGE_changeo
		tar -cJf ../new_IMGT_IGE_first_seq_of_clone.txz *
		
		rm -rf $outdir/new_IMGT_IGE_changeo
		
		cd $outdir/change_o
	else
		echo "No IGE sequences" > "$outdir/change_o/change-o-db-defined_clones-IGE.txt"
		echo "No IGE sequences" > "$outdir/change_o/change-o-defined_clones-summary-IGE.txt"
	fi

	cd "$tmp"
	
	rm -rf $outdir/new_IMGT
	rm -rf $outdir/new_IMGT_IGA/
	rm -rf $outdir/new_IMGT_IGA1/
	rm -rf $outdir/new_IMGT_IGA2/
	rm -rf $outdir/new_IMGT_IGG/
	rm -rf $outdir/new_IMGT_IGG1/
	rm -rf $outdir/new_IMGT_IGG2/
	rm -rf $outdir/new_IMGT_IGG3/
	rm -rf $outdir/new_IMGT_IGG4/
	rm -rf $outdir/new_IMGT_IGM/
	rm -rf $outdir/new_IMGT_IGE/

	echo "<div class='tabbertab' title='Clonal Relation' style='width: 7000px;'>" >> $output #clonality tab

	function clonality_table {
		local infile=$1
		local outfile=$2
		
		echo "<table class='pure-table pure-table-striped'>" >> $outfile
		echo "<thead><tr><th>Clone size</th><th>Nr of clones</th><th>Nr of sequences</th></tr></thead>" >> $outfile
		
		first='true'
		
		while read size clones seqs
		do
			if [[ "$first" == "true" ]]; then
				first="false"
				continue
			fi
			echo "<tr><td>$size</td><td>$clones</td><td>$seqs</td></tr>" >> $outfile
		done < $infile
		
		echo "</table>" >> $outfile
	}
	echo "<div class='tabber'>" >> $output

	echo "<div class='tabbertab' title='All'>" >> $output
	clonality_table $outdir/change_o/change-o-defined_clones-summary.txt $output
	echo "</div>" >> $output

	echo "<div class='tabbertab' title='IGA'>" >> $output
	clonality_table $outdir/change_o/change-o-defined_clones-summary-IGA.txt $output
	echo "</div>" >> $output

	echo "<div class='tabbertab' title='IGG'>" >> $output
	clonality_table $outdir/change_o/change-o-defined_clones-summary-IGG.txt $output
	echo "</div>" >> $output

	echo "<div class='tabbertab' title='IGM'>" >> $output
	clonality_table $outdir/change_o/change-o-defined_clones-summary-IGM.txt $output
	echo "</div>" >> $output

	echo "<div class='tabbertab' title='IGE'>" >> $output
	clonality_table $outdir/change_o/change-o-defined_clones-summary-IGM.txt $output
	echo "</div>" >> $output

	echo "<div class='tabbertab' title='Overlap' style='width: 7000px;'>" >> $output
	cat "$outdir/sequence_overview/index.html" | sed -e 's:</td>:</td>\n:g' | sed "s:href='\(.*\).html:href='sequence_overview/\1.html:g" >> $output # rewrite href to 'sequence_overview/..."
	echo "</div>" >> $output
	
	echo "</div>" >> $output #clonality tabber end
	
	echo "<br />" >> $output
	cat $dir/shm_clonality.htm >> $output
	
	echo "</div>" >> $output #clonality tab end

fi

# Use python's zipfile utility to prevent needing another dependency in the
# container.
current_dir=$(pwd)
cd $outdir
python -m zipfile -c all_outputs.zip \
  merged.txt filtered.txt unmatched.txt shm_overview.txt motif_per_seq.txt \
  mutation_by_id.txt base_overview.html aid_motives.txt relative_mutations.txt \
  absolute_mutations.txt tandems_by_id.txt scatter.txt frequency_ranges_class.txt \
  frequency_ranges_subclasses.txt transitions_all_sum.txt transitions_IGA_sum.txt \
  transitions_IGA1_sum.txt transitions_IGA2_sum.txt transitions_IGG_sum.txt \
  transitions_IGG1_sum.txt transitions_IGG2_sum.txt transitions_IGG3_sum.txt \
  transitions_IGG4_sum.txt transitions_IGM_sum.txt transitions_IGE_sum.txt \
  aa_id_mutations.txt absent_aa_id.txt aa_histogram_sum.txt \
  aa_histogram_sum_IGA.txt aa_histogram_sum_IGG.txt aa_histogram_sum_IGM.txt \
  aa_histogram_sum_IGE.txt baseline.txt baseline_IGA.pdf baseline_IGA.txt \
  baseline_IGG.pdf baseline_IGG.txt baseline_IGM.pdf baseline_IGM.txt \
  baseline_IGE.pdf baseline_IGE.txt IGA_pie.txt IGG_pie.txt \
  sequence_overview/index.html change_o/change-o-db-defined_clones*.txt \
  *.txz

cd $current_dir


echo "<div class='tabbertab' title='Downloads'>" >> $output

echo "<table class='pure-table pure-table-striped'>" >> $output
echo "<thead><tr><th>info</th><th>link</th></tr></thead>" >> $output
echo "<tr><td>All output files in a zip file</td><td><a href='all_outputs.zip' download='all_outputs.zip' >Download</a></td></tr>" >> $output
echo "<tr><td>The complete dataset</td><td><a href='merged.txt' download='merged.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The filtered dataset</td><td><a href='filtered.txt' download='filtered.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The alignment info on the unmatched sequences</td><td><a href='unmatched.txt' download='unmatched.txt' >Download</a></td></tr>" >> $output

echo "<tr><td colspan='2' style='background-color:#E0E0E0;'>SHM Overview</td></tr>" >> $output
echo "<tr><td>The SHM Overview table as a dataset</td><td><a href='shm_overview.txt' download='shm_overview.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>Motif data per sequence ID</td><td><a href='motif_per_seq.txt' download='motif_per_seq.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>Mutation data per sequence ID</td><td><a href='mutation_by_id.txt' download='mutation_by_id.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>Base count for every sequence</td><td><a href='base_overview.html'>View</a></td></tr>" >> $output
echo "<tr><td>The data used to generate the percentage of mutations in AID and pol eta motives plot</td><td><a href='aid_motives.txt' download='aid_motives.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data used to generate the relative mutation patterns plot</td><td><a href='relative_mutations.txt' download='relative_mutations.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data used to generate the absolute mutation patterns plot</td><td><a href='absolute_mutations.txt' download='absolute_mutations.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>Data about tandem mutations by ID</td><td><a href='tandems_by_id.txt' download='tandems_by_id.txt' >Download</a></td></tr>" >> $output

echo "<tr><td colspan='2' style='background-color:#E0E0E0;'>SHM Frequency</td></tr>" >> $output
echo "<tr><td>The data  generate the frequency scatter plot</td><td><a href='scatter.txt' download='scatter.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data used to generate the frequency by class plot</td><td><a href='frequency_ranges_classes.txt' download='frequency_ranges_classes.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for frequency by subclass</td><td><a href='frequency_ranges_subclasses.txt' download='frequency_ranges_subclasses.txt' >Download</a></td></tr>" >> $output

echo "<tr><td colspan='2' style='background-color:#E0E0E0;'>Transition Tables</td></tr>" >> $output
echo "<tr><td>The data for the 'all' transition plot</td><td><a href='transitions_all_sum.txt' download='transitions_all_sum.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for the 'IGA' transition plot</td><td><a href='transitions_IGA_sum.txt' download='transitions_IGA_sum.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for the 'IGA1' transition plot</td><td><a href='transitions_IGA1_sum.txt' download='transitions_IGA1_sum.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for the 'IGA2' transition plot</td><td><a href='transitions_IGA2_sum.txt' download='transitions_IGA2_sum.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for the 'IGG' transition plot</td><td><a href='transitions_IGG_sum.txt' download='transitions_IGG_sum.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for the 'IGG1' transition plot</td><td><a href='transitions_IGG1_sum.txt' download='transitions_IGG1_sum.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for the 'IGG2' transition plot</td><td><a href='transitions_IGG2_sum.txt' download='transitions_IGG2_sum.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for the 'IGG3' transition plot</td><td><a href='transitions_IGG3_sum.txt' download='transitions_IGG3_sum.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for the 'IGG4' transition plot</td><td><a href='transitions_IGG4_sum.txt' download='transitions_IGG4_sum.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for the 'IGM' transition plot</td><td><a href='transitions_IGM_sum.txt' download='transitions_IGM_sum.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for the 'IGE' transition plot</td><td><a href='transitions_IGE_sum.txt' download='transitions_IGE_sum.txt' >Download</a></td></tr>" >> $output

echo "<tr><td colspan='2' style='background-color:#E0E0E0;'>Antigen Selection</td></tr>" >> $output
echo "<tr><td>AA mutation data per sequence ID</td><td><a href='aa_id_mutations.txt' download='aa_id_mutations.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>Presence of AA per sequence ID</td><td><a href='absent_aa_id.txt' download='absent_aa_id.txt' >Download</a></td></tr>" >> $output

echo "<tr><td>The data used to generate the aa mutation frequency plot</td><td><a href='aa_histogram_sum.txt' download='aa_histogram_sum.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data used to generate the aa mutation frequency plot for IGA</td><td><a href='aa_histogram_sum_IGA.txt' download='aa_histogram_sum_IGA.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data used to generate the aa mutation frequency plot for IGG</td><td><a href='aa_histogram_sum_IGG.txt' download='aa_histogram_sum_IGG.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data used to generate the aa mutation frequency plot for IGM</td><td><a href='aa_histogram_sum_IGM.txt' download='aa_histogram_sum_IGM.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data used to generate the aa mutation frequency plot for IGE</td><td><a href='aa_histogram_sum_IGE.txt' download='aa_histogram_sum_IGE.txt' >Download</a></td></tr>" >> $output

echo "<tr><td>Baseline PDF (<a href='http://selection.med.yale.edu/baseline/'>http://selection.med.yale.edu/baseline/</a>)</td><td><a href='baseline.pdf' download='baseline.pdf' >Download</a></td></tr>" >> $output
echo "<tr><td>Baseline data</td><td><a href='baseline.txt' download='baseline.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>Baseline IGA PDF</td><td><a href='baseline_IGA.pdf' download='baseline_IGA.pdf' >Download</a></td></tr>" >> $output
echo "<tr><td>Baseline IGA data</td><td><a href='baseline_IGA.txt' download='baseline_IGA.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>Baseline IGG PDF</td><td><a href='baseline_IGG.pdf' download='baseline_IGG.pdf' >Download</a></td></tr>" >> $output
echo "<tr><td>Baseline IGG data</td><td><a href='baseline_IGG.txt' download='baseline_IGG.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>Baseline IGM PDF</td><td><a href='baseline_IGM.pdf' download='baseline_IGM.pdf' >Download</a></td></tr>" >> $output
echo "<tr><td>Baseline IGM data</td><td><a href='baseline_IGM.txt' download='baseline_IGM.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>Baseline IGE PDF</td><td><a href='baseline_IGE.pdf' download='baseline_IGE.pdf' >Download</a></td></tr>" >> $output
echo "<tr><td>Baseline IGE data</td><td><a href='baseline_IGE.txt' download='baseline_IGE.txt' >Download</a></td></tr>" >> $output

echo "<tr><td colspan='2' style='background-color:#E0E0E0;'>CSR</td></tr>" >> $output
echo "<tr><td>The data for the IGA subclass distribution plot</td><td><a href='IGA_pie.txt' download='IGA_pie.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The data for the IGG subclass distribution plot</td><td><a href='IGG_pie.txt' download='IGG_pie.txt' >Download</a></td></tr>" >> $output


echo "<tr><td colspan='2' style='background-color:#E0E0E0;'>Clonal Relation</td></tr>" >> $output
echo "<tr><td>Sequence overlap between subclasses</td><td><a href='sequence_overview/index.html'>View</a></td></tr>" >> $output
echo "<tr><td>The Change-O DB file with defined clones and subclass annotation</td><td><a href='change_o/change-o-db-defined_clones.txt' download='change_o/change-o-db-defined_clones.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The Change-O DB defined clones summary file</td><td><a href='change_o/change-o-defined_clones-summary.txt' download='change_o/change-o-defined_clones-summary.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just just the first sequence of a clone</td><td><a href='new_IMGT_first_seq_of_clone.txz' download='new_IMGT_first_seq_of_clone.txz' >Download</a></td></tr>" >> $output

echo "<tr><td>The Change-O DB file with defined clones of IGA</td><td><a href='change_o/change-o-db-defined_clones-IGA.txt' download='change_o/change-o-db-defined_clones-IGA.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The Change-O DB defined clones summary file of IGA</td><td><a href='change_o/change-o-defined_clones-summary-IGA.txt' download='change_o/change-o-defined_clones-summary-IGA.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just just the first sequence of a clone (IGA)</td><td><a href='new_IMGT_IGA_first_seq_of_clone.txz' download='new_IMGT_IGA_first_seq_of_clone.txz' >Download</a></td></tr>" >> $output

echo "<tr><td>The Change-O DB file with defined clones of IGG</td><td><a href='change_o/change-o-db-defined_clones-IGG.txt' download='change_o/change-o-db-defined_clones-IGG.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The Change-O DB defined clones summary file of IGG</td><td><a href='change_o/change-o-defined_clones-summary-IGG.txt' download='change_o/change-o-defined_clones-summary-IGG.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just just the first sequence of a clone (IGG)</td><td><a href='new_IMGT_IGG_first_seq_of_clone.txz' download='new_IMGT_IGG_first_seq_of_clone.txz' >Download</a></td></tr>" >> $output

echo "<tr><td>The Change-O DB file with defined clones of IGM</td><td><a href='change_o/change-o-db-defined_clones-IGM.txt' download='change_o/change-o-db-defined_clones-IGM.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The Change-O DB defined clones summary file of IGM</td><td><a href='change_o/change-o-defined_clones-summary-IGM.txt' download='change_o/change-o-defined_clones-summary-IGM.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just just the first sequence of a clone (IGM)</td><td><a href='new_IMGT_IGM_first_seq_of_clone.txz' download='new_IMGT_IGM_first_seq_of_clone.txz' >Download</a></td></tr>" >> $output

echo "<tr><td>The Change-O DB file with defined clones of IGE</td><td><a href='change_o/change-o-db-defined_clones-IGE.txt' download='change_o/change-o-db-defined_clones-IGE.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>The Change-O DB defined clones summary file of IGE</td><td><a href='change_o/change-o-defined_clones-summary-IGE.txt' download='change_o/change-o-defined_clones-summary-IGE.txt' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just just the first sequence of a clone (IGE)</td><td><a href='new_IMGT_IGE_first_seq_of_clone.txz' download='new_IMGT_IGE_first_seq_of_clone.txz' >Download</a></td></tr>" >> $output

echo "<tr><td colspan='2' style='background-color:#E0E0E0;'>Filtered IMGT output files</td></tr>" >> $output
echo "<tr><td>An IMGT archive with just the matched and filtered sequences</td><td><a href='new_IMGT.txz' download='new_IMGT.txz' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just the matched and filtered IGA sequences</td><td><a href='new_IMGT_IGA.txz' download='new_IMGT_IGA.txz' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just the matched and filtered IGA1 sequences</td><td><a href='new_IMGT_IGA1.txz' download='new_IMGT_IGA1.txz' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just the matched and filtered IGA2 sequences</td><td><a href='new_IMGT_IGA2.txz' download='new_IMGT_IGA2.txz' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just the matched and filtered IGG sequences</td><td><a href='new_IMGT_IGG.txz' download='new_IMGT_IGG.txz' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just the matched and filtered IGG1 sequences</td><td><a href='new_IMGT_IGG1.txz' download='new_IMGT_IGG1.txz' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just the matched and filtered IGG2 sequences</td><td><a href='new_IMGT_IGG2.txz' download='new_IMGT_IGG2.txz' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just the matched and filtered IGG3 sequences</td><td><a href='new_IMGT_IGG3.txz' download='new_IMGT_IGG3.txz' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just the matched and filtered IGG4 sequences</td><td><a href='new_IMGT_IGG4.txz' download='new_IMGT_IGG4.txz' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just the matched and filtered IGM sequences</td><td><a href='new_IMGT_IGM.txz' download='new_IMGT_IGM.txz' >Download</a></td></tr>" >> $output
echo "<tr><td>An IMGT archive with just the matched and filtered IGE sequences</td><td><a href='new_IMGT_IGE.txz' download='new_IMGT_IGE.txz' >Download</a></td></tr>" >> $output

echo "</table>" >> $output

echo "<br />" >> $output
cat $dir/shm_downloads.htm >> $output

echo "</div>" >> $output #downloads tab end

echo "</div>" >> $output #tabs end 

echo "</html>" >> $output


echo "---------------- naive_output.r ----------------"
echo "---------------- naive_output.r ----------------<br />" >> $log

if [[ "$naive_output" == "yes" ]]
then
	echo "output naive output"
	if [[ "${class_filter}" == "101_101" ]]
	then
		echo "copy new_IMGT.txz to ${naive_output_all}"
		cp $outdir/new_IMGT.txz ${naive_output_all}
	else
		echo "copy for classes"
		cp $outdir/new_IMGT_IGA.txz ${naive_output_ca}
		cp $outdir/new_IMGT_IGG.txz ${naive_output_cg}
		cp $outdir/new_IMGT_IGM.txz ${naive_output_cm}
		cp $outdir/new_IMGT_IGE.txz ${naive_output_ce}
	fi
fi

echo "</table>" >> $outdir/base_overview.html

mv $log $outdir/log.html

echo "<html><center><h1><a href='index.html'>Click here for the results</a></h1>Tip: Open it in a new tab (middle mouse button or right mouse button -> 'open in new tab' on the link above)<br />" > $log
echo "<table border = 1>" >> $log
echo "<thead><tr><th>Info</th><th>Sequences</th><th>Percentage</th></tr></thead>" >> $log
tIFS="$TMP"
IFS=$'\t'
while read step seq perc
	do
		echo "<tr>" >> $log
		echo "<td>$step</td>" >> $log
		echo "<td>$seq</td>" >> $log
		echo "<td>${perc}%</td>" >> $log
		echo "</tr>" >> $log
done < $outdir/filtering_steps.txt
echo "</table>" >> $log
echo "<br />" >> $log
cat $dir/shm_first.htm >> $log
echo "</center></html>" >> $log

IFS="$tIFS"

echo "---------------- remove_files----------------"
echo "---------------- remove_files----------------<br />" >> $log

rm -r -v -f $outdir/baseline
rm -r -v -f $PWD/files
rm -v $PWD/aa.txt
rm -v $PWD/aa_change_stats.txt
rm -v $PWD/gapped_aa.txt
rm -v $PWD/gapped_nt.txt
rm -v $PWD/hotspots.txt
rm -v $PWD/junction.txt
rm -v $PWD/mutationanalysis.txt
rm -v $PWD/mutationstats.txt
rm -v $PWD/sequences.txt
rm -v $PWD/summary.txt
rm -v $PWD/Rplots.pdf

filename="$dir/remove_files.txt"

while read file; do
    rm -v -f $outdir/$file
done < "$filename"

echo "---------------- Done! ----------------"
echo "---------------- Done! ----------------<br />" >> $outdir/log.html






