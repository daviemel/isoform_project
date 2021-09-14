#!/bin/sh

datadir=$1 #directory of where all the data/fastqfiles are stored (absolute path)
outputdir=$2 #directory of where you want the output to be stored
genomePath=$3 # absolute path to processed genome index to be used for STAR alignment
annotationPath=$4 # abs path to genome annotation file (typically .gtf file)
filesofInterest=$5 # abs path to .csv containing .fastq file names

files_of_interest=$(while read -r line; do printf '%s\n' "$line"; done < "${filesofInterest}")

echo "Processing files from $datadir ."

#initial setup/prep
cd $outputdir
mkdir STAR_scripts
mkdir STAR_results
res_fold="${outputdir}STAR_results/"
mkdir $res_fold/coord_bams
mkdir $res_fold/pcrless_bams
mkdir $res_fold/star_logs
mkdir $res_fold/star_logs/PCRLess
mkdir $res_fold/star_logs/Init_QCMetrics
mkdir $res_fold/star_logs/Init_Runmsgs

cd $datadir

for file in $files_of_interest; do
    echo "Processing ${file}.fastq.tar"
    tar -xf ${file}.fastq.tar;
done


tempi=0

for file in $files_of_interest; do
    r1=${file}_R1.fastq.gz
    r2=${file}_R2.fastq.gz
    echo "$r1"
    echo "$r2"

    
    #Writing the script to process each sample
    echo '#!/bin/bash' >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo module load STAR >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mkdir $res_fold$file >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo STAR --genomeDir $genomePath --runThreadN 12 --sjdbGTFfile $annotationPath --readFilesIn "$datadir/$r1" "$datadir/$r2" --outFileNamePrefix "$res_fold$file/$file" --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mv "$res_fold$file/$file"Aligned.sortedByCoord.out.bam "${res_fold}coord_bams/" >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mv "$res_fold$file/$file"Log.final.out "${res_fold}star_logs/Init_QCMetrics/" >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    echo mv "$res_fold$file/$file"Log.out "${res_fold}star_logs/Init_Runmsgs/" >> "${outputdir}/STAR_scripts/STARParaScript${tempi}.sh"
    chmod +x $outputdir/STAR_scripts/STARParaScript$tempi.sh
    echo "A script has been generated for $file . Script has the name STARParaScript$tempi.sh and its results will start with $file ."
    ((tempi=tempi+1))
done

cd $outputdir/STAR_scripts/

for file in $(ls *.sh); do
  echo "${outputdir}STAR_scripts/$file" >> STARParaCom.txt
done

mv STARParaCom.txt $outputdir

cd $outputdir
