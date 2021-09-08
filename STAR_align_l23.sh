#!/bin/sh

datadir=$1 #directory of where all the data/fastqfiles are stored (absolute path)
outputdir=$2 #directory of where you want the output to be stored
genomePath=$3 # absolute path to processed genome index to be used for STAR alignment
annotationPath=$4 # abs path to genome annotation file (typically .gtf file)

files_of_interest="LS-15025_S29_E1-50
SM-D9CY7_S59_E1-50
LS-15059_S69_E1-50
LS-15379_S04_E1-50
LS-15059_S60_E1-50
LS-15003_S28_E1-50
LS-15002_S29_E1-50
LS-15086_S68_E1-50
LS-15921_S95_E1-50
LS-15082_S71_E1-50
LS-15086_S69_E1-50
LS-15082_S42_E1-50
SM-D9CYH_S07_E1-50
SM-D9CY7_S11_E1-50
LS-15381_S66_E1-50
SM-D9E62_S71_E1-50
LS-15031_S06_E1-50
SM-D9E4J_S76_E1-50
LS-14696_S40_E1-50
SM-D9CY7_S28_E1-50
LS-15379_S08_E1-50
SM-D9CYB_S16_E1-50
SM-D9E5J_S55_E1-50
LS-15086_S67_E1-50
LS-15043_S60_E1-50
LS-15029_S16_E1-50
SM-D9E5J_S26_E1-50
LS-15006_S78_E1-50
LS-15921_S67_E1-50
SM-D9CY7_S16_E1-50
LS-15921_S25_E1-50
SM-D9CZK_S28_E1-50
LS-15379_S31_E1-50
SM-D9CYH_S03_E1-50
LS-15006_S52_E1-50
LS-15060_S44_E1-50
LS-15081_S62_E1-50
LS-15038_S41_E1-50
LS-15921_S39_E1-50
LS-15381_S74_E1-50
LS-15921_S80_E1-50
LS-15933_S40_E1-50
LS-15006_S51_E1-50
LS-15020_S89_E1-50
LS-15933_S70_E1-50
SM-D9E63_S83_E1-50
LS-15380_S34_E1-50
LS-15347_S24_E1-50
LS-15002_S28_E1-50
LS-15083_S53_E1-50
LS-15080_S81_E1-50
SM-D9CZF_S41_E1-50
LS-15349_S21_E1-50
LS-15346_S66_E1-50
LS-15381_S58_E1-50
SM-D9CY7_S17_E1-50
SM-D9CYA_S68_E1-50
LS-15024_S76_E1-50
LS-15006_S82_E1-50
LS-15060_S28_E1-50
LS-15346_S68_E1-50
LS-15921_S96_E1-50
LS-15068_S77_E1-50
SM-D9CZT_S30_E1-50
LS-15921_S56_E1-50
SM-DAIFV_S94_E1-50
LS-15921_S22_E1-50
LS-15921_S79_E1-50
LS-15349_S13_E1-50
LS-14696_S46_E1-50
SM-D9E4M_S90_E1-50
LS-15080_S86_E1-50
LS-14696_S42_E1-50
LS-15921_S72_E1-50
LS-15921_S82_E1-50
LS-15024_S88_E1-50
LS-15376_S14_E1-50
LS-15080_S80_E1-50
LS-15052_S42_E1-50
LS-15381_S80_E1-50
LS-15381_S71_E1-50
SM-D9E61_S88_E1-50
LS-15347_S25_E1-50
LS-15912_S52_E1-50
LS-15381_S62_E1-50
LS-15380_S24_E1-50
LS-15083_S63_E1-50
LS-15380_S08_E1-50
SM-D9CY7_S26_E1-50
LS-15921_S47_E1-50
SM-D9EQV_S06_E1-50
LS-15921_S26_E1-50
LS-15067_S81_E1-50
LS-15379_S09_E1-50
LS-15060_S40_E1-50
LS-15059_S70_E1-50
LS-15038_S29_E1-50
LS-15044_S10_E1-50
LS-15082_S52_E1-50
LS-15381_S94_E1-50"

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
