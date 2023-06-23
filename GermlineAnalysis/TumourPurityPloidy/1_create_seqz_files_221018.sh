




cd ~/work/tumor_project/data/Sequenza

mkdir -p seqz.gz.logs/ small.seqz.gz.logs/

#grab tumor and normal sample names and their plates

#download tumor-normal pair manifest
aws s3 cp s3://XXXXX/somatic_variants/pairs.txt .

#download GC content file
aws s3 cp s3://XXXXX/hs37d5.gc50Base.txt.gz .


screen -R Sequenza_last
#here need to check what are the latest samples that are being worked on
awk '$1!~/^#/' pairs.txt | while read l; do

  tumor_plate=$(echo $l | awk '{print $1}')
  tumor_prefix=$(echo $l | awk '{print $2}')
  normal_plate=$(echo $l | awk '{print $3}')
  normal_prefix=$(echo $l | awk '{print $4}')
  sample=$(echo $tumor_prefix | sed 's/_TUM//')

  echo $tumor_prefix $normal_prefix $sample

  #download files
  aws s3 cp s3://XXXXX/Tumor_Plate${tumor_plate}/BQSR/${tumor_prefix}.recalibrated.bam ${tumor_prefix}.bam
  aws s3 cp s3://XXXXX/WES/Tumor_Plate${tumor_plate}/BQSR/${tumor_prefix}.recalibrated.bai ${tumor_prefix}.bai
  aws s3 cp s3://XXXXX/WES/Normal_Plate${normal_plate}/BQSR/${normal_prefix}.recalibrated.bam ${normal_prefix}.bam
  aws s3 cp s3://XXXXX/WES/Normal_Plate${normal_plate}/BQSR/${normal_prefix}.recalibrated.bai ${normal_prefix}.bai

  #create seqz file
  time sequenza-utils bam2seqz --tumor ${tumor_prefix}.bam --normal ${normal_prefix}.bam \
    -gc hs37d5.gc50Base.txt.gz --fasta /home/ubuntu/shared/ref/hg19/with_decoy_genome/hs37d5.fa \
      2>seqz.gz.logs/${tumor_prefix}.seqz.gz.log | gzip > ${tumor_prefix}.seqz.gz
  #takes 3.5hours/sample

  #also do binned seqz file
  time sequenza-utils seqz_binning --window 50 --seqz ${tumor_prefix}.seqz.gz | \
    gzip > ${tumor_prefix}.small.seqz.gz 2>small.seqz.gz.logs/${tumor_prefix}.small.seqz.gz.log

  aws s3 cp ${tumor_prefix}.seqz.gz s3://XXXXX/purity/Sequenza/${tumor_prefix}.seqz.gz \
  --storage-class STANDARD_IA
  aws s3 cp ${tumor_prefix}.small.seqz.gz s3://XXXXX/WES/purity/Sequenza/small_files/${tumor_prefix}.small.seqz.gz \
  --storage-class STANDARD_IA

  rm ${tumor_prefix}.ba* ${normal_prefix}.ba* ${tumor_prefix}.seqz.gz ${tumor_prefix}.small.seqz.gz
  #leave the log file behind

done





#check logs
cd seqz.gz.logs
ls *log | while read l; do
  echo $l
  grep "[mpileup]" $l
  echo "-------------"
done









#########
