

###############################################
#   separate files of the germline variants
###############################################
#use AWS SQS to help with the queueing of the jobs
cd /home/mzabidi/temp/

#first, populate SQS queue
#use the URL as provided at the AWS Console
#populate the SQD queue with the samples to be processed
URL=https://sqs.ap-southeast-1.amazonaws.com/XXXXXXXXXXX
for sample in $(awk 'NR>1' ~/tumor_project/data/germline_variants/pairs.txt | cut -f4 | sort | uniq | sed 's/_MN//'); do
  echo $sample
  aws sqs send-message --queue-url $URL --message-body $sample
done

#spawn workers
#do in multiple Screens
#leave ~20-30% free memory and CPU
screen -R split_by_SQS_screen_1

SQS=$(mktemp)
URL=https://sqs.ap-southeast-1.amazonaws.com/XXXXXXXXXXX
#poll queue to begin
aws sqs receive-message --queue-url $URL > $SQS

#keep polling until the queue exhausts
#also check if no other worker working on the same sample
#use the existence of log file to check
while [ -s $SQS ] && [ ! -f "log/${sample}.log" ]; do
  #parse message
  #grab sample to be processed
  #and the resource handle
  sample=$(cat $SQS | awk '$1=="\"Body\":"{print $2}' | sed 's/\"//g' | sed 's/,//')
  rhandle=$(cat $SQS | awk '$1=="\"ReceiptHandle\":"{print substr($2,2,length($2)-3)}')

  echo $sample $rhandle

  #split samples from the big VCF file
  time /usr/lib/jvm/java-8-oracle/bin/java -jar /opt/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R /home/ubuntu/shared/ref/hg19/with_decoy_genome/hs37d5.fa \
  -V MyBrCa.recal_snps.recal_indels.vcf.gz \
  -sn $sample -o ${sample}.vcf.gz \
  2>log/${sample}.log

  #delete message
  #use the resource resource handle
  #observation: actually the message will still be deleted if the resource handle is not an exact match (?)
  aws sqs delete-message --queue-url $URL --receipt-handle $rhandle

  #poll next message to keep the while loop running
  aws sqs receive-message --queue-url $URL > $SQS

done

#check logs see if all has succeeded
ls log/SD*log | while read l; do
  grep "Done. ---" $l
done
#all done! great!



#######DONE
