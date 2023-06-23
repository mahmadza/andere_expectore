
###############################################
#  call germline variants from WES normal
#     call on the whole WES regions
###############################################

cd ~/work/tumor_project/data/germline_variants

#download tumor-normal pair manifest
aws s3 cp s3://XXXXXX/pairs.txt .

#####################################
#      call germline variants
#####################################

screen -R germline_variants
awk '$0!~/^#/' pairs.txt | awk '!seen[$4]++' | while read l; do

  normal_plate=$(echo $l | awk '{print $3}')
  normal_prefix=$(echo $l | awk '{print $4}')
  sample=$(echo $normal_prefix | sed 's/_MN//')

  echo $normal_prefix $sample

  aws s3 cp s3://XXXXXX/Normal_Plate${normal_plate}/BQSR/${normal_prefix}.recalibrated.bam ${normal_prefix}.bam
  aws s3 cp s3://XXXXXX/Normal_Plate${normal_plate}/BQSR/${normal_prefix}.recalibrated.bai ${normal_prefix}.bai

  #call variants using HaplotypeCaller in GVCF mode
  java -jar /opt/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T HaplotypeCaller \
     -R /home/ubuntu/shared/ref/hg19/with_decoy_genome/hs37d5.fa \
     -I ${normal_prefix}.bam \
     --emitRefConfidence GVCF \
     --dbsnp /home/ubuntu/shared/ref/hg19/BroadInstitute/dbsnp_138.b37.snps.vcf.gz \
     -L /home/ubuntu/shared/kits/nexterarapidcapture_exome_targetedregions_v1.2.37.bed.gz -ip 100 \
     -o raw_variants/${sample}.raw.snps.indels.g.vcf.gz \
     2>raw_variants_log/${sample}.raw.snps.indels.g.vcf.gz.log

   #repair header, put sample name there, and recreate gvcf file
   zcat raw_variants/${sample}.raw.snps.indels.g.vcf.gz | \
     awk -vp=$sample -vOFS="\t" '($1=="#CHROM"){$NF=p}{print}' | bgzip -c -@ $(nproc) > $t
   mv $t raw_variants/${sample}.raw.snps.indels.g.vcf.gz

  aws s3 cp raw_variants/${sample}.raw.snps.indels.g.vcf.gz s3://XXXX/Normal_Plate${normal_plate}/HaplotypeCaller_nextera/ \
    --storage-class STANDARD_IA

  #remove files, and leave log file behind to verify all went OK
  rm ${normal_prefix}.ba* raw_variants/${sample}.raw.snps.indels.g.vcf.gz raw_variants/${sample}.raw.snps.indels.g.vcf.gz.tbi

done

#check logs
ls *raw.snps.indels.g.vcf.gz.log | while read l; do
  echo $l $(grep "Done. --" $l)
done
#all OK!


#######DONE
