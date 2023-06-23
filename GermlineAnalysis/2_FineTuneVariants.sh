
#############################
#   consolidate gVCFs
#############################

#redownload gVCF files from S3
cd ~/work/tumor_project/data/germline_variants
mkdir -p raw_variants
awk '$0!~/^#/' pairs.txt | awk '!seen[$4]++' | while read l; do

  normal_plate=$(echo $l | awk '{print $3}')
  normal_prefix=$(echo $l | awk '{print $4}')
  sample=$(echo $normal_prefix | sed 's/_MN//')

  aws s3 cp s3://XXXX${normal_plate}/HaplotypeCaller_nextera/${sample}.raw.snps.indels.g.vcf.gz \
    raw_variants/
  #no need to download index files
  #just create new one
  tabix raw_variants/${sample}.raw.snps.indels.g.vcf.gz &

done

#combine cohort
java -jar /opt/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T CombineGVCFs \
  -R /home/ubuntu/shared/ref/hg19/with_decoy_genome/hs37d5.fa \
  $(awk '$0!~/^#/' pairs.txt | awk '!seen[$4]++{print $4}' | sed 's/_MN//' | \
    awk '{print "--variant raw_variants/"$1".raw.snps.indels.g.vcf.gz"}' | tr '\n' ' ') \
  -o MyBrCa.g.vcf.gz

#joint-call cohort
java -jar /opt/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt $(nproc) \
  -R /home/ubuntu/shared/ref/hg19/with_decoy_genome/hs37d5.fa \
  --variant MyBrCa.g.vcf.gz \
  -o MyBrCa.vcf.gz

#validate variants
#see if there's any errors
java -jar /opt/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T ValidateVariants \
  -R /home/ubuntu/shared/ref/hg19/with_decoy_genome/hs37d5.fa \
  --variant MyBrCa.vcf.gz \
  --dbsnp /home/ubuntu/shared/ref/hg19/BroadInstitute/dbsnp_138.b37.snps.vcf.gz \
=#all OK


#can remove the combined uncalled variants
rm MyBrCa.g.vcf.gz MyBrCa.g.vcf.gz.tbi


#############################
#       VQSR, SNVs
#############################

#dont invoke DP when working with exome datasets
#since there is extreme variation in the depth to which targets are captured
#In whole genome experiments this variation is indicative of error but that is not the case in capture experiments

#build SNP recalibration model
java -jar /opt/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar \
    -T VariantRecalibrator -nt $(nproc) \
    -R /home/ubuntu/shared/ref/hg19/with_decoy_genome/hs37d5.fa \
    -input MyBrCa.vcf.gz \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
      ~/work/tumor_project/data/germline_variants/resources/hapmap_3.3.b37.vcf.gz \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 \
      ~/work/tumor_project/data/germline_variants/resources/1000G_omni2.5.b37.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 \
      ~/work/tumor_project/data/germline_variants/resources/1000G_phase1.snps.high_confidence.b37.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
      ~/work/tumor_project/data/germline_variants/resources/dbsnp_138.b37.vcf.gz \
    -an QD \
    -an FS \
    -an SOR \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
    -an InbreedingCoeff \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -rscriptFile recalibrate_SNP_plots.R \
    2>recalibrate_SNP.recal.log

#apply the calibration to the SNPs
java -jar /opt/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar \
    -T ApplyRecalibration  -nt $(nproc) \
    -R /home/ubuntu/shared/ref/hg19/with_decoy_genome/hs37d5.fa \
    -input MyBrCa.vcf.gz \
    -mode SNP \
    --ts_filter_level 99.0 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -o MyBrCa.recal_snps.raw_indels.vcf.gz
#--ts_filter_level 99.0
#means will allow us to retrieve 99% of truth training set SNPs from HapMap and Omni
#best practice: take 99.9 instead actually

#compress recalibrate_SNP.recal
gzip recalibrate_SNP.recal
#compress recalibrate_SNP_plots.R
gzip recalibrate_SNP_plots.R


#############################
#       VQSR, indels
#############################

#build indel recalibration model
#CAUTION:Coverage (DP)
#Total (unfiltered) depth of coverage.
#Note that this statistic should not be used with exome datasets; see caveat detailed in the VQSR arguments FAQ doc.
java -jar /opt/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar \
    -T VariantRecalibrator -nt $(nproc) \
    -R /home/ubuntu/shared/ref/hg19/with_decoy_genome/hs37d5.fa \
    -input MyBrCa.recal_snps.raw_indels.vcf.gz \
    -resource:mills,known=false,training=true,truth=true,prior=12.0 \
      ~/work/tumor_project/data/germline_variants/resources/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
      ~/work/tumor_project/data/germline_variants/resources/dbsnp_138.b37.vcf.gz \
    -an QD \
    -an FS \
    -an SOR \
    -an MQRankSum \
    -an ReadPosRankSum \
    -an InbreedingCoeff \
    -mode INDEL \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    --maxGaussians 4 \
    -recalFile recalibrate_INDEL.recal \
    -tranchesFile recalibrate_INDEL.tranches \
    -rscriptFile recalibrate_INDEL_plots.R

#apply model
java -jar /opt/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar \
    -T ApplyRecalibration -nt $(nproc) \
    -R /home/ubuntu/shared/ref/hg19/with_decoy_genome/hs37d5.fa \
    -input MyBrCa.recal_snps.raw_indels.vcf.gz \
    -mode INDEL \
    --ts_filter_level 99.0 \
    -recalFile recalibrate_INDEL.recal \
    -tranchesFile recalibrate_INDEL.tranches \
    -rscriptFile recalibrate_SNP_plots.R \
    -o MyBrCa.recal_snps.recal_indels.vcf.gz
#best practice: take 99.9 actually


#remove intermediate files
rm MyBrCa.recal_snps.raw_indels.vcf.gz MyBrCa.recal_snps.raw_indels.vcf.gz.tbi



#######DONE
