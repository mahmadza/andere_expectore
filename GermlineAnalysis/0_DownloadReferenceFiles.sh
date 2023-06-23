

###############################################
#      download reference files for VQSR
###############################################

cd ~/work/tumor_project/data/germline_variants/resources

#connect via FTP to the Broad Institute
ftp ftp.broadinstitute.org
Name (ftp.broadinstitute.org:mzabidi): gsapubftp-anonymous
#turn on passive mode, otherwise will get 500 Illegal PORT command
ftp> pass

#set local FTP file download directory
lcd /home/MYUSER/work/tumor_project/data/germline_variants/resources

#download files in b37 folder
cd bundle/b37
get hapmap_3.3.b37.vcf.gz
get hapmap_3.3.b37.vcf.gz.md5
get 1000G_omni2.5.b37.vcf.gz
get 1000G_omni2.5.b37.vcf.gz.md5
get 1000G_phase1.indels.b37.vcf.gz
get 1000G_phase1.indels.b37.vcf.gz.md5
get 1000G_phase1.snps.high_confidence.b37.vcf.gz
get 1000G_phase1.snps.high_confidence.b37.vcf.gz.md5
get dbsnp_138.b37.vcf.gz
get dbsnp_138.b37.vcf.gz.md5
get Mills_and_1000G_gold_standard.indels.b37.vcf.gz
get Mills_and_1000G_gold_standard.indels.b37.vcf.gz.md5

#verify checksum
for x in $(ls *md5 | sed 's/.md5//'); do
  checksum=$(md5sum $x | awk '{print $1}')
  toverify=$(awk '{print $1}' $x.md5)
  if [ "$checksum" == "$toverify" ]; then
    echo $x PASSED
  fi
done

#ALL CLEAR!
#remove md5 files
rm *.md5

#expand from gunzip, bgzip and then create proper index
ls *.gz | while read l; do
  echo $l
  gunzip $l
  bgzip ${l%.gz}
  tabix $l &
done
