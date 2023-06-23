



#requires a lot of memory

#probably better to only focus on cellularity and ploidy for now
#and write into a file

cd ~/work/tumor_project/data/Sequenza

#download pairs.txt file
aws s3 cp s3://XXXXX/somatic_variants/pairs.txt .


#download seq files
grep -v "#" pairs.txt | while read l; do

  tumor_plate=$(echo $l | awk '{print $1}')
  tumor_prefix=$(echo $l | awk '{print $2}')
  normal_plate=$(echo $l | awk '{print $3}')
  normal_prefix=$(echo $l | awk '{print $4}')
  sample=$(echo $tumor_prefix | sed 's/_TUM//')

  echo $tumor_prefix $normal_prefix $sample

  #aws s3 cp s3://XXXXX/purity/Sequenza/${tumor_prefix}.seqz.gz .
  #download only small file as it's already enough
  aws s3 cp s3://XXXXX/purity/Sequenza/small_files/${tumor_prefix}.small.seqz.gz .

done

#skip SD0014 for now
ls *small*gz | sed 's/.small.seqz.gz//' | grep -v SD0014 | awk 'NR>=1 && NR<=100' > samples_small_1_to_100.txt
ls *small*gz | sed 's/.small.seqz.gz//' | grep -v SD0014 | awk 'NR>=101 && NR<=200' > samples_small_101_to_200.txt
ls *small*gz | sed 's/.small.seqz.gz//' | grep -v SD0014 | awk 'NR>=201 && NR<=300' > samples_small_201_to_300.txt
ls *small*gz | sed 's/.small.seqz.gz//' | grep -v SD0014 | awk 'NR>=301 && NR<=400' > samples_small_301_to_400.txt
ls *small*gz | sed 's/.small.seqz.gz//' | grep -v SD0014 | awk 'NR>=401 && NR<=500' > samples_small_401_to_500.txt
ls *small*gz | sed 's/.small.seqz.gz//' | grep -v SD0014 | awk 'NR>=501 && NR<=600' > samples_small_501_to_600.txt


mkdir -p cellularity ploidy



screen -R process_sequenza_501_600

R

library("sequenza")

for(sample in c(read.table("samples_small_501_to_600.txt"))$V1)
{

  print(sample)

  #extract sequenza data
  seqz_data=sequenza.extract(paste0(sample,".small.seqz.gz"),chromosome.list=c(seq(1:22),"X"))

  #calculate cellularity and ploidy
  #result as list(ploidy, collularity, lpp)
  seqz_fit=sequenza.fit(seqz_data)

  #get confidence intervals of the parameters
  seqz_conf_int<-get.ci(seqz_fit)

  #output cellularity (at the point of the maximum posterior probability)
  file_cell_out=file(paste0("cellularity/",sample,"_cell.txt"))
  writeLines(as.character(seqz_conf_int$max.cellularity),file_cell_out)
  close(file_cell_out)

  #get ploidy
  file_ploidy_out=file(paste0("ploidy/",sample,"_ploidy.txt"))
  writeLines(as.character(seqz_conf_int$max.ploidy),file_ploidy_out)
  close(file_ploidy_out)

}



#'"
grep -v "#" pairs.txt | while read l; do
  tumor=$(echo $l | awk '{print $2}')
  ls ploidy/${tumor}_ploidy.txt
done
#all OK

grep -v "#" pairs.txt | while read l; do
  tumor=$(echo $l | awk '{print $2}')
  ls cellularity/${tumor}_cell.txt
done
#all OK


#create a summary file and re-upload to S3
grep -v "#" pairs.txt | while read l; do
  tumor=$(echo $l | awk '{print $2}')
  echo -en $tumor"\t"$(cat ploidy/${tumor}_ploidy.txt)"\n"
done > ploidy_all.txt

grep -v "#" pairs.txt | while read l; do
  tumor=$(echo $l | awk '{print $2}')
  echo -en $tumor"\t"$(cat cellularity/${tumor}_cell.txt)"\n"
done > cellularity_all.txt


aws s3 cp ploidy_all.txt s3://XXXXX/purity/Sequenza/small_files/
aws s3 cp cellularity_all.txt s3://XXXXX/purity/Sequenza/small_files/













#########
