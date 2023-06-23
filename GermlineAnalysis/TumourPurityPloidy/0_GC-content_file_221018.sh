
# This bash script is to generate genome-wide GC content profile
# and subsequently upload the profile to S3

# Generate genome-wide GC content file
# Use hg19 with decoy as the genome
sequenza-utils gc_wiggle --fasta /home/ubuntu/shared/ref/hg19/with_decoy_genome/hs37d5.fa -w 50 | \
gzip > hs37d5.gc50Base.txt.gz

# Upload to S3
aws s3 cp hs37d5.gc50Base.txt.gz s3://XXXXX/hs37d5.gc50Base.txt.gz
