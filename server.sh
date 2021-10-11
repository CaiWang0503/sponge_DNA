##########log in#########
ssh  beswcai@genome.ljmu.ac.uk #150.204.78.6
#password
FuzzyL0g1c!
ssh genome2 #go to genome2
screen #go into virtual environment
#screen -r  # detached
#screen -r 88
#control+a+d # exit
#control+l #move to top

#rm -rf 4_L001-ds.3401d0eec9344c9a99065dea0e58da1c.fastq #remove and restart from beginning
###################cope and install#####################################
#make a directory for reference sequences
mkdir db_obitools
mv EMBL_r143* ./db_obitools
#
scp -r 3_L001-ds.815457c2482c4a36a6bfa50baeb62673 beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/tank/ #stay in genome1 cope files to genome2
scp -r 4_L001-ds.3401d0eec9344c9a99065dea0e58da1c beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/tank/ #stay in genome1 cope files to genome2
#upload ngsfilter files
scp -r /Users/wang/Desktop/ob1_mbc/sponge_tank_L3_ngsfilter.txt beswcai@genome.ljmu.ac.uk:/home/beswcai/
scp -r /Users/wang/Desktop/ob1_mbc/sponge_tank_L4_ngsfilter.txt beswcai@genome.ljmu.ac.uk:/home/beswcai/
#log in and mv files to genome2
scp -r sponge* beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/tank/ #stay in genome1 cope files to genome2

mkdir tank
#scp -r * beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/tank/ #stay in genome1 cope files to genome2
#scp -r * beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/applications/

mv /home/genome2/beswcai/tank/3_L001-ds.815457c2482c4a36a6bfa50baeb62673/SpongeDNA-obj1-LIB01_S3_L001_R1_001.fastq   /home/genome2/beswcai/tank
mv /home/genome2/beswcai/tank/3_L001-ds.815457c2482c4a36a6bfa50baeb62673/SpongeDNA-obj1-LIB01_S3_L001_R2_001.fastq  /home/genome2/beswcai/tank
mv /home/genome2/beswcai/tank/4_L001-ds.3401d0eec9344c9a99065dea0e58da1c/SpongeDNA-obj2-LIB02_S4_L001_R1_001.fastq   /home/genome2/beswcai/tank
mv /home/genome2/beswcai/tank/4_L001-ds.3401d0eec9344c9a99065dea0e58da1c/SpongeDNA-obj2-LIB02_S4_L001_R2_001.fastq  /home/genome2/beswcai/tank
###################bioinformatics start here#####################################
##the original obitols tutoria can find here: https://pythonhosted.org/OBITools/wolves.html
##https://github.com/shump2/cobble2012/blob/master/obi_cobble.job
#############1. quality control obitools input data############
mkdir fastqc
fastqc -o fastqc/ --extract -f fastq *.fastq

############2. Recover full sequence reads from forward and reverse partial reads######
# enter virtualenv
source /home/genome2/beswcai/applications/OBITools-venv/bin/activate

#echo Paired-end alignment. Annotate the reads with quality 40 and split the output in two files
illuminapairedend -r SpongeDNA-obj1-LIB01_S3_L001_R1_001.fastq SpongeDNA-obj1-LIB01_S3_L001_R2_001.fastq | obiannotate -S goodali:'"Good_sp12s_s3" if score>40.00 else "Bad_sp12_s3"' | obisplit -t goodali_s3
illuminapairedend -r SpongeDNA-obj2-LIB02_S4_L001_R1_001.fastq SpongeDNA-obj2-LIB02_S4_L001_R2_001.fastq | obiannotate -S goodali:'"Good_sp12s_s4" if score>40.00 else "Bad_sp12s_s4"' | obisplit -t goodali_s4

#the following command only use good sequences
###########3. convert fastq to fasta for demultiplexing in parallel######
seqtk seq -a Good_sp12s_s3.fastq > Good_sp12s_s3.fasta
seqtk seq -a Good_sp12s_s4.fastq > Good_sp12s_s4.fasta

###########4. Assign each sequence record to the corresponding sample/marker combination######
mkdir demulti
ngsfilter -t sponge_tank_L3_ngsfilter.txt --fasta-output -u unidentified_sp12s.fastq Good_sp12s.fasta --DEBUG > sp12stank.filtered.fasta
mv sp12stank.filtered.fasta ./demulti/sp12stank.filtered.fasta

###########7. Filter the seqs with length between 220 and 260 bp and with no 'N'##########
#echo Filter the seqs with length between 140 and 220 bp and with no 'N' # -p 'count>=10'
obigrep -p 'seq_length>140' -p 'seq_length<220' -s '^[ACGT]+$' demulti/sp12stank.filtered_sorted.fasta > sp12stank.filtered_length.fasta

###########8. Get the count statistics##########
#echo Calculate stats per sample
obistat -c sample -a seq_length sp12stank.filtered_length.fasta > sample_stats_sp12stank.length_filter.txt |  \sort -nk1 | head -20

###########9. Dereplicate reads into uniq sequences##########
#echo Group the unique seqs
obiuniq -m sample sp12stank.filtered_length.fasta > sp12stank.unique.fasta






#####################
##### CLUSTERING ####
#####################

#echo swarm using vsearch nonchimeras file
#mkdir swarm
#~/peter/applications/swarm/src/swarm -d 3 -z -t 10 -o swarm/cobble2012_SWARM13_output -s swarm/cobble2012_SWARM13_stats -w swarm/cobble2012_SWARM13_seeds.fasta vsearch/cobble2012.nonchimeras.fasta
