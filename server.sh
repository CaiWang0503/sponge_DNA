##################################log in##################################
ssh beswcai@genome.ljmu.ac.uk #150.204.78.6
#password
FuzzyL0g1c!
ssh genome2 #go to genome2
screen #go into virtual environment
#screen -r  # detached
#screen -r 88
#control+a+d # detached
#control+l #move command line to top

#rm -rf 4_L001-ds.3401d0eec9344c9a99065dea0e58da1c.fastq #remove and restart from beginning #do not run this

###################cope and install#####################################
#make a directory for reference sequences
mkdir db_obitools
mv EMBL_r143* ./db_obitools # move all the files into this folder
#mv to genome2
scp -r db_obitools beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/tank/
#move data to genome2, stay in genome1 cope files to genome2
scp -r 3_L001-ds.815457c2482c4a36a6bfa50baeb62673 beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/tank/
scp -r 4_L001-ds.3401d0eec9344c9a99065dea0e58da1c beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/tank/
#upload ngsfilter files from local
scp -r /Users/wang/Desktop/ob1_mbc/sponge_tank_L3_ngsfilter.txt beswcai@genome.ljmu.ac.uk:/home/beswcai/
scp -r /Users/wang/Desktop/ob1_mbc/sponge_tank_L4_ngsfilter.txt beswcai@genome.ljmu.ac.uk:/home/beswcai/
#log in and mv files to genome2
scp -r sponge* beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/tank/ #stay in genome1 cope files to genome2
#go to genome2 and make a folder for objective1
mkdir tank
#scp -r * beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/tank/ #stay in genome1 cope files to genome2
#scp -r * beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/applications/
#put data into "tank" folder
mv /home/genome2/beswcai/tank/3_L001-ds.815457c2482c4a36a6bfa50baeb62673/SpongeDNA-obj1-LIB01_S3_L001_R1_001.fastq   /home/genome2/beswcai/tank
mv /home/genome2/beswcai/tank/3_L001-ds.815457c2482c4a36a6bfa50baeb62673/SpongeDNA-obj1-LIB01_S3_L001_R2_001.fastq  /home/genome2/beswcai/tank
mv /home/genome2/beswcai/tank/4_L001-ds.3401d0eec9344c9a99065dea0e58da1c/SpongeDNA-obj2-LIB02_S4_L001_R1_001.fastq   /home/genome2/beswcai/tank
mv /home/genome2/beswcai/tank/4_L001-ds.3401d0eec9344c9a99065dea0e58da1c/SpongeDNA-obj2-LIB02_S4_L001_R2_001.fastq  /home/genome2/beswcai/tank

###################bioinformatics start here#####################################
##the original obitols tutoria can find here: https://pythonhosted.org/OBITools/wolves.html
##the original peter shum tutoria can be find here: https://github.com/shump2/cobble2012/blob/master/obi_cobble.job
###########################################################################################
#setup all the PATH here, and then can use directly
vsearch=~/applications/vsearch-2.18.0/bin/vsearch #usage: $vsearch
swarm=~/applications/swarm/src/swarm #usage: $swarm
export PATH=$PATH:~/applications/seqtk/
#############1. quality control raw data############
mkdir fastqc
fastqc -o fastqc/ --extract -f fastq *.fastq

############2. trimming & Recover full sequence reads from forward and reverse partial reads######
# enter OBITools virtualenv
source /home/genome2/beswcai/applications/OBITools-venv/bin/activate

#loop obicut trimming
find * -maxdepth 0 -name "*_001.fastq" > samplelist.txt  # find all files ending with _1.fq.gz
sample_info=samplelist.txt # put samplelist.txt into variable
sample_names=($(cut -f1,2,3,4 -d"_" "${sample_info}" | uniq)) # convert variable to array
    echo "${sample_names[@]}" # echo all array elements
    echo ${sample_names[1]} # echo first array element
echo "There are" "${#sample_names[@]}" "samples that will be processed" # echo number of elements in the array

for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
      sample_prefix="$( basename $sample "_001.fastq")"
      echo ${sample_prefix}
      obicut -e 150 ${sample_prefix}_001.fastq > ${sample_prefix}_trim150.fastq
done

#trimming separated
#obicut -e 150 SpongeDNA-obj1-LIB01_S3_L001_R1_001.fastq > SpongeDNA-obj1-LIB01_S3_L001_R1_trim150.fastq
#obicut -e 150 SpongeDNA-obj1-LIB01_S3_L001_R2_001.fastq > SpongeDNA-obj1-LIB01_S3_L001_R2_trim150.fastq
#obicut -e 150 SpongeDNA-obj2-LIB02_S4_L001_R1_001.fastq > SpongeDNA-obj2-LIB02_S4_L001_R1_trim150.fastq
#obicut -e 150 SpongeDNA-obj2-LIB02_S4_L001_R2_001.fastq > SpongeDNA-obj2-LIB02_S4_L001_R2_trim150.fastq

#check files size
ls -lht SpongeDNA-obj1-LIB01_S3_L001_R2_trim150.fastq
#echo Paired-end alignment. Annotate the reads with quality 40 and split the output in two files
#illuminapairedend -r cobb12_trim180.R2.fastq cobb12_trim250.R1.fastq | obiannotate -S goodali:'"Good_cobbCOI" if score>40.00 else "Bad_cobbCOI"' | obisplit -t goodali
illuminapairedend -r SpongeDNA-obj1-LIB01_S3_L001_R1_trim150.fastq SpongeDNA-obj1-LIB01_S3_L001_R2_trim150.fastq | obiannotate -S goodali:'"Good_sp12s.s3" if score>40.00 else "Bad_sp12s.s3"' | obisplit -t goodali
illuminapairedend -r SpongeDNA-obj2-LIB02_S4_L001_R1_trim150.fastq SpongeDNA-obj2-LIB02_S4_L001_R2_trim150.fastq | obiannotate -S goodali:'"Good_sp12s.s4" if score>40.00 else "Bad_sp12s.s4"' | obisplit -t goodali
ls -lht Good_sp12s.s4.fastq #check files size, make sure it's no empty

#the following command only use good sequences
###########3. convert fastq to fasta for demultiplexing in parallel######
seqtk seq -a Good_sp12s.s3.fastq > Good_sp12s_s3.fasta
seqtk seq -a Good_sp12s.s4.fastq > Good_sp12s_s4.fasta

#clean folder before moving on
mkdir intermediate
mv *.fastq ./intermediate
#rm -rf 4_L001-ds.3401d0eec9344c9a99065dea0e58da1c
#rm -rf 3_L001-ds.815457c2482c4a36a6bfa50baeb62673

###########4. Assign each sequence record to the corresponding sample/marker combination######
mkdir demulti
ngsfilter -t sponge_tank_L3_ngsfilter.txt --fasta-output -u unidentified_sp12s_s3.fasta Good_sp12s_s3.fasta --DEBUG > sp12stank_s3.filtered.fasta
ngsfilter -t sponge_tank_L4_ngsfilter.txt --fasta-output -u unidentified_sp12s_s4.fasta Good_sp12s_s4.fasta --DEBUG > sp12stank_s4.filtered.fasta
ls -lht sp12stank_s3.filtered.fasta
mv sp12stank* ./demulti
mv unidentified* ./demulti
#concatenate all the separate *.filtered.fasta into a single file
ngsfilter_results=~/tank/demulti
cat $(find $ngsfilter_results -name '*.filtered.fasta' | xargs)> demulti/spongetank.filtered.fasta
cat $(find $ngsfilter_results -name '*unidentified*' | xargs)> demulti/unidentified.spongetank.fasta
cd demulti
ls -lht spongetank.filtered.fasta # make sure it concatenates successful

###########5. Filter seqs##########
#echo Filter the seqs with length between 140 and 220 bp and with no 'N' #-p 'count>=10'
#obigrep -p 'seq_length>140' -p 'seq_length<220' -s '^[ACGT]+$' demulti/spongetank.filtered.fasta > spongetank.filtered_length.fasta
obigrep -p 'seq_length>140' -p 'seq_length<220' spongetank.filtered.fasta > spongetank.filtered_length.fasta
#'^[ACGT]+$' do not work! why?
#obigrep -s '^[ACGT]+$' spongetank.filtered_length.fasta > spongetank.filtered_length_noN.fasta
#ls -lht spongetank.filtered_length_noN.fasta

###########6. Get the count statistics##########
#echo Calculate stats per sample
obistat -c sample -a seq_length spongetank.filtered_length.fasta > sample_stats_spongetank.length_filter.txt

###########7. Dereplicate reads into uniq sequences##########
#echo Group the unique seqs
obiuniq -m sample spongetank.filtered_length.fasta > spongetank.unique.fasta
head -5 spongetank.unique.fasta

###########8. Exchange the identifier to a short index##########
#if i use "%10d" here, swarm will show error
obiannotate --seq-rank spongetank.unique.fasta | obiannotate --set-identifier '"'tank'_%0*d" %(9,seq_rank)' > spongetank.new9.fasta
head -5 spongetank.new9.fasta

###########9. convert to vsearch format#####
#Rscript ~/applications/R_scripts_metabarpark/owi_obifasta2vsearch -i spongetank.new10.fasta -o spongetank.vsearch.fasta
Rscript /Users/wang/Desktop/ob1_mbc/R_scripts_metabarpark/owi_obifasta2vsearch -i spongetank.new9.fasta -o spongetank.vsearch.fasta
head -5  spongetank.vsearch.fasta
sed 's/ ;/;/g' spongetank.vsearch.fasta > spongetank.vsearch.mod.fasta
#head -5  spongetank.vsearch.mod.fasta
wc -l spongetank.vsearch.mod.fasta

###########10. CHIMERA DETECTION ##########
#echo Run UCHIME de novo in VSEARCH
mkdir vsearch_output
#vsearch --uchime_denovo spongetank.vsearch.mod.fasta --sizeout --nonchimeras vsearch_output/spongetank.nonchimeras.fasta --chimeras vsearch_output/spongetank.chimeras.fasta --threads 28 --uchimeout vsearch_output/spongetank.uchimeout2.txt &> vsearch_output/log.spongetank_chimeras
#sed 's/;/ ;/g' vsearch_output/spongetank.nonchimeras.fasta |grep -e ">" | awk 'sub(/^>/, "")' | awk '{print $1}' > vsearch_output/spongetank.nonchimeras.txt # text file used for owi_recount_sumaclust step
wc -l ./vsearch_output/spongetank.nonchimeras.fasta
head -5 ./vsearch_output/spongetank.nonchimeras.fasta

###########11. CLUSTERING ##########
#echo swarm using vsearch nonchimeras file
mkdir swarm_output
$swarm -d 3 -z -t 10 -o swarm_output/spongetank_SWARM3_output -s swarm_output/spongetank_SWARM3_stats -w swarm_output/spongetank_SWARM3_seeds.fasta vsearch_output/spongetank.nonchimeras.fasta
#wc -l ./swarm_output/spongetank_SWARM3_seeds.fasta
head -5 ./swarm_output/spongetank_SWARM3_seeds.fasta

################################
##### TAXONOMIC ASSIGNMENT #####
################################
scp -r beswcai@genome.ljmu.ac.uk:/home/beswcai/db_obitools/  /Users/wang/Desktop/
#Use ecoPCR to simulate an in silico` PCR for teleo2 primer
ecoPCR -d ./taxo_peter20211011/EMBL_r143 -e 3 -l 140 -L 220 AAACTCGTGCCAGCCACC GGGTATCTAATCCCAGTTTG > tele02.ecopcr
mkdir
#Clean the database
obigrep -d ./taxo_peter20211011/EMBL_r143 --require-rank=species --require-rank=genus --require-rank=family tele02.ecopcr > tele02_clean.fasta
obiuniq -d ./taxo_peter20211011/EMBL_r143 tele02_clean.fasta > tele02_clean_uniq.fasta
obigrep -d ./taxo_peter20211011/EMBL_r143 --require-rank=family tele02_clean_uniq.fasta > tele02_clean_uniq_clean.fasta
obiannotate --uniq-id tele02_clean_uniq_clean.fasta > db_tele02.fasta

####################################################################
ecotag -d ./taxo_peter20211011/EMBL_r143 -R db_tele02.fasta --sort=count -r ./demulti/swarm_output/spongetank_SWARM3_seeds.fasta > spongetank_SWARM3.ecotag.fasta
sed 's/;s/; s/g' spongetank_SWARM3.ecotag.fasta > spongetank_SWARM3.ecotag_NEW.fasta

#mkdir ecotag
#cd ecotag
#echo here use the script "submit_parallel_ecotag.sh" to parallelize the ecotag command, you need to edit the paths to you sumaclust generated
#echo fasta file, e.g. cobble2012.sumaclust95.centers.fasta, and ecopcr database. This script generates files and folders by splitting the
#echo cobble2012.sumaclust95.centers.fasta file into 100 sequences per file and 100 files per folder. It will generate  as many files/folders.
#echo You can set the split files into hatever you like and it will run each file as a separate job. Typically 100 sequneces per file run in about ~10-15 minutes.
#sh submit_parallel_ecotag.sh
#echo Once the previous step has complete, e.g. overnight, concatenate all *.ecotag.fasta from all folders
#ecotag_results=~/stanford/Cobble_final/ecotag_all/
#cat $(find $ecotag_results -name '*.ecotag.fasta' | xargs)> ecotag_all/cobble2012.ecotag.fasta

#echo sort ecotag.fasta
#grep ">" spongetank_SWARM3.ecotag.fasta | sed 's/>//g' | sort -k1.6n > spongetank.ecotag_idlist.txt
#cdbfasta spongetank_SWARM3.ecotag.fasta -o spongetank.ecotag.fasta.index
#cat spongetank.ecotag_idlist.txt | cdbyank spongetank.ecotag.fasta.index > spongetank_SWARM3.ecotag_sorted.fasta
#rm cobble2012.ecotag.fasta.index
#rm spongetank.ecotag_idlist.txt

######################
##R scripts for reformatting metabarcoding databases CREDIT: OWEN WANGENSTEEN Find R scripts here: https://github.com/metabarpark/R_scripts_metabarpark
#echo Add taxa above order level

#Rscript ~/applications/R_scripts_metabarpark/owi_add_taxonomy ecotag_all/cobble2012.ecotag_sorted.fasta cobble2012.ecotag.fasta.annotated.csv
#change  dir_taxo <- "/Users/wang/Desktop/ob1_mbc/R_scripts_metabarpark/dir_taxo/"
Rscript ./R_scripts_metabarpark/owi_add_taxonomy spongetank_SWARM3.ecotag_NEW.fasta spongetank_SWARM3.ecotag.annotated.csv
sed 's/;";/";/g' spongetank_SWARM3.ecotag.annotated.csv > spongetank_SWARM3.ecotag.annotated_id.csv
rm spongetank_SWARM3.ecotag.annotated.csv
#echo recount abundance by sample
obitab -o ./demulti/spongetank.new9.fasta > spongetank.new.tab
#Rscript ~/applications/R_scripts_metabarpark/owi_recount_swarm ./demulti/swarm_output/spongetank_SWARM3_output spongetank.new.tab
Rscript ./R_scripts_metabarpark/owi_recount_swarm ./demulti/swarm_output/spongetank_SWARM3_output spongetank.new.tab
mv ./demulti/swarm_output/spongetank_SWARM3_output.counts.csv  spongetank_SWARM3_output.counts.csv

#echo combine ecotag and abundance files
#Rscript ~/peter/applications/R_scripts_metabarpark/owi_combine -i cobble2012.ecotag.fasta.annotated.csv -a swarm/cobble2012_SWARM13_output.counts.csv -o cobble2012_all_SWARM_FINAL_MOTUs.csv
Rscript ./R_scripts_metabarpark/owi_combine -i spongetank_SWARM3.ecotag.annotated_id.csv -a spongetank_SWARM3_output.counts.csv -o spongetank_all_SWARM_FINAL_MOTUs.csv
sed 's/;/,/g' spongetank_all_SWARM_FINAL_MOTUs.csv > spongetank_all_SWARM_FINAL_OTU.csv

#echo collapse MOTUs
Rscript ./R_scripts_metabarpark/owi_collapse -s 16 -e 116 -t 0.50 -i spongetank_all_SWARM_FINAL_OTU.csv
#-s 14 Sample columns start; -e sample columns end. Default = 98; -t 0.50 Threshold for collapsing
