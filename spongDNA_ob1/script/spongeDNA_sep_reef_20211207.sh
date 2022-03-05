###################bioinformatics start here#####################################
##the original obitols tutorial can find here: https://pythonhosted.org/OBITools/wolves.html
##the original reference tutorial we used can be find here: https://github.com/shump2/cobble2012/blob/master/obi_cobble.job
###########################################################################################
#setup all the PATH here, and then can use directly
#home:/Users/wang/src/sponge_DNA/spongDNA_ob1/sep_reef

#trimming separated
obicut -e 145 SpongeDNA-obj1-LIB01_S3_L001_R1_001.fastq > SpongeDNA-obj1-LIB01_S3_L001_R1_trim145.fastq
obicut -e 145 SpongeDNA-obj1-LIB01_S3_L001_R2_001.fastq > SpongeDNA-obj1-LIB01_S3_L001_R2_trim145.fastq
obicut -e 145 SpongeDNA-obj2-LIB02_S4_L001_R1_001.fastq > SpongeDNA-obj2-LIB02_S4_L001_R1_trim145.fastq
obicut -e 145 SpongeDNA-obj2-LIB02_S4_L001_R2_001.fastq > SpongeDNA-obj2-LIB02_S4_L001_R2_trim145.fastq

#echo Paired-end alignment. Annotate the reads with quality 40 and split the output in two files
illuminapairedend -r SpongeDNA-obj1-LIB01_S3_L001_R1_trim145.fastq SpongeDNA-obj1-LIB01_S3_L001_R2_trim145.fastq | obiannotate -S goodali:'"Good_sp12s.s3" if score>40.00 else "Bad_sp12s.s3"' | obisplit -t goodali
wc -l Good_sp12s.s3.fastq #check files size, make sure it's no empty
illuminapairedend -r SpongeDNA-obj2-LIB02_S4_L001_R1_trim145.fastq SpongeDNA-obj2-LIB02_S4_L001_R2_trim145.fastq | obiannotate -S goodali:'"Good_sp12s.s4" if score>40.00 else "Bad_sp12s.s4"' | obisplit -t goodali
wc -l Good_sp12s.s4.fastq #check files size, make sure it's no empty

#the following command only use good sequences
###########3. convert fastq to fasta for demultiplexing in parallel######
seqtk seq -a Good_sp12s.s3.fastq > Good_sp12s_s3.fasta
seqtk seq -a Good_sp12s.s4.fastq > Good_sp12s_s4.fasta

#clean folder before moving on
mkdir intermediate
mv *.fastq ./intermediate

###########4. Assign each sequence record to the corresponding sample/marker combination######
#in sep_reef folder i remove the reef samples in the "sponge_tank_L3/L4_ngsfilter.txt"
mkdir demulti
cd demulti
ngsfilter -t /Users/wang/src/sponge_DNA/spongDNA_ob1/sep_reef/sponge_tank_L3_ngsfilter.txt --fasta-output -u unidentified_sp12s_s3.fasta /Users/wang/src/sponge_DNA/spongDNA_ob1/sep_reef/Good_sp12s_s3.fasta --DEBUG > sp12stank_s3.filtered.fasta
ngsfilter -t /Users/wang/src/sponge_DNA/spongDNA_ob1/sep_reef/sponge_tank_L4_ngsfilter.txt --fasta-output -u unidentified_sp12s_s4.fasta /Users/wang/src/sponge_DNA/spongDNA_ob1/sep_reef/Good_sp12s_s4.fasta --DEBUG > sp12stank_s4.filtered.fasta
wc -l sp12stank_s3.filtered.fasta #544590
wc -l unidentified_sp12s_s3.fasta #620905
wc -l sp12stank_s4.filtered.fasta #2267708
wc -l unidentified_sp12s_s4.fasta #115625

#concatenate all the separate *.filtered.fasta into a single file
ngsfilter_results=/Users/wang/src/sponge_DNA/spongDNA_ob1/sep_reef/demulti
cat $(find $ngsfilter_results -name '*.filtered.fasta' | xargs)> spongetank.filtered.fasta
cat $(find $ngsfilter_results -name '*unidentified*' | xargs)> unidentified.spongetank.fasta
#cd demulti
seqkit --threads 3 stats $ngsfilter_results/spongetank.filtered.fasta > mergedtrimmed_read_counts_filtered.txt #718,030
seqkit --threads 3 stats $ngsfilter_results/unidentified.spongetank.fasta > mergedtrimmed_read_counts_unidentified.txt#71,883

###########5. Filter seqs##########
#echo Filter the seqs with length between 140 and 220 bp and with no 'N' #-p 'count>=10'
#'^[ACGT]+$' do not work! if i paste directly, i have to type by hand, because atom do not recognize this is a python expression.
obigrep -p 'seq_length>130' -p 'seq_length<190' -s '^[ACGT]+$' spongetank.filtered.fasta > spongetank.filtered_length_noN.fasta

#ls -lht spongetank.filtered_length_noN.fasta
seqkit --threads 3 stats spongetank.filtered_length_noN.fasta > spongetank.filtered_length_noN.txt #660,201
###########6. Get the count statistics##########
#echo Calculate stats per sample
obistat -c sample -a seq_length spongetank.filtered_length_noN.fasta > sample_stats_spongetank.length_filter.txt
###########7. Dereplicate reads into uniq sequences##########
#echo Group the unique seqs
obiuniq -m sample spongetank.filtered_length_noN.fasta > spongetank.unique.fasta
head -5 spongetank.unique.fasta
obistat -c count spongetank.unique.fasta | sort -nk1 | head -50
seqkit --threads 3 stats spongetank.unique.fasta > spongetank.unique.txt #80,829
###########8. Exchange the identifier to a short index##########
obiannotate --seq-rank spongetank.unique.fasta | obiannotate --set-identifier '"'tank'_%0*d" %(9,seq_rank)' > spongetank.new9.fasta
head -2 spongetank.new9.fasta
###########9. convert to vsearch format#####
#i annotate the install command in the file  "owi_obifasta2vsearch", because the install code is different from the old R version.
Rscriptpath=/Users/wang/src/sponge_DNA/spongDNA_ob1/R_scripts_metabarpark
Rscript $Rscriptpath/owi_obifasta2vsearch -i spongetank.new9.fasta -o spongetank.vsearch.fasta
sed 's/ ;/;/g' spongetank.vsearch.fasta > spongetank.vsearch.mod.fasta
seqkit --threads 3 stats spongetank.vsearch.mod.fasta > spongetank.vsearch.mod.txt #

###########10. CHIMERA DETECTION ##########
#echo Run UCHIME de novo in VSEARCH
mkdir vsearch_output
vsearch --uchime_denovo spongetank.vsearch.mod.fasta --sizeout --nonchimeras vsearch_output/spongetank.nonchimeras.fasta --chimeras vsearch_output/spongetank.chimeras.fasta --threads 28 --uchimeout vsearch_output/spongetank.uchimeout.txt &> vsearch_output/log.spongetank_chimeras
seqkit --threads 3 stats ./vsearch_output/spongetank.nonchimeras.fasta > spongetank.nonchimeras.txt #79,646
###########11. CLUSTERING ##########
#echo swarm using vsearch nonchimeras file
mkdir swarm_output
# -d 1 # default and recommended# # -f # fastidious# # -z # use usearch size= for abundance
swarm -d 3 -z -t 20 -o swarm_output/spongetank_SWARM3_output -s swarm_output/spongetank_SWARM3_stats -w swarm_output/spongetank_SWARM3_seeds.fasta vsearch_output/spongetank.nonchimeras.fasta
################################
##### TAXONOMIC ASSIGNMENT #####
################################
#scp -r beswcai@genome.ljmu.ac.uk:/home/beswcai/db_obitools/  /Users/wang/Desktop/
##Use ecoPCR to simulate an in silico` PCR for teleo2 primer
#ecoPCR -d ./taxo_peter20211011/EMBL_r143 -e 3 -l 130 -L 190 AAACTCGTGCCAGCCACC GGGTATCTAATCCCAGTTTG > tele02.ecopcr
#mkdir
##Clean the database
#obigrep -d ./taxo_peter20211011/EMBL_r143 --require-rank=species --require-rank=genus --require-rank=family tele02.ecopcr > tele02_clean.fasta
#obiuniq -d ./taxo_peter20211011/EMBL_r143 tele02_clean.fasta > tele02_clean_uniq.fasta
#obigrep -d ./taxo_peter20211011/EMBL_r143 --require-rank=family tele02_clean_uniq.fasta > tele02_clean_uniq_clean.fasta
#obiannotate --uniq-id tele02_clean_uniq_clean.fasta > db_tele02.fasta

#############ecotag run on genome######################
#######################################################
#
scp -r /Users/wang/src/sponge_DNA/spongDNA_ob1/sep_reef/demulti/swarm_output/spongetank_SWARM3_seeds.fasta beswcai@genome.ljmu.ac.uk:/home/beswcai/
#log in and mv files to genome2
scp -r spongetank* beswcai@genome2.ljmu.ac.uk:/home/genome2/beswcai/tank/ #stay in genome1 cope files to genome2
#run code peter has a parallel script to boost running (for large dataset), please find the command on the original script.
source /home/genome2/beswcai/applications/OBITools-venv/bin/activate
ecotag -d ./db_obitools/taxo_peter20211011/EMBL_r143 -R db_tele02_embl143_20211011.fasta --sort=count -r spongetank_SWARM3_seeds.fasta > spongetank_SWARM3.ecotag.sepreef.20211207.fasta
#download files
scp -r spongetank_SWARM3.ecotag.sepreef* beswcai@genome.ljmu.ac.uk:/home/beswcai/ #stay in genome2 cope files to genome1
scp -r beswcai@genome.ljmu.ac.uk:/home/beswcai/spongetank_SWARM3.ecotag.sepreef* /Users/wang/Desktop/ #stay in local
home=/Users/wang/src/sponge_DNA/spongDNA_ob1/sep_reef
cd $home
sed 's/;s/; s/g' spongetank_SWARM3.ecotag.sepreef.20211207.fasta > spongetank_SWARM3.ecotag.sepreef.new.20211207.fasta
######################
##R scripts for reformatting metabarcoding databases CREDIT: OWEN WANGENSTEEN Find R scripts here: https://github.com/metabarpark/R_scripts_metabarpark
###important!!######change  dir_taxo <- "/Users/wang/Desktop/ob1_mbc/R_scripts_metabarpark/dir_taxo/"
Rscript $Rscriptpath/owi_add_taxonomy spongetank_SWARM3.ecotag.sepreef.new.20211207.fasta spongetank_SWARM3.ecotag.sepreef.annotated.csv
sed 's/;";/";/g' spongetank_SWARM3.ecotag.sepreef.annotated.csv > spongetank_SWARM3.ecotag.sepreef.annotated_id.csv
rm spongetank_SWARM3.ecotag.sepreef.annotated.csv

#echo recount abundance by sample
obitab -o ./demulti/spongetank.new9.fasta > spongetank.new.tab
#Rscript ~/applications/R_scripts_metabarpark/owi_recount_swarm ./demulti/swarm_output/spongetank_SWARM3_output spongetank.new.tab
Rscript $Rscriptpath/owi_recount_swarm ./demulti/swarm_output/spongetank_SWARM3_output spongetank.new.tab

#echo combine ecotag and abundance files
#Rscript ~/peter/applications/R_scripts_metabarpark/owi_combine -i cobble2012.ecotag.fasta.annotated.csv -a swarm/cobble2012_SWARM13_output.counts.csv -o cobble2012_all_SWARM_FINAL_MOTUs.csv
Rscript $Rscriptpath/owi_combine -i spongetank_SWARM3.ecotag.sepreef.annotated_id.csv -a ./demulti/swarm_output/spongetank_SWARM3_output.counts.csv -o spongetank_SWARM3_FINAL_MOTUs.csv
sed 's/;/,/g' spongetank_SWARM3_FINAL_MOTUs.csv > spongetank_SWARM3_FINAL_OTUs.csv

rm spongetank_SWARM3_FINAL_MOTUs.csv

# i delete "cut" column manually for "spongetank_all_SWARM_FINAL_OTU.csv"

#clean up
mkdir intermediate
mv *.fastq ./intermediate
mv *.fasta ./intermediate
mv *_id.csv ./intermediate
mv *.tab ./intermediate

#echo collapse MOTUs
#-s 14 Sample columns start; -e sample columns end. Default = 98; -t 0.50 Threshold for collapsing
#Rscript $Rscriptpath/owi_collapse -s 16 -e 116 -t 0.50 -i spongetank_10reads_SWARM3_OTU_LULU_20220220.csv

Rscript $Rscriptpath/owi_collapse -s 16 -e 104 -t 0.50 -i spongetank_SWARM3_OTU_LULU_20220220.csv

#local blastn
cat Tetrapoda_UK_species_12sgenbank_nogap_edit.fas fish_harper_refdata_20190726.fa > Uk.species_local.database.fas
mkdir local_blast
mv Gramma_loreto.fasta local_blast/Gramma_loreto.fasta
cd local_blast
makeblastdb -in Gramma_loreto.fasta -dbtype nucl # make the uk ref dataset BLASTABLE
blastn -db /Users/wang/src/sponge_DNA/tank_experiment_20211015/local_blast/Gramma_loreto.fasta -query /Users/wang/src/sponge_DNA/tank_experiment_20211015/spongetank_all_SWARM3_FINAL_OTU.fasta -num_threads 3  -max_target_seqs 1 -outfmt 6 -out sponge_query_Gramma_loreto20211019.txt
