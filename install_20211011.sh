#########installing history#########
#seqtk cdbfasta obitools ecoPCR VSEARCH
mkdir applications
ls ~/applications/

#install fastqc in local mac
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip
chmod 755 fastqc
/Users/wang/src/FastQC/fastqc -help

#install cdbfasta
#https://github.com/gpertea/cdbfasta
#Clone this repository and run 'make' in the cdbfasta directory.
#Make sure you have the zlib library and development files installed on your linux system.
#Running 'make' should produce the binaries 'cdbfasta' (the indexer program)
git clone https://github.com/gpertea/cdbfasta.git
cd cdbfasta
#~/src/cdbfasta-master
make

#install seqtk
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
/home/genome2/beswcai/applications/seqtk/seqtk -h

# install OBITools
virtualenv -python=python2 OBITools-virtualenv
source OBITools/bin/activate

#install ecoPCR
git clone https://git.metabarcoding.org/obitools/ecopcr.git
cd ecopcr/src
make
sudo cp ecoPCR /usr/local/bin/.
sudo cp ecofind /usr/local/bin/.
sudo cp ecogrep /usr/local/bin/.

#install ecoprimers
git clone https://git.metabarcoding.org/obitools/ecoprimers.git
cd ecoprimers/src
make

# Install VSEARCH
#VSEARCH is a set of tools for working with metabarcoding data
wget https://github.com/torognes/vsearch/archive/v2.18.0.tar.gz
tar xzf v2.18.0.tar.gz
cd vsearch-2.18.0
./autogen.sh
./configure CFLAGS="-O3" CXXFLAGS="-O3"
make

#Install swarm
git clone https://github.com/torognes/swarm.git
cd swarm/src/
make

###########path################
export PATH=$PATH:~/applications/seqtk/
