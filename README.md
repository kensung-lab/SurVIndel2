# SurVIndel2

An deletion and tandem duplication caller for Illumina paired-end WGS data.

## Installation

In order to compile the code, the following are required:
- A C and a C++ compiler are required. If the GCC suite is used, version 4.9.3 or above is required.
- CMake (2.8 or above)

The following commands should be sufficient:

```
git clone https://github.com/kensung-lab/SurVIndel2
cd SurVIndel2/
./build_htslib.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

If you are compiling on the same platform as where you will execute it, you can use -DNATIVE=ON to create faster executables
```
cmake -DCMAKE_BUILD_TYPE=Release -DNATIVE=ON . && make
```

Python is necessary to run SurVIndel2. Libraries NumPy (http://www.numpy.org/), PyFaidx (https://github.com/mdshw5/pyfaidx) and PySam (https://github.com/pysam-developers/pysam) are required. If 
Python 2 is used, numpy 1.16.6, pyfaidx 0.5.9 and pysam 0.16.0.1 are the recommended (i.e., tested) versions. If Python 3 is used, then numpy 1.21.2, pyfaidx 0.5.9.1 and pysam 0.16.0.1 were 
tested.

## Building the machine learning model

scikit-learn (https://scikit-learn.org) must be installed. Download the file ml-training-data.zip from the latest release, and place it within the SurVIndel2 folder.
Then, you can run

```
unzip ml-training-data.zip
mkdir ml-model
python3 train_classifier.py ml-training-data/HG00096,ml-training-data/HG00171,ml-training-data/HG00512,ml-training-data/HG00513,ml-training-data/HG00514,ml-training-data/HG00731,ml-training-data/HG00732,ml-training-data/HG00733,ml-training-data/HG00864,ml-training-data/HG01114,ml-training-data/HG01505,ml-training-data/HG01596,ml-training-data/HG02011,ml-training-data/HG02492,ml-training-data/HG02587,ml-training-data/HG02818,ml-training-data/HG03009,ml-training-data/HG03065,ml-training-data/HG03125,ml-training-data/HG03371,ml-training-data/HG03486,ml-training-data/HG03683,ml-training-data/HG03732,ml-training-data/NA12329,ml-training-data/NA12878,ml-training-data/NA18534,ml-training-data/NA18939,ml-training-data/NA19238,ml-training-data/NA19239,ml-training-data/NA19240,ml-training-data/NA19650,ml-training-data/NA19983,ml-training-data/NA20509,ml-training-data/NA20847 ALL ml-model/
```

This will build the model in the a folder called ml-model. The process may take a while.

## Running

SurVIndel2 needs a BAM/CRAM file, a (possibly empty) working directory and reference genome in FASTA format.
The BAM/CRAM file must be coordinate-sorted and indexed. Furthermore, the MC and the MQ tag must be present for all primary alignments, when applicable.

Recent versions of BWA MEM (0.7.17) will add the MC tag. The easiest (but probably not the fastest) way to add the MQ tag is to use Picard FixMateInformation 
(http://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation) 
```
java -jar picard.jar FixMateInformation I=file.bam
```

The basic command to run SurVIndel2 is
```
python survindel2.py --threads N_THREADS BAM_FILE WORKDIR REFERENCE_FASTA
```

For other parameters, please see the help with
```
python survindel2.py -h
```

## Output

The output is a standard VCF file. It will be placed under WORKDIR/out.pass.vcf.gz. These are the deletions and duplications that SurVIndel2 deemed confident enough. 

The file WORKDIR/out.vcf.gz contains all of the deletions and duplcations, including those that did not pass the filters. Most of them will be false positives. It is not recommend to use this file unless for specific situations (e.g., you are looking for something specific).


If you built the machine learning model, you can use it to produce a more accurate set of calls:
```
python run_classifier.py WORKDIR/out.vcf WORKDIR/out.pass-ml.vcf.gz WORKDIR/stats.txt ALL ml-model/
```

This will generate a file WORKDIR/out.pass-ml.vcf.gz which contains the calls that the machine learning model predicted as real.
