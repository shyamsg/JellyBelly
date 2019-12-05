# JellyBelly
Toolset for vectorizing sequence data via spaced kmers

# JellyBelly is a new project that will be in constant change. This master branch will always contain a working version. If you want to follow the latest changes and new functionality, be sure to add the "dev" branch. The dev branch will probably be a huge mess.

# Graphic summary

![graphsum](/misc/graphsum.svg)

# Dependencies
  JellyBelly has no mayor dependencies other than zlib.
  
  JellyBelly uses [klib](https://github.com/attractivechaos/klib) (included in this repo) for handling sequence data and hashing of kmer sequences.
  
  JellyBelly uses cmake and make for compilation.

# Get and compile
    git clone --recursive https://github.com/7PintsOfCherryGarcia/JellyBelly.git
    cd JellyBelly/build
    cmake ..
    make


# Run

     Usage:

	 JellyBelly [options] -f <sequence file> -s <spacedkmer file>

     Options:

    -f <filename>	Sequence file. fasta/q file with sequences or
	  		sequencing reads. Can be gziped. provide "-" if
	  		reading from stdin. If using sequencing data eg.
	  		illumina reads, provide the -C flag to use canonical
	  		form.

    -s <filename>	Spaced kmer file. Binary file containing the masks
	  		to be used. Refer to manual for detailed documentation.

    -C 		Canonical mode. Lexicographically smallest kmer is counted.
	  		Set this flag when analyzing sequencing reads.

    -h 		        This help  message.

# Usage

## Input Data
  JellyBelly can be executed on any type of DNA nucleotide data in the FASTA/FASTQ format. There is only **one** thing to consider. If you are runing JellyBelly on **sequencing reads**, make shure to set the **-C** flag. This will have a two fold efect:
  
  1. The entire sequencing library will be encoded in a single spaced kmer vector.
  1. Only lexicographically smallest kmers will be computed. Reads coming from random fragments can be from either the forward or reverse strand.
  
  For example:
  
  Assuming the sequencing library "my_seq.fq.gz", and the spaced kmer file "SpacedKmer_K10_S5.bin" you can run JellyBelly in the following way:
    
    JellyBelly -f my_seq.fq.gz -s SpacedKmer_K10_S5.bin -C > output.txt
  
  output.txt will be a tab separated text file with a count between 0 and 1 inclusive for each spaced kmer. The number of counts (number of spaced kmers) depends on the length of the number of ones in the spaced kmer mask and is equal to 4^S where S is the the number of ones in the mask.
  
  JellyBelly's output values are scaled spaced kmer counts ranging from 0 to 1 includisve. I will add another option for raw output (actual counts). Spaced kmer count values are sorted lexicographcally eg. first value corresponds to AAA..A second to AAA..T and so on.
  
  output.txt
  
      1	  2	  3	  	4^S
    
    0.245	0.014	0.547	...    0.436
  
  If you are running JellyBelly on assembled contigs or reference genomes. Every sequence in a fasta file will be encoded into a spaced kmer vector. If you want a set of sequences to be encoded into a single vector you will have to concatenate them. I will soon add another option for computing a single vector from a set of sequences.
  
  


# spaced kmer file and how to compute them
  A core concept of JellyBelly is the use of a "mask" to only consider specific positions of a kmer. For example given the following kmer of length 10:
  
                                                   ACGGTGCAAT
						   
  and the following mask:
  
                                                   1001100001
						   
  The corresponding spaced kmer will be:
  
                                                   A  GT    T == AGTT

Brieafly, given a kmer size K and a spaced kmer size S, there are S choose K different ways of selecting S positions out of a kmer of length K. JellyBelly includes another program to compute these masks called create_spaces. create_spaces will do a brute force search of all possible mask configurations. Subsequently, only the most entropic(dissorderly) masks are kept.

To create a spaced kmer file simply run create_spaces with your desired kmer size and spaced kmer file.

    create_spaces <kmer length> <smer length>
	
This will create a .bin file. This is a binary file that contains the entropic masks. This file should be used as the -s option for JellyBelly.

JellyBelly is hardcoded to use the first mask in the .bin file as of now. Options will be added to make this costumizable.

# Utilities

python/\*.py

Series of helper python scripts for ploting PCA and UMAP. You will need numpy, sklearn, plotly, and colorlover. You can install these with:

pip install numpy sklearn plotly colorlover

These scripts plot html file that you can open and interact with in a web browser.

build/distMatrix

This is a small C program that takes spaced kmer vectors and computes all pairwise euclidean distances.

# Upcoming changes

* Spaced kmer files
    * Currewnt version of create_spaces keeps only the spaced kmer masks with maximum entropy. Future versions will keep all spaced kmer masks that a within a certain fraction of the most entropic one.
    * JellyBelly uses only the first spaced kmer mask. Future versions will let the user choose which mask to use. This will require another utility program to list all the masks in a spaced kmer file.
* Output of JellyBelly
    * Currently JellyBelly outputs spaced kmer scaled counts in text format. Future versions will output this information in binary format to a filename specified as an option. Converting binary output to text output will be implemented separately.
* JellyBelly's execution
    * Currently JellyBelly has only one way of executing it. Future versions will introduce 'commands'. For example, there will be a comand to vectorize sequence data, another comand to print binary vector files, and another comand scale raw counts.
