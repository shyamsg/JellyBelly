# JellyBelly
Toolset for vectorizing sequence data via spaced kmers

# If you are here because of my talk at the embl PhD symposium I promise to push JellyBelly's code and manual very soon. :)

# Dependencies
  JellyBelly has no mayor dependencies other than zlib.
  
  JellyBelly uses [klib](https://github.com/attractivechaos/klib) (included in this repo) for handling sequence data.
  
  JellyBelly uses cmake and make for compilation

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
  JellyBelly can be executed on any type of DNA nucleotide data. There is only one thing to consider. If you are runing JellyBelly on sequencing reads, make shure to set the -C flag. This will have a two fold efect:
  
  1. The entire sequencing library will be encoded in a single spaced kmer vector.
  2. Only lexicographically smallest kmers will be computed. Reads coming from random fragments can be from either the forward or reverse strand.
  
  If you are running JellyBelly on assembled contigs or reference genomes. Every sequence in a fasta file will be encoded into a spaced kmer vector. If you want a set of sequences to be encoded into a single vector you will have to concatenate them. I will soon add another option for computing a single vector from a set of sequences.
  
  JellyBelly will output spaced kmer vectors to stdout. Each value in the output vector is a scaled spaced kmer count. I will add another option for raw output. Kmer count values are sorted lexicographcally eg. first value corresponds to AAA..A second to AAA..T and so on.


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
