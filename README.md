# JellyBelly
Toolset for vectorizing sequence data via spaced kmers

# If you are here because of my talk at the embl PhD symposium I promise to push JellyBelly's code and manual very soon. :)

# Dependencies
  JellyBelly has no mayor dependencies other than zlib.
  
  JellyBelly uses [klib](https://github.com/attractivechaos/klib) (included in this repo) for handling sequence data.
  
  JellyBelly uses cmake and make for compilation

# get and compile
  git clone https://github.com/7PintsOfCherryGarcia/JellyBelly.git
  
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

    -C 		       Canonical mode. Lexicographically smallest kmer is counted.
	  		Set this flag when analyzing sequencing reads.

    -h 		        This help  message.

# Usage
  JellyBelly can be executed on any type of DNA nucleotide data. There is only one thing to consider. If you are runing JellyBelly on sequencing reads, make shure to set the -C flag. This will have a two fold efect:
  
  1. The entire sequencing library will be encoded in a single spaced kmer vector.
  2. Only lexicographically smallest kmers will be computed. Reads coming from random fragments can be from either the forward or reverse strand.
  
  If you are running JellyBelly on assembled contigs or reference genomes. Every sequence in a fasta file will be encoded into a spaced kmer vector. If you want a set of sequences to be encoded into a single vector you will have to concatenate them. I will soon add another option for computing a single vector from a set of sequences.
  
  JellyBelly will output spaced kmer vectors to stdout. Each value in the output vector is a scaled spaced kmer count. I will add another option for raw output. Kmer count values are sorted lexicographcally eg. first value corresponds to AAA..A second to AAA..T and so on.

