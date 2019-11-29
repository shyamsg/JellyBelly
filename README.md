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

    -C 		        Canonical mode. Lexicographically smallest kmer is counted.
	  		Set this flag when analyzing sequencing reads.

   -h 		        This help  message.
