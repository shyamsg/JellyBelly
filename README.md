# JellyBelly

## Toolset for vectorizing sequence data.

# JellyBelly is a new project that will be in constant change. This master branch will always contain a working version. If you want to follow the latest changes and new functionality, be sure to add the "dev" branch. The dev branch will probably be a huge mess.

# Graphic summary

![graphsum](/misc/graphsum.svg)

# What is JellyBelly?

  JellyBelly is a tool to encode sequence data as vectors of constant length (See graphic summary). The main concept of JellyBelly is the **spaced kmer**. A spaced kmer is a sequence of length S derived from a sequence of length K with K > S where only a subset of the positions in sequence K are taken into account (see **spaced kmer file and how to compute them**). The positions that are considered to build sequence S are based on a mask. This mask consists of a string of 1's and 0's. Masks are always of length K and have S 1's. Only the positions with 1's in the mask are extracted from sequence K to build S.
  
  JellyBelly scans sequence data and extracts all spaced kmers given a mask, a kmer length and a spaced kmer length. JellyBelly then builds a vector of spaced kmer counts. These vectors can then be used to do sequence comparison via vector calculus.
  
# Dependencies
  JellyBelly has no mayor dependencies other than zlib.
  
  JellyBelly uses [klib](https://github.com/attractivechaos/klib) (included in this repo) for handling sequence data and hashing of kmer sequences.
  
  JellyBelly uses cmake and make for compilation.

# Get and compile
    git clone --recursive https://github.com/7PintsOfCherryGarcia/JellyBelly.git
    cd JellyBelly
    mkdir build; cd build
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

	-q <int>	Number of spaced kmer vectors to keep in memory before
	  		writting them to the outputfile. (-q 100)

	-o <filename>	Output filename. (-o /dev/stdout)

	-b 		Binary output. (OFF)

	-r 		Raw output: spaced kmer counts (OFF)

	-l 		Genome mode. A single spaced kmer vector is
	  		computed for all input sequences. (OFF)

	-h 		This help  message.
# Usage

## Input Data
  JellyBelly can be executed on any type of DNA nucleotide data in the FASTA/FASTQ format. In addition a "spaced kmer file" (explained later) is also needed. 
  
  For example:
  
  Assuming the sequencing library "my_seq.fq.gz", and the spaced kmer file "SpacedKmer_K10_S5.bin" you can run JellyBelly in the following way:
    
    JellyBelly -f my_seq.fq.gz -s SpacedKmer_K10_S5.bin > output.txt

output.txt will be a tab separated text file with a count between 0 and 1 inclusive for each spaced kmer. The number of counts (number of spaced kmers) depends on the ones in the spaced kmer mask and is equal to 4^S where S is the the number of ones in the mask.

## Output Data
  
  JellyBelly's output values are scaled spaced kmer counts ranging from 0 to 1. Spaced kmer count values are sorted lexicographcally eg. first value corresponds to AAA..A second to AAA..T and so on.
  
  output.txt
  
      1	  2	  3	  	4^S
    
    0.245	0.014	0.547	...    0.436


  JellyBelly writes its output to either stdout (default) or to a file indicated by the -o option. Vectors are written row wise (one vector per row) and each vector corresponds to each input sequence in the same order. If the -l(genome mode) flag is set, only a sinlge vector is written as output.
  
  You can set the -b (binary output) option to make JellyBelly write the output vectors in binary form. This will make the output file much smaller. If the -b option is set. JellyBelly will write additional information to the output file. This additional information consists of a header and a tail following this format:

header | 26 bytes
------------ | --------
7 zero bytes  | All of JellyBelly's output binary files start with 7 bytes set to 0.
4 bytes (int) | The next 4 bytes encode the spaced kmer length.
3 zero bytes | After the spaced kmer length, an additional 3 bytes are set to 0.
4 bytes (int) | The next 4 bytes encode the kmer length.
2 zero bytes | After the kmer length, an additional 2 bytes are set to 0.
1 byte (char) | The output type: scaled(0) or raw(1) is coded in the value of this byte.
5 bytes (5 chars) | The header ends in 5 bytes set to the values 74,69,76,76,89 respectively.

After the header. You will find the vector values followed by a tail.

tail | 17 bytes
------------ | --------
4 zero bytes | After vector data, 4 zero bytes are written.
8 bytes (long)  | Number of samples in the current file.
5 bytes (5 chars) | The tail end in 5 bytes set to the values 66,69,76,76,89 respectively.

Check [python](https://github.com/7PintsOfCherryGarcia/JellyBelly/tree/master/python) for detailed instructions on how to load JellyBelly's output vectors into numpy arrays for downstream analysis.

## Genome mode  
  If you are running JellyBelly on assembled contigs or reference genomes. Every sequence in a fasta file will be encoded into a spaced kmer vector. If you want a set of sequences to be encoded into a single vector make sure to use the -l flag (genome mode).
  
  
# spaced kmer file and how to compute them
  A core concept of JellyBelly is the use of a "mask" to only consider specific positions of a kmer. For example given the following kmer of length 10:
  
                                                   ACGGTGCAAT
						   
  and the following mask:
  
                                                   1001100001
						   
  The corresponding spaced kmer will be:
  
                                                   A  GT    T == AGTT

Briefly, given a kmer size K and a spaced kmer size S, there are S choose K different ways of selecting S positions out of a kmer of length K. JellyBelly includes another program to compute these masks called create_spaces. create_spaces will do a brute force search of all possible mask configurations. Subsequently, only the most entropic(dissorderly) masks are kept.

To create a spaced kmer file simply run create_spaces with your desired kmer size and spaced kmer file.

    create_spaces <kmer length> <smer length>
	
This will create a .bin file. This is a binary file that contains the entropic masks. This file should be used as the -s option for JellyBelly.

JellyBelly uses the same mask for all input sequences.

JellyBelly is hardcoded to use the first mask in the .bin file as of now. Options will be added to make this costumizable.

# Memory management

JellyBelly is designed to have a low system memory footprint. For most use cases JellyBelly can run on your average computer. 4Gb of system memory can handle most use cases. There are 2 main aspects that affect memory consumption: spaced kmer size and the buffer size (-q flag). 

As the lenght of the spaced kmer increases so does the memory required for the hash table and the output buffer except these grows exponentially. More specifically, given a spaced kmer size of length S, 4^S sequences must be stored plus the additional space needed for the hash table.

The buffer size controls how many sequences are vectorized before writing the output to disk. If you don;t have much memory in your computer, try to keep this parameter low.

# Other tools

tool location | function
------------ | --------
python/compute_EUdist.py | Given two spaced kmer vectors, computes the eucledian distance between them
python/create_PCAplot.py | Given a set of spaced kmer vectors and sample infomation, computes a PCA and plots a projection into the first two compionents
python/create_UMAPplot.py| Given a set of spaced kmer vectors and sample information, computes a UMAP transformation and plots the first two axes
build/distMatrix | Given a set of spaced kmer vecotrs, computes all pairwise euclidean distances and prints an nxn matrix with n being the number of vectors provided. The diagonal is filled with zeroes.

Check [python](https://github.com/7PintsOfCherryGarcia/JellyBelly/tree/master/python) for more information.
    


# Upcoming changes

* Spaced kmer files
    * Currewnt version of create_spaces keeps only the spaced kmer masks with maximum entropy. Future versions will keep all spaced kmer masks that a within a certain fraction of the most entropic one.
    * JellyBelly uses only the first spaced kmer mask. Future versions will let the user choose which mask to use. This will require another utility program to list all the masks in a spaced kmer file.
* Output of JellyBelly
    * Currently JellyBelly outputs spaced kmer scaled counts in text format. Future versions will output this information in binary format to a filename specified as an option. Converting binary output to text output will be implemented separately.
* JellyBelly's execution
    * Currently JellyBelly has only one way of executing it. Future versions will introduce 'commands'. For example, there will be a comand to vectorize sequence data, another comand to print binary vector files, and another comand scale raw counts.

# Other stuff

* protein sequences
    * Currently not supported
* Use of multiple threads
    * The speed advantages of multithreading do not justify the resources I would have to invest in implementing it. Will be implemented in the near future.
* Use in windows
    * No
* I used this tool and don't know how to cite this work
    * Me neither, not planning on writing a paper so github link should suffice. Alternatively you could cite a presentation I gave at the [EMBL PhD symposium](http://phdsymposium.embl.org/) where I presented this work. In a talk titled "Deep Autoencoders for Sequence Similarity Exploration"
