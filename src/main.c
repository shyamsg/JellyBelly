/*
MIT License

Copyright (c) 2019 Julian Regalado Perez

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*
################################################################################
JetBelly - Extract and bin kmers from fasta sequences of fastq sequencing files
           into spaced kmers
*/


#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>

#include "main.h"


int main(int argc, char **argv)
{
    fprintf(stderr,"%d\n", argc);
    JellyOpts opts;

    // Check for correct number of arguments
    if (argc < 5) {
      fprintf(stderr,"Not enough parameters.\n");
      usage();
    }

    // Read options
    read_opts(argc, argv, &opts);

    gzFile fp;
    // Read from stdin
    if (strcmp( opts.seqfilename, "-") == 0) {
      opts.seqfilename = "/dev/stdin";
    }
    fp = gzopen(opts.seqfilename, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: Can not open input file.\n");
        usage();
    }

    FILE *smer_file = fopen(opts.smerfilename,"rb");
    if (! smer_file) {
        fprintf(stderr,"ERROR: Can not open smer kmer file.\n");
        usage();
    }

    // Run
    int ret = belly_start(fp, smer_file, opts);
    fclose(smer_file);
    return(ret);
}


/*

*/
void usage()
{
    fprintf(stderr,"Usage:\n");
    fprintf(stderr,"\tJellyBelly [options] ");
    fprintf(stderr,"-f <sequence file> -s <spacedkmer file>\n\n");
    fprintf(stderr,"Options:\n\n");
    fprintf(stderr,"\t-f <filename>\tSequence file. fasta/q file with sequences or\n");
    fprintf(stderr,"\t  \t\tsequencing reads. Can be gziped. provide \"-\" if\n");
    fprintf(stderr,"\t  \t\treading from stdin. If using sequencing data eg.\n");
    fprintf(stderr,"\t  \t\tillumina reads, provide the -C flag to use canonical\n");
    fprintf(stderr,"\t  \t\tform.\n\n");
    fprintf(stderr,"\t-s <filename>\tSpaced kmer file. Binary file containing the masks\n");
    fprintf(stderr,"\t  \t\tto be used. Refer to manual for detailed documentation.\n\n");
    fprintf(stderr,"\t-C \t\tCanonical mode. Lexicographically smallest kmer is ");
    fprintf(stderr,"counted.\n\t  \t\tSet this flag when analyzing sequencing reads.\n\n");
    fprintf(stderr,"\t-b <int>\tNumber of spaced kmer vectors to keep in memory before\n");
    fprintf(stderr,"\t  \t\twritting them to the outputfile.\n\n");
    fprintf(stderr,"\t-h \t\tThis help  message.\n\n");
    exit(1);
}


/*

*/
void read_opts(int argc, char **argv, JellyOpts *opts)
{
    opts->seqfilename = NULL;
    opts->smerfilename = NULL;
    opts->mode = 0;
    opts->buffersize = 100;
    int elem;
    while (( elem = getopt(argc, argv, "f:s:Chb:") ) >= 0) {
        switch(elem) {
        case 'f':
            opts->seqfilename = optarg;
            break;
        case 's':
            opts->smerfilename = optarg;
            break;
        case 'C':
            opts->mode = 1;
            break;
        case 'b':
            opts->buffersize = atoi(optarg);
            break;
        case 'h':
            usage();
        }
    }
    if (!opts->seqfilename || !opts->smerfilename) {
        fprintf(stderr,"\t ERROR: Please provide a sequence file (-f) and\n");
        fprintf(stderr,"\t        a spaced kmer file file (-s)\n\n");
        usage();
    }

    if (opts->buffersize <= 0) {
        fprintf(stderr,"\tERROR: Please specify a buffersize >= 0\n");
        usage();
    }
    print_opts(*opts);
}


/*

*/
void print_opts(JellyOpts opts)
{
    fprintf(stderr,"\t sequence file: %s\n", opts.seqfilename);
    fprintf(stderr,"\t spaced kmer file: %s\n", opts.smerfilename);
    fprintf(stderr,"\t kmer mode: %d\n", opts.mode);
    fprintf(stderr,"\t buffer size: %d\n", opts.buffersize);
}

