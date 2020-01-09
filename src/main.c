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
#include <unistd.h>
#include <poll.h>
#include "main.h"


int main(int argc, char **argv)
{
    print_info();

    jellyopts opts;
    // Read options
    read_opts(argc, argv, &opts);

    // Check for correct number of arguments
    //if (argc < 5) {
    //    fprintf(stderr,"\tERROR: Not enough parameters.\n");
    //    usage();
    //}


    gzFile fp;
    FILE *smer_file;

    // Read from stdin
    if (strcmp( opts.seqfilename, "-") == 0) {
        opts.seqfilename = "/dev/stdin";
        if (check_stdin(opts.seqfilename) == 0) {
            fprintf(stderr,"\tERROR: -f - set but no data from stdin was detected.\n");
            usage();
        }
    }
    fp = gzopen(opts.seqfilename, "r");
    if (!fp) {
        fprintf(stderr, "\tERROR: Can not open input file.\n");
        usage();
    }

    smer_file = fopen(opts.smerfilename, "rb");
    if (! smer_file) {
        fprintf(stderr,"\tERROR: Can not open smer kmer file.\n");
        usage();
    }

    // Run
    if (!belly_start(fp, smer_file, opts)) {
        fprintf(stderr, "\tERROR: JelyBelly ran unsuccesfully.\n");
        fclose(smer_file);
        return -1;
    }
    fclose(smer_file);
    return 0;
}


void usage()
{
    fprintf(stderr,"\n\nUsage:\n");
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
    //fprintf(stderr,"\t-C \t\tCanonical mode. Lexicographically smallest kmer is counted.\n");
    //fprintf(stderr,"\t  \t\tSet this flag when analyzing sequencing reads. (OFF)\n\n");
    fprintf(stderr,"\t-q <int>\tNumber of spaced kmer vectors to keep in memory before\n");
    fprintf(stderr,"\t  \t\twritting them to the outputfile. (-q 100)\n\n");
    fprintf(stderr,"\t-o <filename>\tOutput filename. (-o /dev/stdout)\n\n");
    fprintf(stderr,"\t-r \t\tRaw output: spaced kmer counts (OFF)\n\n");
    fprintf(stderr,"\t-l \t\tGenome mode. A single spaced kmer vector is\n");
    fprintf(stderr,"\t  \t\tcomputed for all input sequences. (OFF)\n\n");
    fprintf(stderr,"\t-h \t\tThis help  message.\n\n");
    exit(-1);
}


void read_opts(int argc, char **argv, jellyopts *opts)
{
    opts->seqfilename = NULL;
    opts->smerfilename = NULL;
    opts->outputfilename = NULL;
    opts->mode = 0;
    opts->buffersize = 100;
    opts->binout = 0;
    opts->raw = 0;
    opts->gmode = 0;
    int elem;
    while (( elem = getopt(argc, argv, "f:s:Chq:o:brl") ) >= 0) {
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
        case 'q':
            opts->buffersize = atoi(optarg);
            break;
        case 'o':
            opts->outputfilename = optarg;
            break;
        case 'b':
            opts->binout = 1;
            break;
        case 'r':
            opts->raw = 1;
            break;
        case 'l':
            opts->gmode = 1;
            break;
        case 'h':
            usage();
        }
    }
    if (!opts->seqfilename || !opts->smerfilename) {
        fprintf(stderr,"\tERROR: Please provide a sequence file (-f) and\n");
        fprintf(stderr,"\t        a spaced kmer file file (-s)\n\n");
        usage();
    }

    if (opts->buffersize <= 0) {
        fprintf(stderr,"\tERROR: Please specify a buffersize >= 0\n");
        usage();
    }

    if (!opts->outputfilename) opts->outputfilename = "/dev/stdout";
    print_opts(*opts);
}


void print_opts(jellyopts opts)
{
    fprintf(stderr,"\tINFO: sequence file:\t\t%s\n", opts.seqfilename);
    fprintf(stderr,"\tINFO: spaced kmer file:\t\t%s\n", opts.smerfilename);
    fprintf(stderr,"\tINFO: output file:\t\t%s\n", opts.outputfilename);
    fprintf(stderr,"\tINFO: kmer mode:\t\t%s\n", (opts.mode==0)?"non-canonical":"canonical");
    fprintf(stderr,"\tINFO: buffer size:\t\t%d\n", opts.buffersize);
    fprintf(stderr,"\tINFO: genome mode:\t\t%s\n", (opts.gmode==1)?"ON":"OFF");
    fprintf(stderr,"\tINFO: output format:\t\t%s\n", (opts.binout==0)?"text":"binary");
    fprintf(stderr,"\tINFO: raw output:\t\t%s\n\n", (opts.raw==1)?"ON":"OFF");
}


int check_stdin(char *filename)
{
    FILE *fp = fopen("/dev/stdin", "r");
    if (!fp) return 0;
    struct pollfd *fpoll;
    fpoll = malloc(sizeof(fpoll));
    fpoll[0].fd = fileno(fp);
    fpoll[0].events = POLLIN;
    fpoll[0].revents = 0;
    int stdinpoll = poll(fpoll, 1, 100);
    free(fpoll);
    fclose(fp);
    return stdinpoll;
}

void print_info()
{
    fprintf(stderr, "JellyBelly v0.0\n");
    fprintf(stderr, "\tMIT License\n");
    fprintf(stderr, "\tCopyright (c) 2019 Julian Regalado [julian.perez@bio.ku.dk]\n");
    fprintf(stderr, "\thttps://github.com/7PintsOfCherryGarcia/JellyBelly\n\n");
}
