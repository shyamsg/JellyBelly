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
*/
#include <zlib.h>

// Li's klib
#include "khash.h"
#include "kseq.h"

#include "main.h"

//TODO Document
#define pair_lt(a, b) ((strcmp((a).key,(b).key) < 0)?1:0)
#define pair_sc(a, b) ((a) < (b))


//Initialize Heng Li's khash
KHASH_MAP_INIT_STR(smer, unsigned int)

// initialize Heng Li's kseq
KSEQ_INIT(gzFile, gzread)



/*
  Struct storing spaced kmers in order to access a sorted array of
  them and have consitent printing
*/
typedef struct {
  char *key;
} SpKMER;


/*
 Struct for storing kmer length value spaced kmer length value, and the mask to be used
*/
typedef struct {
    //Sequence variables
    char *kmer_seq; //kmer sequence
    int kmerlength; //kmer length
    char *smer_seq; //spaced kmer sequence
    int smerlength; //spaced kmer length
    unsigned long int smernum;
    //Output variables
    int *mask; //spaced kmer mask positions
    int buffersize;
    float *routput;
    float *soutput;
    FILE *output;
} jellydata;


/*
 Struct for storing hash variables
*/
typedef struct {
    khash_t(smer) *h;
    khint_t k;
    unsigned long int hashsize;
} jellyhash;


int belly_start(gzFile fp,
                FILE *smer_file,
                jellyopts opts);


int belly_jellyinit(jellydata *jdata, jellyopts opts, FILE *smer_file);


int belly_read_header(FILE *file,
                      jellydata *info);


int check_zeros(char *zero_vector,
                int length);


/*Sets "mask" to corresponding values given "kmerlength" kmer length and
  "spacelength" number of spaces.
  input:
  int kmerlength -  Length of kmer to use
  int smernum - Number of spaced kmers contained in spaced kmer file
  FILE *smer_file - file pointer to spaced kmer file
  int *mask - Positions of the kmer to use
  returns:
  void
*/
int belly_get_mask(FILE *smer_file,
                   jellydata *info,
                   int mode);


char *get_smer(char *array,
               int kmerlen,
               int smernum);


int belly_hashinit(SpKMER **smerlist, jellyhash *smerhash, jellydata jdara);


int randint(int max);


unsigned int belly_hash_fill(jellyhash *smerhash,
                             jellydata jdata,
                             SpKMER *smerlist);


void belly_exit(gzFile fp,
                kseq_t *seq,
                jellyhash smerhash,
                jellydata *info,
                SpKMER *smerlist);


void belly_fill_smer_list(SpKMER *smerlist,
                          jellydata info,
                          int idx,
                          unsigned int *smer_idx);


unsigned long belly_extract_spaces(kseq_t *seq,
                                   jellydata *info,
                                   jellyhash *smerhash,
                                   SpKMER *smerlist,
                                   jellyopts mode);


unsigned long int belly_count(char *seq,
                              int l,
                              jellydata *info,
                              jellyhash *smerhash);


void belly_extract_vector(jellyhash *smerhash,
                          unsigned long vector_length,
                          SpKMER *smerlist,
                          unsigned int *vector);


int belly_scale(unsigned int *vec,
                unsigned int len,
                float *scvec,
                unsigned long int *svidx);


unsigned int belly_max(unsigned int *vector,
                       int length);


unsigned int belly_min(unsigned int *vector);


int belly_allocateinfo(jellydata *info,
                       int mode);


void belly_vectorout(float *vector,
                     unsigned int size,
                     FILE *outfile);

int belly_minmax(unsigned int *vector, int length, int *min, int *max);
