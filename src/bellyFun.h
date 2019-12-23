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
#include "ksort.h"

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
    char *kmer_seq;
    int kmerlength;
    char *smer_seq;
    int smerlength;
    unsigned long int smernum;
    int *mask;
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
} jellyhash;


//TODO Document
int belly_start(gzFile fp,
                FILE *smer_file,
                JellyOpts opts);


int belly_read_header(FILE *file,
                      jellydata *info,
                      int mode);


int check_zeros(char *zero_vector, int length);

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
int belly_get_mask(FILE *smer_file, jellydata *info, int mode);

char *get_smer(char *array, int kmerlen, int smernum);

int randint(int max);


//TODO Document
unsigned int belly_hash_init(jellyhash *smerhash, jellydata *info, SpKMER *smerlist);

//TODO Document
void belly_exit(gzFile fp,
                kseq_t *seq,
                jellyhash smerhash,
                jellydata *info,
                SpKMER *smerlist,
                unsigned int hashsize);


void belly_fill_smer_list(SpKMER *smerlist, jellydata *info ,int idx, unsigned int *smer_idx);


unsigned long belly_extract_spaces(kseq_t *seq,
                                   jellydata *info,
                                   jellyhash *smerhash,
                                   SpKMER *smerlist,
                                   JellyOpts mode);


unsigned long int belly_count(char *seq,
                              int l,
                              jellydata *info,
                              jellyhash *smerhash);


void belly_extract_vector(jellyhash *smerhash,
                          unsigned long vector_length,
                          SpKMER *smerlist,
                          unsigned int *vector);

int belly_scale(unsigned int *vector,
                unsigned int hashsize,
                float *scaled_vector,
                int length,
                unsigned long int *sv_idx,
                unsigned long int numkmers);


unsigned int belly_max(unsigned int *vector, int length);


unsigned int belly_min(unsigned int *vector);


int belly_allocateinfo(jellydata *info, int mode);


void belly_vectorout(float *vector, unsigned int size);
