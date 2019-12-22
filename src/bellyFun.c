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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>

// Functions
#include "bellyFun.h"
#include "bellyHash.h"
#include "bellyMisc.h"


//Initialize ksort for scaled vector
KSORT_INIT(pairSC, unsigned int, pair_sc)
KSORT_INIT_GENERIC(int)


int belly_start(gzFile fp, FILE *smer_file, JellyOpts opts)
{
    srand(time(NULL));

    jellydata info;
    info.kmer_seq = NULL;
    info.smer_seq = NULL;
    info.mask = NULL;
    info.routput = NULL;
    info.soutput = NULL;
    info.output = NULL;
    info.buffersize = opts.buffersize;

    //Check for output file
    info.output = fopen(opts.outputfilename, opts.omode);
    if (!info.output) {
      fprintf(stderr, "ERROR: Could not open output file: %s.\n",
              opts.outputfilename);
      if (errno == 13) fprintf(stderr, "Permission denied\n");
      return 1;
    }

    //Read spaced kmer file
    if (!belly_read_header(smer_file, &info, opts.mode)) {
        fprintf(stderr, "ERROR: Corrupted bin file.\n");
        return 1;
    }
    fprintf(stderr,"kmer %d smer %d snum %lu\n",info.kmerlength,
                                                info.smerlength,
                                                info.smernum);

    //if (!belly_allocateinfo(&info, opts.mode)) {
    //    fprintf(stderr,"ERROR: Could not allocate enough memory.\n");
    //    return 1;
    //}

    //Get spaced kmer mask
    //if (!belly_get_mask(smer_file, &info)) {
    //    fprintf(stderr, "ERROR: Could not read spaced kmer mask.\n");
    //    return 1;
    //}

    //Create spaced kmer struct and sort lexicographically
    unsigned int hashsize = pow(4, info.smerlength);
    fprintf(stderr, "Generating %u spaced kmers\n", hashsize);
    SpKMER *smerlist;
    smerlist = malloc(hashsize*sizeof(SpKMER));
    if (!smerlist) {
        fprintf(stderr, "ERROR: Failed to allocate enough memory.\n");
        return 1;
    }
    unsigned int smer_idx = 0;
    belly_fill_smer_list(smerlist, &info, 0, &smer_idx);
    fprintf(stderr, "smer_index: %u\n", smer_idx);
    fprintf(stderr,"first: %s\nlast: %s\n",smerlist[0].key, smerlist[hashsize-1].key);


    jellyhash smerhash;
    hash_declare(&smerhash);
    //After initialized hash table, fill with spaced kmers
    fprintf(stderr,"building hash table\n");
    hashsize = belly_hash_init(&smerhash, &info, smerlist);
    if(!hashsize) {
        fprintf(stderr,"ERROR: Could not initilize hash\n");
        return 1;
    }

    //Open sequence file
    kseq_t *seq;
    seq = kseq_init(fp);
    if (!seq) {
        fprintf(stderr,"ERROR: Could not open sequence file.\n");
        return 1;
    }
    fprintf(stderr,"binning kmers\n");

    unsigned long bases =  belly_extract_spaces(seq, &info, &smerhash, smerlist, opts);

    fprintf(stderr,"%lu bases\n",bases);
    belly_exit(fp, seq, smerhash, &info, smerlist, hashsize);
    return 0;
}


int belly_read_header(FILE *file, jellydata *info, int mode)
{
    int bytes_read;
    char zero_vector[19];
    void *header = malloc(48);
    if (!header) return 0;
    bytes_read = fread(header, 1, 48, file);
    if (bytes_read != 48) {
        free(header);
        return 0;}

    //Read int
    memcpy(&(info->kmerlength), header + 7, 4);
    //Read int
    memcpy(&(info->smerlength), header + 14, 4);
    //Read long
    memcpy(&(info->smernum), header + 21, 8);

    memcpy(zero_vector, header, 7);
    if (!check_zeros(zero_vector, 7)) {
        return 0;
    }

    memcpy(zero_vector, header + 11, 3);
    if (!check_zeros(zero_vector, 3)) {
        return 0;
    }

    memcpy(zero_vector, header + 18, 3);
    if (!check_zeros(zero_vector, 3)) {
        return 0;
    }

    memcpy(zero_vector, header + 29, 19);
    if (!check_zeros(zero_vector, 19)) {
        return 0;
    }
    free(header);

    //Get spaced kmer mask
    if (!belly_get_mask(file, info, mode)) {
      fprintf(stderr, "ERROR: Could not read spaced kmer mask.\n");
      return 1;
    }

    return 1;
}


int check_zeros(char *zero_vector, int length)
{
    for (int i = 0; i < length; i++) {
        if (zero_vector[i] != 0) {
            return 0;
        }
    }
    return 1;
}


int belly_get_mask(FILE *smer_file, jellydata *info, int mode)
{
    if (!belly_allocateinfo(info, mode)) {
        fprintf(stderr,"ERROR: Could not allocate enough memory.\n");
        return 1;
    }
    info->mask = malloc((info->smerlength + 1)*sizeof(int));
    if (!info->mask) return 0;

    char *smer_array = malloc(info->smernum*info->kmerlength);
    if (!smer_array) return 0;

    int bytes_read = fread(smer_array, 1, info->smernum*info->kmerlength, smer_file);
    if (bytes_read != (info->smernum*info->kmerlength)) {
        fprintf(stderr, "ERROR: Spaced kmer file corrupted.\n");
        free(smer_array);
        return 0;
    }

    char *smer = get_smer(smer_array, info->kmerlength, info->smernum);
    if (!smer) return 0;
    fprintf(stderr,"smer is: %s\n",smer);

    int j = 0;
    for (int i = 0; i < info->kmerlength; i++) {
        if (smer[i] == '1') {
        info->mask[j] = i;
        ++j;
        }
    }
    info->mask[j] = -1;
    free(smer);
    free(smer_array);

    //Get spaced kmer mask
    if (!belly_get_mask(smer_file, info, mode)) {
      fprintf(stderr, "ERROR: Could not read spaced kmer mask.\n");
      return 0;
    }

    return 1;
}


char *get_smer(char *array, int kmerlength, int smernum)
{

    //TODO Change name to belly_getsmer
    char *smer = malloc(kmerlength+1);
    if (!smer) return NULL;
    smer[kmerlength] = '\0';

    //smernum - 1 as this number is used as an index for char *array
    int rand_smernum = randint(smernum - 1);
    rand_smernum = 0;
    int kmeridx = rand_smernum*kmerlength;
    memcpy(smer, array + kmeridx, kmerlength);

    if(strlen(smer) != kmerlength) {
        fprintf(stderr, "ERROR: kmer length and spaced kmer mask of different lenght.\n");
        return NULL;
    }
    return smer;
}


int randint(int max)
{
    return 0 + rand() / (RAND_MAX / (max - 0 + 1) + 1);
}


unsigned int belly_hash_init(jellyhash *smerhash,
                             jellydata *info,
                             SpKMER *smerlist)
{
    int absent;
    for (unsigned long int i = 0; i < pow(4, info->smerlength); i++) {
        //fprintf(stderr,"|%s\n",smerlist[i].key);
        smerhash->k = kh_put(smer, smerhash->h, smerlist[i].key, &absent);
        // Set value to 0
        kh_value(smerhash->h, smerhash->k) = 0;
    }

    // Report hash size (number of elements in hash, must be equal to 4^smerlen)
    unsigned int size = kh_size(smerhash->h);
    fprintf(stderr,"hash size is %d\n",size);
    if (size != pow(4, info->smerlength)) {
        fprintf(stderr,"Error at generating spaced kmer hash.\n");
        return(1);
    }
    return(size);
}


void belly_fill_smer_list(SpKMER *smerlist,
                          jellydata *info,
                          int idx,
                          unsigned int *smer_idx)
{
    int i;
    if (idx < (info->smerlength - 1)) {
        for (i = 0; i < strlen(BASES); i++) {
            info->smer_seq[idx] = BASES[i];
            belly_fill_smer_list(smerlist, info, idx + 1, smer_idx);
        }
    }
    else {
        for (i = 0; i < strlen(BASES); i++) {
            info->smer_seq[idx] = BASES[i];
            // Add spaced kmer
            // There is your smer!!!
            smerlist[*smer_idx].key = malloc((info->smerlength + 1)*sizeof(char));
            smerlist[*smer_idx].key = memcpy(smerlist[*smer_idx].key,
                                             info->smer_seq,
                                             (info->smerlength + 1));
            *smer_idx += 1;
        }
    }
}


unsigned long belly_extract_spaces(kseq_t *seq,
                                   jellydata *info,
                                   jellyhash *smerhash,
                                   SpKMER *smer_list,
                                   JellyOpts opts)
{
    //TODO exit properly
    unsigned long hash_size = kh_size(smerhash->h);
    unsigned long int total_seq = 0;
    int sequence_length;
    unsigned long int n = 0;
    unsigned long int total_kmers = 0;
    unsigned int *hash_vector = malloc(hash_size*sizeof(unsigned int));
    if (!hash_vector) {
        total_seq = 0;
        goto exit;
    }

    //Vector holding scaled spaced kmer values for input sequences.
    float *scale_vector = malloc((hash_size*info->buffersize)*sizeof(float));
    if (!scale_vector) {
        total_seq = 0;
        goto exit;
    }
    //Index for scaled vector
    unsigned long int sv_idx = 0;
    unsigned long int numkmers = 0;

    //Read sequences in loop
    while ((sequence_length = kseq_read(seq)) >= 0) {
        //Several sanity checks
        //Check for valid sequence file format (fasta/fastq)
        if (sequence_length < 0) {
            fprintf(stderr,"ERROR: Could not read in sequence file.\n");
            total_seq = 0;
            goto exit;
        }
        //Check for sequence >= kmerlength
        else if (sequence_length < info->kmerlength) continue;
        //Check for empty sequence
        else if (sequence_length == 0) continue;
        //Check for valid kmers (Only containing A C G T)
        if (!(numkmers = belly_count(seq->seq.s, sequence_length, info, smerhash))) {
            fprintf(stderr, "WARNING: No valid kmers in sequence.\n");
            continue;
        }

        total_seq += sequence_length;
        total_kmers += numkmers;
        //Output spaced kmer vector for current sequence
        if (!opts.mode) {
            belly_extract_vector(smerhash, hash_size, smer_list, hash_vector);
            if (!belly_scale(hash_vector,
                             hash_size,
                             scale_vector,
                             hash_size,
                             &sv_idx,
                             numkmers)) {
                total_seq = 0;
                goto exit;
            }
            hash_reset(smerhash->h, smerhash->k);
        }
        n++;
    }

    if (n == 0) {
        fprintf(stderr,"Error in sequence file\n");
        total_seq = 0;
        goto exit;
    }
    //Output spaced kmer vector for file
    if (opts.mode) {
        sv_idx = 0;
        belly_extract_vector(smerhash, hash_size, smer_list, hash_vector);
        if (!belly_scale(hash_vector,
                         hash_size,
                         scale_vector,
                         hash_size,
                         &sv_idx,
                         total_kmers)) {
            total_seq = 0;
            goto exit;
        }
    }

    //Log info
    fprintf(stderr,"%lu total reads.\n",n);
    fprintf(stderr, "%lu total kmers\n", total_kmers);
    exit:
        free(scale_vector);
        free(hash_vector);
        return(total_seq);
}


/*
  Count spaced kmers in a given sequence
*/
unsigned long int belly_count(char *seq,
                              int l,
                              jellydata *info,
                              jellyhash *smerhash)
{
    unsigned long int kmer_num = 0;
    unsigned long int valid_kmer_num = 0;
    //Loop over possible kmers in sequence
    for (; kmer_num < l - info->kmerlength + 1; kmer_num++) {
        //Get kmer from sequence
        //TODO Instead of copying kmer sequence to new array, pass address directly to
        //hash_genSpacemer.
        info->kmer_seq = memcpy(info->kmer_seq, &seq[kmer_num], info->kmerlength);
        //Extract positions stored in mask
        hash_genSpacemer(info->kmer_seq, info->mask, info->smer_seq);
        //Check if spaced kmer already exists and increase number
        smerhash->k = kh_get(smer, smerhash->h, info->smer_seq);
        if (smerhash->k == kh_end(smerhash->h)) {
            continue;
        }
        else {
            kh_val(smerhash->h, smerhash->k) += 1;
            valid_kmer_num += 1;
        }
    }
    return valid_kmer_num;
}


void belly_extract_vector(jellyhash *smerhash,
                          unsigned long vector_length,
                          SpKMER *smerlist,
                          unsigned int *vector)
{
    for (unsigned long i = 0; i < vector_length; i++) {
        smerhash->k = kh_get(smer, smerhash->h, smerlist[i].key);
        vector[i] = kh_val(smerhash->h, smerhash->k);
    }
}


int belly_scale(unsigned int *vec,
                unsigned int hsize,
                float *scvec,
                int len,
                unsigned long int *svidx,
                unsigned long int nkmers)
{
    //TODO bundle diff, min and tmp into struct and pass by reference
    unsigned int diff;
    unsigned int min;
    unsigned int *tmp = malloc(len*sizeof(unsigned int));
    if (!tmp) return 0;
    memcpy(tmp, vec, len*sizeof(unsigned int));

    ks_mergesort(pairSC, len, vec, 0);
    min = belly_min(vec);
    diff = belly_max(vec, len) - min;
    fprintf(stderr,"min: %d\t max: %d\n", min, belly_max(vec, len));
    for (unsigned long i = 0; i < len; i++)
        scvec[(*svidx*hsize) + i] = (((float)tmp[i] - (float)min) / (float)diff);

    //TODO do the actual printing
    if (*svidx == 0) {
        *svidx=0;
        fprintf(stderr,"Dumping\n");
        belly_vectorout(scvec, hsize);
    }
    else {
        *svidx += 1;
    }
    free(tmp);
    return 1;
}


void belly_exit(gzFile fp,
                kseq_t *seq,
                jellyhash smerhash,
                jellydata *info,
                SpKMER *smerlist,
                unsigned int hashsize)
{
  // Close zlib file handler
  gzclose(fp);
  // Close kseq pointer
  kseq_destroy(seq);
  // Free hash
  kh_destroy(smer, smerhash.h);
  // free spaced kmer arrays
  for(int i = 0; i < hashsize; i++) {
    free(smerlist[i].key);
  }
  free(smerlist);
  free(info->kmer_seq);
  free(info->smer_seq);
  free(info->mask);
  fclose(info->output);
}


unsigned int belly_max(unsigned int *vector, int length)
{
    return vector[length - 1];
}


unsigned int belly_min(unsigned int *vector) {
    return vector[0];
}


void belly_vectorout(float *vector, unsigned int size)
{
    unsigned int pos = 0;
    //for (unsigned int j = 0; j < 1; j++) {
    //pos = j*size;
    fprintf(stdout,"%.8f",vector[pos]);
    for (unsigned long int i = 1; i < size; i++) {
        //Change to print to file and change between binary and text
        fprintf(stdout,"\t%.8f",vector[pos + i]);
    }
    fprintf(stdout, "\n");
    //}
}


int belly_allocateinfo(jellydata *info, int mode)
{
    if (!mode) {
        fprintf(stderr,"Running on contig sequence mode.\n");
    }
    else {
        fprintf(stderr,"Running on read sequence mode.\n");
    }

    // Declare kmer and spaced kmer array
    info->kmer_seq = malloc(((info->kmerlength)+1)*sizeof(char));
    if (!info->kmer_seq) {
        fprintf(stderr,"ERROR: unable to allocate memory for kmer string.\n");
        return 0;
    }
    info->kmer_seq[info->kmerlength] = '\0';

    info->smer_seq = malloc(((info->smerlength)+1)*sizeof(char));
    if (!info->smer_seq) {
        fprintf(stderr,"ERROR: unable to allocate memory for spaced kmer string.\n");
        return 0;
    }
    info->smer_seq[info->smerlength] = '\0';

    return 1;
}
