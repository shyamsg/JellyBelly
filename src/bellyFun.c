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

//TODO Memory leak when using -r(raw output option)

int belly_start(gzFile fp, FILE *smer_file, jellyopts opts)
{
    fprintf(stderr, "Running JellyBelly.\n");
    srand(time(NULL));
    int ret;
    //Initialize JellyBelly's data structures
    jellydata jdata;
    if (!belly_jellyinit(&jdata, opts, smer_file)) {
        fprintf(stderr, "\tERROR: Could not initialize JellyBelly.\n");
        ret = -1;
        goto exit;
    }


    //Load spaced kmers and initialize hash table
    SpKMER *smerlist;
    jellyhash smerhash;
    if (!belly_hashinit(&smerlist, &smerhash, jdata)) {
        fprintf(stderr, "\tERROR: Could not initialize hash table.\n");
        ret = -1;
        goto exit;
    }

    //Open sequence file
    kseq_t *seq;
    seq = kseq_init(fp);
    if (!seq) {
        fprintf(stderr,"\tERROR: Could not open sequence file.\n");
        ret = -1;
        goto exit;
    }

    if ( !belly_extract_spaces(seq, &jdata, &smerhash, smerlist, opts)) {
      fprintf(stderr, "\tERROR: Failed to count kmers.\n");
      ret = -1;
      goto exit;
    }
    ret = 1;
    exit:
      belly_exit(fp, seq, smerhash, &jdata, smerlist);
      return ret;
}


int belly_jellyinit(jellydata *jdata, jellyopts opts, FILE *smer_file)
{
    jdata->kmer_seq = NULL;
    jdata->smer_seq = NULL;
    jdata->mask = NULL;
    jdata->routput = NULL;
    jdata->soutput = NULL;
    jdata->output = NULL;
    jdata->buffersize = opts.buffersize;

    //Read spaced kmer file
    if (!belly_read_header(smer_file, jdata)) {
        fprintf(stderr, "\tERROR: Corrupted bin file.\n");
        return 0;
    }
    fprintf(stderr,"\tINFO: kmer size: %d, smer size: %d snum %lu\n",
            jdata->kmerlength,
            jdata->smerlength,
            jdata->smernum);

    //TODO Add error checking for pow() operation
    jdata->hashsize = (unsigned int)pow(4, jdata->smerlength);

    if (!belly_allocateinfo(jdata, opts.raw, opts.gmode)) {
        fprintf(stderr,"\tERROR: Could not allocate enough memory.\n");
        return 0;
    }

    //Get spaced kmer mask
    if (!belly_get_mask(smer_file, jdata, opts.mode)) {
        fprintf(stderr, "\tERROR: Could not read spaced kmer mask.\n");
        return 0;
    }

    if (opts.binout) {
        jdata->output = fopen(opts.outputfilename, "wb");
        if (!belly_checkfile(jdata->output)) return 0;
        if (!belly_ofileinit(jdata)) {
            fprintf(stderr, "\tERROR: Failed to initialize output file.\n");
            return 0;
        }
    }
    else {
        jdata->output = fopen(opts.outputfilename, "w");
        if (!belly_checkfile(jdata->output)) return 0;
    }


    return 1;
}


int belly_read_header(FILE *file,
                      jellydata *jdata)
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
    memcpy(&(jdata->kmerlength), header + 7, 4);
    //Read int
    memcpy(&(jdata->smerlength), header + 14, 4);
    //Read long
    memcpy(&(jdata->smernum), header + 21, 8);

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

    return 1;
}


int check_zeros(char *zero_vector,
                int length)
{
    for (int i = 0; i < length; i++) {
        if (zero_vector[i] != 0) {
            return 0;
        }
    }
    return 1;
}


int belly_allocateinfo(jellydata *jdata, int raw, int gmode)
{

    unsigned int size = jdata->hashsize;
    if (gmode) {
        if (raw) {
            jdata->routput = malloc(size*sizeof(unsigned int));
            if (!jdata->routput) {
                fprintf(stderr, "\tERROR: Could not allocate data for output array.\n");
                fprintf(stderr, "\t\ttry decreasing buffer size (-q)\n");
                return 0;
            }
        }
        else {
            jdata->soutput = malloc(size*sizeof(float));
            if (!jdata->soutput) {
                fprintf(stderr, "\tERROR: Could not allocate data for output array.\n");
                fprintf(stderr, "\t\ttry decreasing buffer size (-q)\n");
                return 0;
            }
        }
    }
    else {
        if (raw) {
            jdata->routput = malloc((jdata->buffersize*size)*sizeof(unsigned int));
            if (!jdata->routput) {
                fprintf(stderr, "\tERROR: Could not allocate data for output array.\n");
                fprintf(stderr, "\t\ttry decreasing buffer size (-q)\n");
                return 0;
            }
        }
        else {
            jdata->soutput = malloc((jdata->buffersize*size)*sizeof(float));
            if (!jdata->soutput) {
                fprintf(stderr, "\tERROR: Could not allocate data for output array.\n");
                fprintf(stderr, "\t\ttry decreasing buffer size (-q)\n");
                return 0;
            }
        }
    }

    jdata->kmer_seq = malloc(((jdata->kmerlength)+1)*sizeof(char));
    if (!jdata->kmer_seq) {
        fprintf(stderr,"\tERROR: unable to allocate memory for kmer string.\n");
        return 0;
    }
    jdata->kmer_seq[jdata->kmerlength] = '\0';

    jdata->smer_seq = malloc(((jdata->smerlength)+1)*sizeof(char));
    if (!jdata->smer_seq) {
        fprintf(stderr,"\tERROR: unable to allocate memory for spaced kmer string.\n");
        return 0;
    }
    jdata->smer_seq[jdata->smerlength] = '\0';

    return 1;
}


int belly_get_mask(FILE *smer_file,
                   jellydata *jdata,
                   int mode)
{

    jdata->mask = malloc((jdata->smerlength + 1)*sizeof(int));
    if (!jdata->mask) return 0;

    char *smer_array = malloc(jdata->smernum*jdata->kmerlength);
    if (!smer_array) return 0;

    int bytes_read = fread(smer_array, 1, jdata->smernum*jdata->kmerlength, smer_file);
    if (bytes_read != (jdata->smernum*jdata->kmerlength)) {
        fprintf(stderr, "\tERROR: Spaced kmer file corrupted.\n");
        free(smer_array);
        return 0;
    }

    char *smer = get_smer(smer_array, jdata->kmerlength, jdata->smernum);
    if (!smer) return 0;
    fprintf(stderr,"\tINFO: smer is: %s\n",smer);

    int j = 0;
    for (int i = 0; i < jdata->kmerlength; i++) {
        if (smer[i] == '1') {
            jdata->mask[j] = i;
            ++j;
        }
    }
    jdata->mask[j] = -1;
    free(smer);
    free(smer_array);

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


int belly_hashinit(SpKMER **smerlist, jellyhash *smerhash, jellydata jdata)
{
    smerhash->hashsize = pow(4, jdata.smerlength);
    fprintf(stderr, "\tINFO: Generating %lu spaced kmers\n", smerhash->hashsize);

    *smerlist = malloc((smerhash->hashsize)*sizeof(SpKMER));
    if (!smerlist) {
        fprintf(stderr, "ERROR: Failed to allocate enough memory.\n");
        return 0;
    }
    unsigned int smer_idx = 0;
    belly_fill_smer_list(*smerlist, jdata, 0, &smer_idx);
    //Initialize hash table and fill with spaced kmers
    hash_declare(smerhash);
    fprintf(stderr,"\tINFO: Building hash table\n");
    if (smerhash->hashsize != belly_hash_fill(smerhash, jdata, *smerlist)) {
        fprintf(stderr,"ERROR: Could not initilize hash\n");
        return 0;
    }

    return 1;
}


unsigned int belly_hash_fill(jellyhash *smerhash,
                             jellydata jdata,
                             SpKMER *smerlist)
{

    int absent;
    //Add spaced kmer sequence to hash and set value to 0
    for (unsigned long int i = 0; i < pow(4, jdata.smerlength); i++) {
        smerhash->k = kh_put(smer, smerhash->h, smerlist[i].key, &absent);
        kh_value(smerhash->h, smerhash->k) = 0;
    }

    // Report hash size (number of elements in hash, must be equal to 4^smerlen)
    unsigned long size = jdata.hashsize;
    fprintf(stderr,"\tINFO: Hash table size is %lu\n",size);
    if (size != pow(4, jdata.smerlength)) {
        fprintf(stderr,"Error at generating spaced kmer hash.\n");
        return 0;
    }

    return(size);
}


void belly_fill_smer_list(SpKMER *smerlist,
                          jellydata info,
                          int idx,
                          unsigned int *smer_idx)
{
    if (idx < (info.smerlength - 1)) {
        for (int i = 0; i < strlen(BASES); i++) {
            info.smer_seq[idx] = BASES[i];
            belly_fill_smer_list(smerlist, info, idx + 1, smer_idx);
        }
    }
    else {
        for (int i = 0; i < strlen(BASES); i++) {
            info.smer_seq[idx] = BASES[i];
            // Add spaced kmer
            // There is your smer!!!
            smerlist[*smer_idx].key = malloc(info.smerlength + 1);
            memcpy(smerlist[*smer_idx].key, info.smer_seq, (info.smerlength + 1));
            *smer_idx += 1;
        }
    }
}


unsigned long belly_extract_spaces(kseq_t *seq,
                                   jellydata *jdata,
                                   jellyhash *smerhash,
                                   SpKMER *smerlist,
                                   jellyopts opts)
{
    unsigned long int total_seq = 0;
    int sequence_length;
    unsigned long int n = 0;
    unsigned long int total_kmers = 0;
    unsigned int size = (jdata->hashsize)*(jdata->buffersize);
    void *output;
    if (opts.raw) {
        output = jdata->routput;
    }
    else {
        output = jdata->soutput;
    }
    //Index for scaled vector
    unsigned long int v_idx = 0;
    unsigned long int numkmers = 0;

    fprintf(stderr, "\tINFO: Vectorizing sequence data.\n");
    //Read sequences in loop
    while ((sequence_length = kseq_read(seq)) >= 0) {
        //Several sanity checks
        //Check for valid sequence file format (fasta/fastq)
        if (sequence_length < 0) {
            fprintf(stderr,"\tERROR: Could not read in sequence file.\n");
            total_seq = 0;
            goto exit;
        }
        //Check for sequence >= kmerlength
        else if (sequence_length < jdata->kmerlength) continue;
        //Check for empty sequence
        else if (sequence_length == 0) continue;

        //Check for valid kmers (Only containing A C G T)
        if (!(numkmers = belly_count(seq->seq.s, sequence_length, jdata, smerhash))) {
            fprintf(stderr, "\tWARNING: No valid kmers in sequence: %s.\n", seq->name.s);
            n++;
            continue;
        }

        total_seq += sequence_length;
        total_kmers += numkmers;

        //Output spaced kmer vector for current sequence
        if (!opts.gmode) {
            //Move sequence to output buffer
            if (!belly_loadbuff(smerhash, opts.raw, v_idx, smerlist, output)) {
                fprintf(stderr, "\tERROR: Failed to load output buffer.\n");
                total_seq = 0;
                goto exit;
            }
            v_idx += jdata->hashsize;
            hash_reset(smerhash->h, smerhash->k);
            //If buffer is full, dump to output file
            if (!((n+1)%jdata->buffersize)) {
                fprintf(stderr, "\tINFO: Dumping vectors.\n");
                if (!belly_dump(jdata, opts, v_idx/(jdata->hashsize), size)) {
                    fprintf(stderr, "\tERROR: Failed to write output vectors.\n");
                    total_seq = 0;
                    goto exit;
                }
                v_idx = 0;
            }
        }
        n++;
    }

    if (n == 0) {
        fprintf(stderr,"\tERROR: Sequence file empty\n");
        total_seq = 0;
        goto exit;
    }
    //TODO: Better control of output writing logic. Not understandable at all

    fprintf(stderr, "\tINFO: Dumping vectors.\n");
    if (opts.gmode) {
        if (!belly_loadbuff(smerhash, opts.raw, v_idx, smerlist, output)) {
            fprintf(stderr, "\tERROR: Failed to load output buffer.\n");
            total_seq = 0;
            goto exit;
        }
        if (!belly_dump(jdata, opts, 1, kh_size(smerhash->h))) {
            fprintf(stderr, "\tERROR: Failed to write output vectors.\n");
            total_seq = 0;
            goto exit;
        }
    }
    else {
        if (!belly_dump(jdata, opts, v_idx/(jdata->hashsize), v_idx)) {
            fprintf(stderr, "\tERROR: Failed to write output vectors.\n");
            total_seq = 0;
            goto exit;
        }
    }

    if (opts.binout) {
        if (!belly_ofiletail(jdata->output, &n)) {
            fprintf(stderr, "\tERROR: Failed to finish output file.\n");
            return 0;
        }
    }

    //Log info
    fprintf(stderr, "\tINFO: Finished vectorizing sequence data.\n");
    fprintf(stderr,"\t\t%lu total reads.\n",n);
    fprintf(stderr, "\t\t%lu total kmers\n", total_kmers);
    fprintf(stderr, "\t\t%lu total bases\n", total_seq);
    exit:
        return(total_seq);
}


unsigned long int belly_count(char *seq,
                              int l,
                              jellydata *jdata,
                              jellyhash *smerhash)
{
    unsigned long int kmer_num = 0;
    unsigned long int valid_kmer_num = 0;
    //Loop over possible kmers in sequence
    for (; kmer_num < l - jdata->kmerlength + 1; kmer_num++) {
        //Get kmer from sequence
        //TODO Instead of copying kmer sequence to new array, pass address directly to
        //hash_genSpacemer.
        jdata->kmer_seq = memcpy(jdata->kmer_seq, &seq[kmer_num], jdata->kmerlength);
        //Extract positions stored in mask
        hash_genSpacemer(jdata->kmer_seq, jdata->mask, jdata->smer_seq);
        //Check if spaced kmer already exists and increase number
        smerhash->k = kh_get(smer, smerhash->h, jdata->smer_seq);
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


int belly_minmax(jellyhash *smerhash,
                 unsigned int *min,
                 unsigned int *diff,
                 SpKMER *smerlist)
{
    *min = -1;
    unsigned int max = 0;
    unsigned int val;
    for (unsigned long int i = 0; i < kh_size(smerhash->h); i++) {
        smerhash->k = kh_get(smer, smerhash->h, smerlist[i].key);
        val = kh_val(smerhash->h, smerhash->k);
        if (val < *min) *min  = val;
        if (val > max) max = val;
    }
    //TODO Check max > min
    *diff = max - *min;
    if (*diff == 0) {
        *diff = 1;
        *min = 0;
    }
    return 1;
}


void belly_exit(gzFile fp,
                kseq_t *seq,
                jellyhash smerhash,
                jellydata *jdata,
                SpKMER *smerlist)
{
  // Close zlib file handler
  gzclose(fp);
  // Close kseq pointer
  kseq_destroy(seq);
  // Free hash
  kh_destroy(smer, smerhash.h);
  // free spaced kmer arrays
  for(int i = 0; i < smerhash.hashsize; i++) {
    free(smerlist[i].key);
  }
  free(smerlist);
  free(jdata->kmer_seq);
  free(jdata->smer_seq);
  free(jdata->mask);
  free(jdata->soutput);
  free(jdata->routput);
  fclose(jdata->output);
}


int belly_loadbuff(jellyhash *smerhash,
                   int raw,
                   unsigned long int v_idx,
                   SpKMER *smerlist,
                   void *buff)
{
    if (raw) {
        unsigned int *tmp = buff;
        for (unsigned long int i = 0; i < kh_size(smerhash->h); i++) {
            smerhash->k = kh_get(smer, smerhash->h, smerlist[i].key);
            tmp[i + v_idx] = kh_val(smerhash->h, smerhash->k);
        }
    }
    else {
        float *tmp = buff;
        unsigned int min;
        unsigned int diff;
        unsigned int val;
        belly_minmax(smerhash, &min, &diff, smerlist);
        for (unsigned long int i = 0; i < kh_size(smerhash->h); i++) {
            smerhash->k = kh_get(smer, smerhash->h, smerlist[i].key);
            val = kh_val(smerhash->h, smerhash->k);
            tmp[i + v_idx] =  ((float)val - (float)min) / (float)diff;
        }
    }
    return 1;
}


int belly_dump(jellydata *jdata,
               jellyopts opts,
               unsigned long int nsamples,
               unsigned int size)
{
    if (opts.binout) {
        if (opts.raw) {
            fwrite(jdata->routput,
                   sizeof(unsigned int),
                   size,
                   jdata->output);
        }
        else {
            int i = fwrite(jdata->soutput,
                   sizeof(float),
                   size,
                   jdata->output);
        }
        //TODO check for fwrite errors
        return 1;
    }

    if (opts.raw) {
        belly_writeraw(jdata, nsamples);
    }
    else {
        belly_writescale(jdata, nsamples);
    }
        //}
    return 1;
}


int randint(int max)
{
    return 0 + rand() / (RAND_MAX / (max - 0 + 1) + 1);
}


int belly_checkfile(FILE *fp)
{
    //Check for output file
    if (!fp) {
        fprintf(stderr, "\tERROR: Could not open output file.\n");
        if (errno == 13 || errno == 21) {
            fprintf(stderr, "\t\tPermission denied %d\n", errno);
        }
        else if (errno == 21) {
            fprintf(stderr, "\t\tFilename is a directory %d\n", errno);
        }
        return 0;
    }
    return 1;
}


int belly_ofileinit(jellydata *jdata)
{
    char *head_start[] = {0,0,0,0,0,0,0};
    char *spacer[] = {0,0,0};
    char tail[] = {5,4,3,2,1};
    int bytenum = 0;
    bytenum += fwrite(head_start, 1, 7, jdata->output);
    bytenum += fwrite(&(jdata->smerlength), sizeof(int), 1, jdata->output);
    bytenum += fwrite(spacer, 1, 3, jdata->output);
    bytenum += fwrite(&(jdata->kmerlength), sizeof(int), 1, jdata->output);
    bytenum += fwrite(spacer, 1, 3, jdata->output);
    bytenum += fwrite(&tail, 1, 5, jdata->output);
    if (bytenum != 20) {
        return 0;
    }
    return 1;
}


int belly_ofiletail(FILE *fp, unsigned long int *n)
{
    char *spacer[] = {0,0,0,0};
    char tail[] = {0,0,0,0,1};
    int bytenum = 0;
    bytenum += fwrite(spacer, 1, 4, fp);
    bytenum += fwrite(n, sizeof(unsigned long int), 1, fp);
    bytenum += fwrite(&tail, 1, 5, fp);
    if (bytenum != 10) {
        return 0;
    }
    return 1;
}


int belly_writescale(jellydata *jdata, unsigned long nsamples)
{
    int idx = 0;
    for (size_t i = 0; i < nsamples; i++) {
        fprintf(jdata->output, "%.6f", jdata->soutput[i+idx]);
        for (size_t j = 1; j < jdata->hashsize; j++) {
            fprintf(jdata->output, "\t%.6f", jdata->soutput[j+idx]);
        }
        fprintf(jdata->output, "\n");
        idx += jdata->hashsize;
    }
}


int belly_writeraw(jellydata *jdata, unsigned long nsamples)
{
    int idx = 0;
    for (size_t i = 0; i < nsamples; i++) {
        fprintf(jdata->output, "%u", jdata->routput[i+idx]);
        for (size_t j = 1; j < jdata->hashsize; j++) {
            fprintf(jdata->output, "\t%u", jdata->routput[j+idx]);
        }
        fprintf(jdata->output, "\n");
        idx += jdata->hashsize;
    }
}
