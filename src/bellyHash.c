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
#include "khash.h"
#include "bellyFun.h"
#include "bellyHash.h"

void hash_declare(jellyhash *smerhash) {
    smerhash->h = kh_init(smer);
    smerhash->k = kh_end(smerhash->h);
}

/*
  Fill hash with all posible combinations of spaced kmers
*/
void hash_fill(char *str,
               int idx,            // index tracking str position
               int len,            // length of spaced kmer
               khash_t(smer) *h,   // hash table
               khint_t k) {        // hash index
    int absent = 0;
    int i;
    if (idx < (len - 1)) {
        for (i = 0; i < strlen(BASES); i++) {
            str[idx] = BASES[i];
            hash_fill(str, idx + 1, len, h, k);
        }
    }
    else {
        for (i = 0; i < strlen(BASES); i++) {
            str[idx] = BASES[i];
            // Add spaced kmer
            k = kh_put(smer, h, str, &absent);
            kh_key(h, k) = strdup(str);
            // Set value to 0
            kh_value(h, k) = 0;
        }
    }
}

/*
  free memory alocated in hash table
*/
//void hash_free(khash_t(smer) *h, khint_t k)
void hash_free(jellyhash smerhash)
{
    for (smerhash.k = 0; smerhash.k < kh_end(smerhash.h); ++(smerhash.k)) {
    if (kh_exist(smerhash.h, smerhash.k)) {
      free((char*)kh_key(smerhash.h, smerhash.k));
    }
  }
}

/*
  Generate all possible spaced kmers
  TODO - implement, not finished yet
*/
void hash_genSpacemer(const char *kmer, int *mask, char *smer) {
    for (; *mask != -1; ++mask) {
        *(smer++) = kmer[*mask];
    }
    //*smer = '\0';
}

// Prints contents of hash in key\tvalue format
void hash_print(khash_t(smer) *h,
                khint_t k,
                SpKMER *smerList,
                unsigned int size)
{
  //Loop over hash integers and print key value pair
  for (int i = 0; i < size; i++) {
    k = kh_get(smer,h,smerList[i].key);
    printf("%s\t%d\n",kh_key(h,k),kh_val(h,k));
  }
}

// Reset hash tabelo to zeros
void hash_reset(khash_t(smer) *h, khint_t k) {
  for (k = 0; k < kh_end(h); ++k) {
    if (kh_exist(h, k)) {
      kh_value(h, k) = 0;
    }
  }
}

