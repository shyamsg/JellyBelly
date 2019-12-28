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

#define BASES "ACGT"

void hash_declare(jellyhash *smerhash);

// Fill hash with all posible combinations of spaced kmers
/*
  Fill hash with all posible combinations of spaced kmers
*/
void hash_fill(char *str,
               int idx,            // index tracking str position
               int len,            // length of spaced kmer
               khash_t(smer) *h,   // hash table
               khint_t k);          // hash index

/*
  free memory alocated in hash table
*/
void hash_free(jellyhash smerhash);


/*
  Generate all possible spaced kmers
*/
void hash_genSpacemer(const char *kmer, int *mask, char *smer);


// Prints contents of hash in key\tvalue format
void hash_print(khash_t(smer) *h,
                khint_t k,
                SpKMER *smerList,
                unsigned int size);


// Reset hash tabelo to zeros
void hash_reset(khash_t(smer) *h, khint_t k);

