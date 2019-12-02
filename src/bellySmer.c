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


void exitFree(khash_t(smer) *h,
              khint_t k,
              kseq_t *seq,
              gzFile fp,
              char *kmer,
              char *spacedKmer) {
  // Free memory alocated inside hash table
  hashfree(h, k);
  // Delete hash table as per khash
  kh_destroy(smer, h);
  // Free memory from kseq
  kseq_destroy(seq);
  // Close zlib file handler
  gzclose(fp);
  // Free pointers to spaced kmer and kmer
  free(kmer);
  free(spacedKmer);
}


int readMaskfromFile(int klen, int slen, char *filename, int *mask, FILE *smerFILE) {
  if (smerFILE = fopen(filename,"r")) {
    //gensmers(klen, slen, filename);
    return 0;
  }
  else {
    //gensmers(klen, slen, filename);
    return 1;
  }
}


//Compute greatest common denominator
unsigned long long gcd(unsigned long long x, unsigned long long y) {
  while (y != 0) {
    unsigned long long t = x % y;
    x = y;
    y = t;
  }
  return x;
}


//Computes how many ways can n object be chosen from a set of m objects m > n
unsigned long long mchoosen(int m, int n) {
  if (m < n) {
    printf("Error n should be smaller than m\n");
    exit(1);
  }
  unsigned long long r = 1;
  for (unsigned long long d = 1; d <= n; ++d, --m) {
    unsigned long long g = gcd(r, d);
    r /= g;
    unsigned long long t = m / (d / g);
    r *= t;
  }
  return r;
}

