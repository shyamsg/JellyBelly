#include <stdio.h>
#include <stdlib.h>
#include "khash.h"
#define LETTERS "ABCDE"

KHASH_MAP_INIT_STR(test, unsigned int)

int main(){
  khint_t k;
  int absent;
  khash_t(test) *h;
  h = kh_init(test);

  char *a = malloc(25);
  int b = 0;
  for (int i = 0; i < 25; i+=5) {
      for (int j = 0; j < 4; j++) {
          a[j+i] = LETTERS[i-b];
      }
      b += 4;
      a[i+4] = '\0';
      kh_put(test, h, a+i, &absent);
  }

  for (k = kh_begin(h); k != kh_end(h); ++k)  // traverse
      if (kh_exist(h, k)) {            // test if a bucket contains data
          fprintf(stderr, "| %s\n",kh_key(h, k));
          kh_val(h, k) = 1;
          fprintf(stderr, "val: %u\n", kh_val(h, k));
      }

  free(a);
  kh_destroy(test, h);
}
