#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "khash.h"

#define BASES "ACGT"

KHASH_MAP_INIT_STR(single, unsigned int)

void fill(char *smers, int smerlen, unsigned long int num, khash_t(single) *h);
void setidx(int *idxmask, int len);
void printmask(int *smermask, int len);
void algo2(int *idxmask, int len);


int main(int argc, char *argv[])
{
  khash_t(single) *h = kh_init(single);
  //khint_t k;

  int smerlength = atoi(argv[1]);
  unsigned long int numsmers = pow(4, smerlength);


  char *smers = malloc(numsmers * (smerlength + 1));
  if (!smers) {
    fprintf(stderr, "ERROR\n");
    exit(-1);
  }
  fprintf(stderr, "Allocated %ld bytes.\n", numsmers * (smerlength + 1));

  fill(smers, smerlength, numsmers * (smerlength + 1), h);
  fprintf(stderr, "Filled.\n");
  for (int i = 0; i < 60; i+=15) {
      fprintf(stderr, "smers[0]: %s\n", smers + i);
  }

  fprintf(stderr, "hash size: %u\n", kh_size(h));
  kh_destroy(single, h);
  free(smers);

  fprintf(stderr, "Mem free.\n");
}


void fill(char *smers, int smerlen, unsigned long int num, khash_t(single) *h)
{
    int absent = 0;
    khint_t k;
    int *smermask = malloc(smerlen*sizeof(int));
    memset(smermask, 0, smerlen);

    fprintf(stderr, "Filling %lu.\n", num);
    unsigned long int smrcount = 0;
    int idx = 0;
    int maskid = 0;
    for (unsigned long int i = 0; i < num; i++) {
        if (idx == smerlen) {
            smers[i] = '\0';
            //setidx(smermask, smerlen);
            algo2(smermask, smerlen);
            //Insert into hash
            //k = kh_put(single, h, smers + (smrcount * (smerlen + 1)), &absent);
            //kh_value(h, k) = 0;
            smrcount += 1;
            idx = 0;
            maskid = 0;
            //fprintf(stderr, "\n");
            continue;
        }
        //fprintf(stderr, "%ld |%d|", i - (smrcount*(smerlen + 1)), maskid);
        //smers[i] = BASES[smermask[i - ( smrcount * (smerlen + 1) ) ] ];
        smers[i] = BASES[smermask[maskid]];
        idx += 1;
        maskid += 1;
    }
    free(smermask);
}


void setidx(int *idxmask, int len)
{
  int flag = 0;

  if (idxmask[len - 1] == 4) {
    idxmask[len - 1] = 0;
    flag = 1;
  }
  else {
    idxmask[len - 1] += 1;
  }

  for (int idx = len - 1; idx >= 0; idx--) {
    if (flag) idxmask[idx] += 1;
    if (idxmask[idx] == 4) {
      idxmask[idx] = 0;
      flag = 1;
    }
    else flag = 0;

  }
}


void algo2(int *idxmask, int len)
{
    int flag = 1;
    for (int idx = len - 1; idx >= 0; idx--) {
        if (idxmask[idx] == 3) {
            idxmask[idx] = 0;
            flag = 1;
        }
        else if (flag) {
                idxmask[idx] += 1;
                flag = 0;
                break;
            }
    }
}

void printmask(int *smermask, int len)
{
  for (int i = 0; i < len; i++) {
    fprintf(stderr, "%d ", smermask[i]);
  }
  fprintf(stderr, "\n");
}
