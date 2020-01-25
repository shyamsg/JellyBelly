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
void fill2(char *smers, int smerlen, unsigned long int num, khash_t(single) *h);

int main(int argc, char *argv[])
{
  khash_t(single) *h = kh_init(single);
  int smerlength = atoi(argv[1]);
  unsigned long int numsmers = pow(4, smerlength);
  char *smers = malloc(numsmers * (smerlength + 1));
  if (!smers) {
    fprintf(stderr, "ERROR\n");
    exit(-1);
  }
  fprintf(stderr, "Allocated %ld bytes.\n", numsmers * (smerlength + 1));
  fill2(smers, smerlength + 1, numsmers * (smerlength + 1), h);
  fprintf(stderr, "hash size: %u\n", kh_size(h));
  kh_destroy(single, h);
  free(smers);
}



void algo2(int *idxmask, int len)
{
    for (int idx = len - 1; idx >= 0; idx--) {
        if (idxmask[idx] == 3) {
            idxmask[idx] = 0;
        }
        else {
            idxmask[idx] += 1;
            break;
            }
    }
}



void fill2(char *smers, int len, unsigned long int num, khash_t(single) *h)
{
    int absent = 0;
    khint_t k;
    int *smermask = malloc((len-1)*sizeof(int));
    memset(smermask, 0, len-1);

    //int idx = 0;
    int j = 0;

    fprintf(stderr, "Filling %lu.\n", num);
    for (unsigned long int i = 0; i < num; i+=len) {
        while (j < len - 1) {
            smers[i+j] = BASES[smermask[j]];
            j++;
        }
        j = 0;

        k = kh_put(single, h, smers + i, &absent);
        kh_value(h, k) = 0;
        algo2(smermask, len - 1);
    }
    fprintf(stderr, "Filled.\n");
    free(smermask);
}


/*
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

void printmask(int *smermask, int len)
{
    for (int i = 0; i < len; i++) {
    fprintf(stderr, "%d ", smermask[i]);
    }
    fprintf(stderr, "\n");
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
    //int maskid = 0;
    unsigned int smeridx = 0;
    for (unsigned long int i = 0; i < num; i++) {
        fprintf(stderr, "%lu %d ", i, idx);
        if (idx == smerlen) {
            smers[i] = '\0';
            //setidx(smermask, smerlen);
            algo2(smermask, smerlen);
            //Insert into hash
            //k = kh_put(single, h, smers + (smrcount * (smerlen + 1)), &absent);
            k = kh_put(single, h, smers + smeridx, &absent);
            kh_value(h, k) = 0;
            smrcount += 1;
            idx = 0;
            //maskid = 0;
            smeridx += (smerlen + 1);
            fprintf(stderr, "\n");
            continue;
        }
        //smers[i] = BASES[smermask[i - ( smrcount * (smerlen + 1) ) ] ];
        smers[i] = BASES[smermask[idx]];
        idx += 1;
        //maskid += 1;
    }
    free(smermask);
}

*/
