#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "khash.h"

KHASH_MAP_INIT_STR(single, unsigned int)
//unsigned char a[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//            0,0,0,0,0,2,0,4,0,0,0,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

//void algo3(unsigned char *idxmask, int len);
void algo3(unsigned char *idxmask, int len, unsigned char *a);
void fill3(char *smers, int len, unsigned long int num, khash_t(single) *h);

int main(int argc, char *argv[]) {
    khash_t(single) *h = kh_init(single);
    int len = atoi(argv[1]);
    unsigned long int num = pow(4, len);
    char *smers = malloc(num * (len + 1));
    if (!smers) {
        fprintf(stderr, "ERROR\n");
        exit(-1);
    }

    fprintf(stderr, "Allocated %ld bytes.\n", num * (len + 1));
    fill3(smers, len + 1, num*(len+1), h);
    fprintf(stderr, "hash size: %u\n", kh_size(h));
    kh_destroy(single, h);
    free(smers);
}

void algo3(unsigned char *idxmask, int len, unsigned char *a)
{
    for (int idx = len - 1; idx >= 0; idx--) {
        if (idxmask[idx] == 84) {
            idxmask[idx] = 65;
        }
        else {
            idxmask[idx] += a[idxmask[idx]];
            break;
        }
    }
}


void fill3(char *smers, int len, unsigned long int num, khash_t(single) *h)
{
    int absent = 0;
    khint_t k;
    unsigned char *mask = malloc(len);
    memset(mask, 65, len-1);
    mask[len - 1] = '\0';
    static unsigned char a[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                              0,0,0,0,0,2,0,4,0,0,0,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    fprintf(stderr, "Filling %lu.\n", num);
    for (unsigned long int i = 0; i < num; i+=len) {
        memcpy(smers + i, mask, len);
        k = kh_put(single, h, smers + i, &absent);
        kh_value(h, k) = 0;
        algo3(mask, len - 1, a);
    }
    fprintf(stderr, "Filled.\n");
    free(mask);
}
