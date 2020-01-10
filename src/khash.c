#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "klib/khash.h"

unsigned char seq_table[256] = {
	  0, 1,	2, 3, 4, 5,	6, 7, 8, 9, 10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 3, 66, 2, 68, 69, 70, 1, 72, 73, 74, 75, 76, 77, 4, 79,
	80, 81, 82, 83, 0, 85, 86, 87, 88, 89, 90,	91,	 92,  93,  94,	95,
	 96, 3, 98, 2, 100, 101, 102, 1, 104, 105, 106, 107, 108, 109, 4, 111,
	112, 113, 114, 115, 0, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};

unsigned char compl_table[256] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 66, 'G', 68, 69, 70, 'C', 72, 73, 74, 75, 76, 77, 'N', 79,
	80, 81, 82, 83, 'A', 85, 86, 87, 88, 89, 90,	91,	 92,  93,  94,	95,
	 96, 't', 98, 'g', 100, 101, 102, 'c', 104, 105, 106, 107, 108, 109, 110, 111,
	112, 113, 114, 115, 'a', 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};

KHASH_MAP_INIT_STR(smers, unsigned int)

typedef struct {
    char *key;
} SpKMER;

char *revcomp(int len, char *seq);
int compSeq(char *seq, char *revcomp, int len);


int main()
{
    khash_t(smers) *h;
    khint_t k, revk;
    int absent;
    clock_t start, end;
    double cpu_time_used;
    h = kh_init(smers);
    SpKMER *mers = malloc(5*sizeof(SpKMER));
    mers[0].key = "ACGT";
    mers[1].key = "GACG";
    mers[2].key = "CGTC";
    mers[3].key = "GTAG";
    mers[4].key = "AAAA";

    char *revseq;
    for (int i = 0; i< 5; i++) {
        fprintf(stderr, "%s\t%s\n", mers[i].key, revseq = revcomp(4,mers[i].key));
        fprintf(stderr, "%d\n", compSeq(mers[i].key, revseq, 4));
        fprintf(stderr, "%d\n", strncmp(mers[i].key, revseq, 4));
        k = kh_put(smers, h, mers[i].key, &absent);
        kh_value(h, k) = 10;
        free(revseq);
    }

    fprintf(stderr,"printing values and deleting keys\n");
    unsigned int keyval, revkeyval;
    for (int i = 0; i< 5; i++) {
        k = kh_get(smers, h, mers[i].key);
        if (k != kh_end(h)) {
            keyval = kh_val(h, k);
            revseq = revcomp(4,mers[i].key);
            if (strcmp(mers[i].key,revseq)) {
                revk = kh_get(smers, h, revseq);
                if (revk != kh_end(h)) {
                    revkeyval = kh_val(h, revk);
                    kh_val(h, k) = keyval + revkeyval;
                    kh_del(smers, h, revk);
                }
            }
            else {
                //printf("%s\t%d\t%s\n",kh_key(h,k),kh_val(h,k), revseq);
            }
            free(revseq);
        }
    }
    for (int i = 0; i< 5; i++) {
        k = kh_get(smers, h, mers[i].key);
        if (k != kh_end(h)) {
            printf("%s\t%d\n", kh_key(h,k), kh_val(h, k));
        }
    }


    free(mers);
    kh_destroy(smers,h);

    start = clock();
    int a;
    for (long int i = 0; i < 60000000; i++) {
        a = compSeq("ACGTACGT", "CCGTACGT", 8);
    }
    end = clock();
    fprintf(stderr, "%d\n",a);
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    fprintf(stderr, "compSeq took: %.10f seconds.\n", cpu_time_used);
    start = clock();
    for (long int i = 0; i < 60000000; i++) {
        a = strcmp("ACGTACGT", "CCGTACGT");
    }
    end = clock();
    fprintf(stderr, "%d\n",a);
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    fprintf(stderr, "strcmp took: %.10f seconds.\n", cpu_time_used);
}


char *revcomp(int len, char *seq) {
    int i;
    char *rev;
    rev = malloc(len + 1);
    rev[len] = 0;
    for (i = 0; i < len; ++i)
        rev[len - i - 1] = compl_table[(int)seq[i]];
    return rev;
}

int compSeq(char *seq, char *revcomp, int len) {
    for(int i = 0; i < len; i++) {
        //Continue if same base
        if(!(seq_table[(int)seq[i]]-seq_table[(int)revcomp[i]])) continue;
        //Return difference in ASCII value of base
        return seq_table[(int)seq[i]]-seq_table[(int)revcomp[i]];
    }
    return 0;
}
