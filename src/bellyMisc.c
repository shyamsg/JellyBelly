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


#include <stdlib.h>
#include <stdio.h>

#include "bellyFun.h"
#include "khash.h"

//KHASH_MAP_INIT_STR(smer, unsigned int)

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


char seq_table[256] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 0, 66, 1, 68, 69, 70, 2, 72, 73, 74, 75, 76, 77, 4, 79,
	80, 81, 82, 83, 3, 85, 86, 87, 88, 89, 90,	91,	 92,  93,  94,	95,
	 96, 0, 98, 1, 100, 101, 102, 2, 104, 105, 106, 107, 108, 109, 110, 111,
	112, 113, 114, 115, 3, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};

//From Heng Li's minimap see file cmappy.h
void belly_revcomp(int len, char *seq, char *rev)
{
    //unsigned char *rev;
	  //rev = malloc(len + 1);
    //if (!rev) return NULL;
    for (int i = 0; i < len; i++)
  	    rev[len - i - 1] = compl_table[(int)seq[i]];
	  rev[len] = 0;
	  //return rev;
}


void belly_reducehash(jellyhash *smerhash,
                      SpKMER *smerlist,
                      unsigned int hashsize,
                      int smerlength)
{
    khint_t revk;
    unsigned int keyvalue, revkeyvalue;
    char *revcomp = malloc(smerlength + 1);
    //Loop over spaced kmer strings
    for (unsigned int i = 0; i < hashsize; i++) {
        //Extract hash value (k)
        smerhash->k = kh_get(smer, smerhash->h, smerlist[i].key);
        //Check if it still exists
        if (smerhash->k != kh_end(smerhash->h)) {
            //Obtain key value
            keyvalue = kh_val(smerhash->h, smerhash->k);
            //Obtain reverse compliment
            belly_revcomp(smerlength, smerlist[i].key, revcomp);
            //Chack for palindromic spaced kmer sequence
            if (strcmp(smerlist[i].key, revcomp)) {
                //Add values and rewrite to to original spaced kmer sequence
                revk = kh_get(smer, smerhash->h, revcomp);
                revkeyvalue = kh_val(smerhash->h, revk);
                kh_val(smerhash->h, smerhash->k) = keyvalue + revkeyvalue;
                //Delete reverse compliment's hash entry
                kh_del(smer, smerhash->h, revk);
            }
        }
    }
    free(revcomp);
    //Check for correct hash size (4**S)/2 if S is odd, (4**S)/2 + 4**S/2 is S is even
    /*
    unsigned int newsize =  kh_size(smerhash->h);
    fprintf(stderr, "Size is now %u\n", newsize);
    for (unsigned int i = 0; i < hashsize; i++) {
        smerhash->k = kh_get(smer, smerhash->h, smerlist[i].key);
        if (smerhash->k != kh_end(smerhash->h)) {
            fprintf(stderr,"%s\t%u\n",smerlist[i].key,kh_val(smerhash->h, smerhash->k));
        }
    }
    */
}


