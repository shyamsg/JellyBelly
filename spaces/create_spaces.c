#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void computeSpaces(int kmerlen, int smerlen);

void usage(){
  printf("spaces 0.0\nUsage:\nspaces <kmer length> <mask lenght>\n");
}

int main(int argc, char *argv[]) {

  //Check for correct number of parameters
  if(argc < 3) {
    usage();
    exit(-1);
  }

  //Read data input
  char *p;
  int kmerlen = (int)strtol(argv[1],&p,10);
  if (*p != 0 || kmerlen == 0) {
    printf("Invalid input for kmer lenght\n");
    usage();
    exit(-1);
  }

  int smerlen = (int)strtol(argv[2],&p,10);
  if (*p != 0 || smerlen == 0) {
    printf("Invalid input for mask lenght\n");
    usage();
    exit(-1);
    }

  if(smerlen >= kmerlen) {
    printf("mask length must be less than kmer length\n");
    usage();
    exit(-1);
  }

  //Start mask computation
  computeSpaces(kmerlen - 2, smerlen - 2);

}
