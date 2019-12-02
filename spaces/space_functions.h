#include "../klib/ksort.h"
#include "../klib/khash.h"


/* Structure holding for any given subkmer,
   it's mask, entropy value and subkmer at
   which entropy was computed */
typedef struct {
  char *mask;
  double ent;
} spacedMer;

//Initialize Heng Li's Ksort
#define pair_lt(a, b) (((a).ent < (b).ent))
#define pair_lt2(a, b) ((a) < (b))
KSORT_INIT(pair, spacedMer, pair_lt)
KSORT_INIT(pair2, float, pair_lt2)
KSORT_INIT_GENERIC(float)

//Initialize khash
KHASH_MAP_INIT_STR(entHash, unsigned int)

//Computes greatest comon denominator
unsigned long long mchoosen(int n, int k);


//Computes how many ways can n object be chosen from a set of m objects m > n
unsigned long long gcd(unsigned long long x, unsigned long long y);


//Computes all posible spaced kmers. This is not the actaul mask but the positions where mask has ones
unsigned long int computeSmers(int *binSpace,
                               int *space,
                               int smerlen,
                               int kmerlen);


//Given a kmer lenth and a spaced kmer length, move ones in mask
void moveOnes(int *space, int idx, int kmerlen, int smerlen);


//Based on spaced kmer positions create mask
void createMask(int *space, int kmerlen, int smerlen, char *smer);


//Find where a mask should be '1'
int findPos(int *space, int val, int smerlen);


//Compute entropy of a spaced kmer
double computeEnt(char *smer,
                  int kmerlen,
                  khash_t(entHash) *h,
                  khint_t k);


//free memory allocated to list
void listErraser(spacedMer *spaceList,long tracer);


//Compute entropy based on a subkmer number
double computeEnt2(char *smer,
                   int subkmerNum,
                   khash_t(entHash) *h,
                   khint_t k,
                   int kmerlen);


//Get index of spaced kmer with max entropy
long listReducer(spacedMer *spaceList, long length);


//Computes entropy of spaced kmers and keeps maximally entropic ones
void entropyCompute(int * binSapce,
                    char *space_list,
                    float *ent_list,
                    int kmerlen,
                    int smerlen,
                    unsigned long long smernum);


//Filter spaced kmers to keep maximally entropic ones
int filterSmers(float *ent_list,
                long array_length);


//Memoty error in case of failed allocations
void memerr();


//Free memory
void freeMem(char *space_list,
             float *ent_list,
             char *mask,
             khash_t(entHash) *h,
             khint_t k);


//Write spaced kmer masks to file
void writeSPacedKmers(char *space_list,
                      long arrayLength,
                      int kmerlen,
                      int smerlen,
                      unsigned long int smernum);


int write_header(FILE *bitmask_file, int kmerlen, int smerlen, unsigned long int smernum);

int write_footer(FILE *bitmask_file);
