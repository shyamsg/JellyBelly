#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "space_functions.h"


void computeSpaces(int kmerlen, int smerlen)
{
    fprintf(stderr,"Settin up\n%d %d\n",kmerlen, smerlen);

    //Array for storing mask positions
    int *space = malloc(smerlen * sizeof(int));
    if (!space) {
        memerr();
    }

    //Compute number of combinations of smerlen ones's in a kmer of kmer length
    //TODO make shure numbers check. min(kmerlen) 3 min(smerlen) 3
    unsigned long long smernum = mchoosen(kmerlen, smerlen);
    if(smernum > 300000000){
        fprintf(stderr,"Excesively large amount of space kmers.\nWill not compute.\n");
        exit(-1);
    }
    fprintf(stderr,"%lld spaced kmers will be generated\n",smernum);

    //Array of spaced kmers.
    int *binSpace = (int *)malloc(smernum * smerlen * sizeof(int));
    if (!binSpace) {
        memerr();
    }

    //Initialize array
    for (int i = 0; i < smerlen; i++) {
        space[i] = i;
    }

    //Miscelaneous counter
    unsigned long int counter = 0;
    //Compute spaced kmers and store number of spaced kmers generated
    counter  = computeSmers(binSpace, space, smerlen, kmerlen);
    if (counter + 1 != smernum) {
        fprintf(stderr,"Count discrepancy in spaced kmer index computation.\n");
        exit(1);
    }
    fprintf(stderr, "Finished: %lu spaced kmer indeces generated\n", counter+1);


    // Compute max entropy of spaced kmers based on subkmers
    // loop over subkmers
    //List of spaced kmers to keep depending on their entropy value
    //spacedMer *spaceList = malloc(10000*sizeof(spacedMer));
    char *space_list = (char *)malloc((10000*(kmerlen+2))*sizeof(char));
    float *ent_list = (float *)malloc(10000*sizeof(float));
    if (!space_list || !ent_list) { 
        memerr();
    }

    entropyCompute(binSpace,
                   space_list,
                   ent_list,
                   kmerlen,
                   smerlen,
                   smernum);

    free(space);
    free(binSpace);
    fprintf(stderr,"Done\n");
}


//Computes how many ways can n object be chosen from a set of m objects m > n
unsigned long long mchoosen(int n, int k)
{
  unsigned long long r = 1;
  for (unsigned long long d = 1; d <= k; ++d, --n) {
    unsigned long long g = gcd(r, d);
    r /= g;
    unsigned long long t = n / (d / g);
    r *= t;
  }
  return r;
}


//Computes greatest comon denominator
unsigned long long gcd(unsigned long long x, unsigned long long y)
{
  while (y != 0) {
    unsigned long long t = x % y;
    x = y;
    y = t;
  }
  return x;
}


//Computes all posible spaced kmers. This is not the actaul mask but the positions where mask has ones
unsigned long int computeSmers(int *binSpace, int *space, int smerlen, int kmerlen)
{

  unsigned long int counter = 0;
  //Compute all spaced kmers
  while(1) {
    // As single array is used, we add to pointer
    memcpy(binSpace + counter, space, smerlen*sizeof(int));
    moveOnes(space, smerlen - 1, kmerlen, smerlen);
    if (space[0] == kmerlen - smerlen) {
      counter += smerlen;
      memcpy(binSpace + counter, space, smerlen*sizeof(int));
      break;
    }
    counter += smerlen;
  }
  // Return right amount of spaced kmers hence index by ho size
  return counter/smerlen;
}


//Given a kmer lenth and a spaced kmer length, move ones in mask
void moveOnes(int *space, int idx, int kmerlen, int smerlen)
{
  if(space[0] == kmerlen - smerlen) {
    return;
  }
  if(space[idx] == kmerlen - smerlen + idx){
    int newidx = idx - 1;
    moveOnes(space, newidx, kmerlen, smerlen);
    space[idx] = space[newidx] + 1;
  }
  else {
    space[idx] += 1;
  }
}


//Create mask string
void createMask(int *space, int kmerlen, int smerlen, char *kmer_string)
{
    for(int i = 1; i < kmerlen + 1; i ++) {
        if(findPos(space, i - 1, smerlen)){
            kmer_string[i] = '1';
        }
        else {
            kmer_string[i] = '0';
        }
    }
}


//Find if kmer position is stored in space
int findPos(int *space, int val, int smerlen)
{
  for( int i = 0; i < smerlen; i ++) {
    if(space[i] == val) {
      return 1;
    }
  }
  return 0;
}


//Computes entropy of a spaced kmer
double computeEnt(char *smer,          //mask string
                  int kmerlen,         //kmer lenght == mask length
                  khash_t(entHash) *h, //hash table to store subkmers
                  khint_t k)           //hash int for hash table
{

  // Value to store entropy computation
  //double ent = 0.0;

  // Value storing last entropy computation. Based on a different subkmer length
  //double prevEnt = 0.0;

  double total_ent = 0.0;
  double ent = 0.0;
  // Loops over subkmer lengths starting at 2 and ending at kmerlen - 1
  for(int i = 2; i <= kmerlen - 1; i ++) {
      ent = computeEnt2(smer, i, h, k, kmerlen);
      total_ent += ent;

  }
  return total_ent;
}


//This is where the magic happens
//TODO DOCUMENT properly
double computeEnt2(char *smer,           //Array for spaced kmer mask string
                   int subkmerNum,       //subkmer length value
                   khash_t(entHash) *h,  //hash table
                   khint_t k,            //hash int
                   int kmerlen)          //kmer length
{
    int absent;
    //Number of subkmers in in kmer
    int numSubkmers = kmerlen - subkmerNum + 1;

    //Array for storing subkmer for entropy computation
    char *submer;
    submer = (char *)malloc(kmerlen*sizeof(char));

    //Fill hash table
    for(int j = 0; j <= kmerlen - subkmerNum; j++) {
      //Copy subkmer from mask into submer
      memcpy(submer,smer+j,subkmerNum*sizeof(char));
      //Add null terminator
      submer[subkmerNum] = '\0';

      // Add subkmer to hash table
      k = kh_get(entHash,h,submer);
      if (k == kh_end(h)) {
        k = kh_put(entHash, h, submer, &absent);
        //This causes memory leak if kh_key not freed at end of this loop.
        kh_key(h, k) = strdup(submer);
        // Set value to 1
        kh_value(h, k) = 1;
      }
      // Increase subkmer count if already added
      else {
        kh_val(h,k) += 1;
      }
    }


    /* Extract subkmers and add them to hash table for entropy computation*/
    double pent = 0.0;
    double ent = 0.0;
    //Loop over hash table and compute entropy
    for(k = 0;k < kh_end(h);k++) {
      if(kh_exist(h,k)) {
        pent = (double)kh_val(h,k)/numSubkmers;
        ent -= pent*log2(pent);
      }
    }
    //reset hash table and free alocated memory for key
    for(k = 0;k < kh_end(h);k++) {
      if (kh_exist(h, k)) {
        //Avoid memory leak cuased by strdup
        free((char*)kh_key(h, k));
      }
      kh_del(entHash,h,k);
    }
    free(submer);
    return ent;
}


void listErraser(spacedMer *spaceList, long length)
{
    for(int i = 0; i < length; i++) {
        free(spaceList[i].mask);
    }
}


//Compute index where list of spaced kmers decreases in entropy
long listReducer(spacedMer *spaceList, long length) {
    long maxIndex;
    float maxEnt = spaceList[length - 1].ent;

    for(maxIndex = length - 1; maxIndex > 0; maxIndex--){
        if(spaceList[maxIndex-1].ent < maxEnt) return maxIndex;
    }

    return 0;
}


//Computes entropy of spaced kmers and keeps maximally entropic ones
void entropyCompute(int *binSpace,
                    char *space_list,
                    float *ent_list,
                    int kmerlen,
                    int smerlen,
                    unsigned long long smernum) {

    //Array for storing mask
    //char *spacedKmer;
    char *mask = malloc((kmerlen + 2) * sizeof(char));;
    mask[0] = '1';
    mask[kmerlen + 1] = '1';


    int limit = 10000;

    //Entropy value
    double ent = 0.0;
    double prevEnt = 0.0;
    unsigned long int arrayLength = 0; //Tracks the number of spacedKmers in the array

    // Start elements for hash table
    khash_t(entHash) *h;
    h = kh_init(entHash);
    khint_t k;
    k = kh_end(h);

    //Loop over all spaced kmers compute entropy and keep all maximaly entropic
    unsigned long int i;
    unsigned long int arrayLimit = (smernum * smerlen) - smerlen;
    unsigned long int mask_index = 0;

    fprintf(stderr, "\nComputing entropy:\n\tLooping over spaced kmers\n");
    for(i = 0; i <= arrayLimit; i += smerlen) {

        //Log
        if ((i/smerlen)%10000 == 0) {
             fprintf(stderr,"\r\t%ld",i/(smerlen));
             fflush(stderr);
         }

        createMask(binSpace + i, kmerlen, smerlen, mask);

         //Compute entropy
        ent = computeEnt(mask, kmerlen + 2, h, k);

        //If a high entropy is encountered, errase list and continue
        if(ent > prevEnt) {
            //listErraser(spaceList,arrayLength);
            prevEnt = ent;
            limit = 10000;
            arrayLength = 0;
            mask_index = 0;
            memcpy(space_list + mask_index, mask, kmerlen + 2);
            ent_list[arrayLength] = (float)ent;
            arrayLength += 1;
            mask_index += kmerlen + 2;
        }
        //Keep spaced kmers within 0% of the max entropy spaced kmer computed so far
        else if(ent >= (prevEnt*1)) {
            if(arrayLength == limit) {
                limit += 10000;
                space_list = realloc(space_list,(limit*(kmerlen + 2))*sizeof(char));
                ent_list = realloc(ent_list, limit*sizeof(float));
            }
            ent_list[arrayLength] = (float)ent;
            memcpy(space_list + mask_index, mask, kmerlen + 2);
            arrayLength += 1;
            mask_index += kmerlen + 2;
        }
    }

    fprintf(stderr,"\n%ld spaced kmers kept\n",arrayLength);

    // Resize array to number of kept spaced kmers
    space_list = realloc(space_list,(arrayLength*(kmerlen + 2))*sizeof(char));
    ent_list = realloc(ent_list, arrayLength*sizeof(float));

    // Filter spaced kemrs to keep maximally entropic ones
    if (filterSmers(ent_list, arrayLength) != 0) {
        fprintf(stderr,"ERROR\n");
    }
    else fprintf(stderr,"PASS\n");


    writeSPacedKmers(space_list,
                     arrayLength,
                     kmerlen,
                     smerlen,
                     arrayLength);

    //Finish everything
    fprintf(stderr,"\nFreeing  mem\n");
    fprintf(stderr,"Before erasing size of array %ld\n",arrayLength);
    freeMem(space_list, ent_list, mask, h, k);
}


//Make shure all kmers have equal entropy values
int filterSmers(float *ent_list,
                long arrayLength) {
    ks_mergesort(pair2, arrayLength, ent_list, 0);
    if (ent_list[0] != ent_list[arrayLength - 1]) return -1;
    else {
        fprintf(stderr,"Entropy: %.10f\n", ent_list[0]);
        return 0;
    }

}

void memerr() {
  fprintf(stderr,"Error: Could not allocate enough memory.\n");
  exit(1);
}


void freeMem(char *space_list,
             float *ent_list,
             char *mask,
             khash_t(entHash) *h,
             khint_t k) {
  free(space_list);
  free(ent_list);
  free(mask);
  for (k = 0; k < kh_end(h); ++k) {
    if (kh_exist(h, k)) {
      free((char*)kh_key(h, k));
    }
  }
  kh_destroy(entHash, h);
}


void writeSPacedKmers(char *space_list,
                      long arrayLength,
                      int kmerlen,
                      int smerlen,
                      unsigned long smernum) {

    //Prepare filename
    char name[99] = "SpacedKmers";
    char skmerlen[99];
    char ssmerlen[99];
    sprintf(skmerlen,"%d", kmerlen + 2);
    sprintf(ssmerlen,"%d", smerlen + 2);
    strcat(name,"_K");
    strcat(name,skmerlen);
    strcat(name,"_S");
    strcat(name,ssmerlen);
    strcat(name,".bin");

    int written_bytes = 0;
    //Write metainformation to beginning of file
    FILE *bitmask_file;
    bitmask_file = fopen(name,"wb");
    if (!bitmask_file) {
        fprintf(stderr, "Error opening file.\n");
        exit(-1);
    }
    //40 byte header signature
    written_bytes += write_header(bitmask_file, kmerlen + 2, smerlen + 2, smernum);

    //write kmers
    unsigned long int amount;
    amount = (unsigned long int)arrayLength * (unsigned long int)(kmerlen +  2);
    written_bytes += fwrite(space_list, 1, amount, bitmask_file);

    //write 24 byte footer
    written_bytes += write_footer(bitmask_file);
    if (written_bytes != (48 + amount + 24)) {
        fprintf(stderr, "ERROR: Failed to write spaced kmer file\n");
        fprintf(stderr, "%d bytes wriyyen\n",written_bytes);
        fclose(bitmask_file);
        exit(-1);
    }
    //exit
    fclose(bitmask_file);
}


int write_header(FILE *bitmask_file, int kmerlen, int smerlen, unsigned long int smernum) {
    char *head_start[] = {0,0,0,0,0,0,0};
    char *spacer[] = {0,0,0};
    char *head_end[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};//,0,0,0,0,0,0,0,0,0,0,0};
    int byte_num = 0;
    int int_write;

    byte_num += fwrite(head_start, 1, 7, bitmask_file);

    int_write = fwrite(&kmerlen, sizeof(int), 1, bitmask_file);
    if (int_write == 1) {
        byte_num += 4;
    }

    byte_num += fwrite(spacer, 1, 3, bitmask_file);

    int_write = fwrite(&smerlen, sizeof(int), 1, bitmask_file);
    if (int_write == 1) {
        byte_num += 4;
    }

    byte_num += fwrite(spacer, 1, 3, bitmask_file);

    int_write = fwrite(&smernum, sizeof(unsigned long int), 1, bitmask_file);
    if (int_write == 1) {
        byte_num += 8;
    }

    byte_num += fwrite(head_end, 1, 19, bitmask_file);
    return byte_num;
}


int write_footer(FILE *bitmask_file) {
    char *foot[] = {0,0,0,0,0,0,0};
    char *signature = "SPACE0KMER0FILE:)";
    int byte_num  = 0;
    byte_num += fwrite(foot, 1, 7, bitmask_file);
    byte_num += fwrite(signature, 1, 17, bitmask_file);
    return byte_num;
}
