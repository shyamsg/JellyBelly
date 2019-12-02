#define BASES "ACGT"

void hash_declare(jellyhash *smerhash);

// Fill hash with all posible combinations of spaced kmers
//TODO Document
void hash_fill(char *str,
               int idx,            // index tracking str position
               int len,            // length of spaced kmer
               khash_t(smer) *h,   // hash table
               khint_t k);          // hash index

// free memory alocated in hash table
//TODO Document
//void hash_free(khash_t(smer) *h, khint_t k);
void hash_free(jellyhash smerhash);

//TODO Document
void hash_genSpacemer(const char *kmer, int *mask, char *smer);


//TODO Document
void hash_print(khash_t(smer) *h,
                khint_t k,
                SpKMER *smerList,
                unsigned int size);


// reset hash tabel to 0
void hash_reset(khash_t(smer) *h, khint_t k);

