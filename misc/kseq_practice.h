
typedef char* j_char;
//Declare klib's dq
KDQ_INIT(j_char)

//Declare klib's kseq
KSEQ_INIT(gzFile, gzread)

typedef struct {
    kseq_t *seq;
    unsigned long n;
    //double time;
    kdq_t(j_char) *q;
    //Syncronization data
    pthread_mutex_t jelly_mtx;
} jelly_t;

void *jelly_read(void *belly_t);
void *jelly_consume(void *data_t);
jelly_t *jelly_init(gzFile fp);
void jelly_kill(jelly_t *belly_t);
