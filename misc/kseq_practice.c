#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <pthread.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include "kseq.h"
#include "kdq.h"
#include "kseq_practice.h"


int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "ERROR: Not enough arguments.\n");
        exit(-1);
    }

    gzFile fp = gzopen(argv[1], "r");
    //unsigned long n;
    if (!fp) {
        fprintf(stderr, "ERROR: Failed to open input file - %s.\n", argv[1]);
        exit(-1);
    }

    //Start jelly_t struct
    jelly_t *belly_t = jelly_init(fp);
    fprintf(stderr, "Jelly object initialized.\n");

    //thread data creation
    pthread_t jelly_readTH;
    pthread_t jelly_consTH;
    pthread_t jelly_cons2TH;
    pthread_t jelly_cons3TH;
    pthread_attr_t attrTH;
    void *status;

    // Initialize and set thread detached attribute
    pthread_attr_init(&attrTH);
    pthread_attr_setdetachstate(&attrTH, PTHREAD_CREATE_JOINABLE);
    //thread creation
    int th_rc = pthread_create(&jelly_consTH, &attrTH, jelly_consume, (void *)belly_t);
    if (th_rc) {
        fprintf(stderr, "ERROR: Thread creation.\n");
    }
    th_rc = pthread_create(&jelly_cons2TH, &attrTH, jelly_consume, (void *)belly_t);
    if (th_rc) {
        fprintf(stderr, "ERROR: Thread creation.\n");
    }
    th_rc = pthread_create(&jelly_cons3TH, &attrTH, jelly_consume, (void *)belly_t);
    if (th_rc) {
        fprintf(stderr, "ERROR: Thread creation.\n");
    }
    th_rc = pthread_create(&jelly_readTH, &attrTH, jelly_read, (void *)belly_t);
    if (th_rc) {
        fprintf(stderr, "ERROR: Thread creation.\n");
    }

    // Free attribute and wait for the other threads
    pthread_attr_destroy(&attrTH);
    int j_rc = pthread_join(jelly_readTH, &status);
    if (j_rc) {
        fprintf(stderr, "Thread error join. Return code: %d\n", j_rc);
    }
    j_rc = pthread_join(jelly_consTH, &status);
    if (j_rc) {
        fprintf(stderr, "Thread error join. Return code: %d\n", j_rc);
    }
    j_rc = pthread_join(jelly_cons2TH, &status);
    if (j_rc) {
        fprintf(stderr, "Thread error join. Return code: %d\n", j_rc);
    }
    j_rc = pthread_join(jelly_cons3TH, &status);
    if (j_rc) {
        fprintf(stderr, "Thread error join. Return code: %d\n", j_rc);
    }

    jelly_kill(belly_t);
    gzclose(fp);
    fprintf(stderr, "Finished.\n");
    pthread_exit(NULL);
    exit(-1);
}

void *jelly_read(void *data_t)
{
    fprintf(stderr, "Reader thread started. %ld\n", syscall(SYS_gettid));
    jelly_t *belly_t = data_t;
    int s_len = 0;
    while ((s_len = kseq_read(belly_t->seq)) >= 0) {
        //Pointer should be freed when extracted from queue
        char *a = malloc(s_len + 1);
        if (!a) fprintf(stderr, "memerr\n");
        memcpy(a, belly_t->seq->seq.s, s_len + 1);
        belly_t->n += s_len;

        //Add to queue
        pthread_mutex_lock(&(belly_t->jelly_mtx));
        kdq_push(j_char, belly_t->q, a);
        pthread_mutex_unlock(&(belly_t->jelly_mtx));

    }

    pthread_mutex_lock(&(belly_t->jelly_mtx));
    kdq_push(j_char, belly_t->q, NULL);
    kdq_push(j_char, belly_t->q, NULL);
    kdq_push(j_char, belly_t->q, NULL);
    pthread_mutex_unlock(&(belly_t->jelly_mtx));
    fprintf(stderr, "%lu bases.\n", belly_t->n);
    pthread_exit(NULL);
}

void *jelly_consume(void *data_t)
{
    sleep(5);
    fprintf(stderr, "Consumer thread started %ld.\n", syscall(SYS_gettid));
    jelly_t *belly_t = data_t;
    char **buff = malloc(300000*sizeof(char *));
    if (!buff) {
        fprintf(stderr, "ERROR: memory.\n");
        pthread_exit(NULL);
    }
    int breakflag = 0;
    int size;
    while ( 1 ) {
        pthread_mutex_lock(&(belly_t->jelly_mtx));
        size = kdq_size(belly_t->q);
        if ( size > 300000) {
            size = 300000;
        }
        else if ( 0 == size ) {
            size = -1;
        }

        for (int i = 0; i < size; i++) {
            buff[i] = *kdq_shift(j_char, belly_t->q);
            if (!buff[i]) {
                fprintf(stderr, "End of queue. Exiting %ld\n", kdq_size(belly_t->q));
                breakflag = 1;
                size = i;
                break;
            }
        }
        pthread_mutex_unlock(&(belly_t->jelly_mtx));
        //work
        for (int i = 0; i < size; i++) {
            free(buff[i]);
        }
        if (breakflag) break;
    }
    free(buff);
    pthread_exit(NULL);
}

jelly_t *jelly_init(gzFile fp)
{
    static jelly_t belly_t;
    /* Using dynamic memory
    jelly_t *belly_t;
    belly_t = malloc(sizeof(jelly_t));
    if (!belly_t) {
        fprintf(stderr, "ERROR: Failed to init jelly_t (mem aloc err).\n");
        return NULL;
    }
    */
    belly_t.seq = kseq_init(fp);
    belly_t.n = 0;
    belly_t.q = kdq_init(j_char);
    //Initialize syncronization data
    pthread_mutex_init(&(belly_t.jelly_mtx), NULL);
    return &belly_t;
}

void jelly_kill(jelly_t *belly_t)
{
    kseq_destroy(belly_t->seq);
    kdq_destroy(j_char, belly_t->q);
    pthread_mutex_destroy(&(belly_t->jelly_mtx));
    //Using dynamic memory
    //free(belly_t);
}
