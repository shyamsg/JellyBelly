#include <zlib.h>

typedef struct {
    char *seqfilename;
    char *smerfilename;
    int mode;
} JellyOpts;

void usage();

void read_opts(int argc, char **argv, JellyOpts *opts);

void print_opts(JellyOpts opts);

int belly_start(gzFile fp, FILE *smer_file, JellyOpts opts);
