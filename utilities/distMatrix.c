#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define HEADER  "#SampleID"

typedef struct {
  char sample[99];
  double *vec;
} sampleVec;


int readVars(char *line) {
  char *var;
  int numVars = 0;
  var = strtok(NULL, "\t");
  while(var != NULL) {
    numVars += 1;
    var = strtok(NULL, "\t");
  }
  return(numVars);
}


void readVec(char *line, sampleVec *sample, int numVars) {
  char *name = strtok(line,"\t");
  strcpy(sample->sample, name);

  sample->vec = malloc(numVars*sizeof(double));
  if(!sample->vec) {
    fprintf(stderr,"failed allocation\n");
  }

  int sampleVars = 0;
  char *var;
  var = strtok(NULL, "\t");
  while(var != NULL) {
    sample->vec[sampleVars] = atof(var);
    sampleVars += 1;
    var = strtok(NULL, "\t");
  }

  if(sampleVars != numVars) {
    fprintf(stderr,"Error in sample %s unmatched number of variables.\n",name);
  }
}

double dist(double *vec1, double *vec2, int len) {
  double S = 0;
  int i;
  for(i = 0; i < len; i++) {
    S += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);
  }
  return(sqrt(S));
}

int main(int argc, char *argv[]) {
  //Open file
  FILE *fp;
  fp = fopen(argv[1],"r");
  if(!fp) {
    fprintf(stderr,"Error opening file.\n");
    exit(1);
  }

  //Store sample names and vector values
  sampleVec *samples = malloc(300*sizeof(sampleVec));

  //Read lines
  char *line = NULL;
  size_t n = 0;
  char *sample;
  int numVars;
  getline(&line, &n, fp);
  sample = strtok(line,"\t");
  if(strcmp(sample,HEADER) != 0) {
    fprintf(stderr,"Error in matrix file. SampleID missing\n");
    exit(1);
  }
  else {
    numVars = readVars(line);
    printf("%d total variables\n",numVars);
  }
  int numSamples = 0;
  while(getline(&line, &n, fp) != -1){
    readVec(line, &samples[numSamples], numVars);
    numSamples += 1;
  }
  double S;
  //Naive printing
  int i;
  int j;
  for(i = 0; i < numSamples; i++){
    for(j = 0; j < numSamples; j++) {
      S = dist(samples[i].vec, samples[j].vec, numVars);
      printf("%f",S);
      if(j == (numSamples - 1)) {
          printf("\n");
        }
      else {
        printf("\t");
      }
    }
  }
  return(0);
}

