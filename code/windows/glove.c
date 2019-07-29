//  GloVe: Global Vectors for Word Representation
//
//  Copyright (c) 2014 The Board of Trustees of
//  The Leland Stanford Junior University. All Rights Reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//
//  For more information, bug reports, fixes, contact:
//    Jeffrey Pennington (jpennin@stanford.edu)
//    GlobalVectors@googlegroups.com
//    http://nlp.stanford.edu/projects/glove/


#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define _FILE_OFFSET_BITS 64
#define MAX_STRING_LENGTH 1000

typedef float real;

struct student{
 float data;
};
typedef struct cooccur_rec {
    int word1;
    int word2;
    real val;
} CREC;

int k = 4; //kmer length
int stop = 1; // 0 (run all iteration) or 1 (stop when cost no decrease )
int num_iter = 40; // Number of full passes through cooccurrence matrix
int vector_size = 30; // Word vector size
int use_binary = 0; // 0: save as text files; 1: save as binary; 2: both. For binary, save both word and context word vectors.
int model = 0; // For text file output only. 0: concatenate word and context vectors (and biases) i.e. save everything; 1: Just save word vectors (no bias); 2: Save (word + context word) vectors (no biases)
real eta = 0.01; // Initial learning rate
real *W, *gradsq;
long long vocab_size;
char *input_file, *save_W_file, *name_file;

/* Efficient string comparison */
int scmp( char *s1, char *s2 ) {
    while (*s1 != '\0' && *s1 == *s2) {s1++; s2++;}
    return(*s1 - *s2);
}

void initialize_parameters() {
	long long a, b;
	vector_size++; // Temporarily increment to allocate space for bias
    
	/* Allocate space for word vectors and context word vectors, and correspodning gradsq */
	W = (real*)_aligned_malloc(2 * vocab_size * (vector_size + 1) * sizeof(real), 128); // Might perform better than malloc
    if (W == NULL) {
        fprintf(stderr, "Error allocating memory for W\n");
        exit(1);
    }
	gradsq =(real*) _aligned_malloc(2 * vocab_size * (vector_size + 1) * sizeof(real), 128); // Might perform better than malloc
	if (gradsq == NULL) {
        fprintf(stderr, "Error allocating memory for gradsq\n");
        exit(1);
    }
    FILE *fout;
    fout = fopen("seed.txt","r");
    struct student s[3];
	for (b = 0; b < vector_size; b++) for (a = 0; a < 2 * vocab_size; a++) {
	 fscanf(fout,"%f",&s[0].data);
     W[a * vector_size + b]=s[0].data;
	//W[a * vector_size + b] = (rand() / (real)RAND_MAX - 0.5) / vector_size;fprintf(fout,"%f",W[a * vector_size + b] );
}
        fclose(fout);
	for (b = 0; b < vector_size; b++) for (a = 0; a < 2 * vocab_size; a++) gradsq[a * vector_size + b] = 1.0; // So initial value of eta is equal to initial learning rate
	vector_size--;
}

/* Train the GloVe model */
real glove_thread() {
    long long a, b ,l1, l2;
    CREC cr;
	real diff, fdiff, temp1, temp2, cost = 0;
    FILE *fin;
    fin = fopen(input_file, "rb");
    
    real* W_updates1 = (real*)malloc(vector_size * sizeof(real));
    real* W_updates2 = (real*)malloc(vector_size * sizeof(real));
    while(true) {
        fread(&cr, sizeof(CREC), 1, fin);
        if (feof(fin)) break;
        if (cr.val <= 1) continue;
        
			/* Get location of words in W & gradsq */
			l1 = ((long long)cr.word1) * (vector_size + 1); // cr word indices start at 1
			l2 = (((long long)cr.word2) + vocab_size) * (vector_size + 1); // shift by vocab_size to get separate vectors for context words
			
			/* Calculate cost, save diff for gradients */
			diff = 0;
			for (b = 0; b < vector_size; b++) diff += W[b + l1] * W[b + l2]; // dot product of word and context word vector
			diff += W[vector_size + l1] + W[vector_size + l2] - log(cr.val); // add separate bias for each word
			fdiff = diff; // multiply weighting function (f) with diff

			// Check for NaN and inf() in the diffs.
			if (isnan(diff) || isnan(fdiff) || isinf(diff) || isinf(fdiff)) {
				fprintf(stderr,"Caught NaN in diff for kdiff for thread. Skipping update");
				continue;
			}

			cost += 0.5 * fdiff * diff; // weighted squared error
			
			/* Adaptive gradient updates */
			fdiff *= eta; // for ease in calculating gradient
			real W_updates1_sum = 0;
			real W_updates2_sum = 0;
			for (b = 0; b < vector_size; b++) {
				// learning rate times gradient for word vectors
				temp1 = fdiff * W[b + l2];
				temp2 = fdiff * W[b + l1];
				// adaptive updates
				W_updates1[b] = temp1 / sqrt(gradsq[b + l1]);
				W_updates2[b] = temp2 / sqrt(gradsq[b + l2]);
				W_updates1_sum += W_updates1[b];
				W_updates2_sum += W_updates2[b];
				gradsq[b + l1] += temp1 * temp1;
				gradsq[b + l2] += temp2 * temp2;
			}
			if (!isnan(W_updates1_sum) && !isinf(W_updates1_sum) && !isnan(W_updates2_sum) && !isinf(W_updates2_sum)) {
				for (b = 0; b < vector_size; b++) {
					W[b + l1] -= W_updates1[b];
					W[b + l2] -= W_updates2[b];
				}
			}

			// updates for bias terms
			W[vector_size + l1] -= fdiff / sqrt(gradsq[vector_size + l1]);
			W[vector_size + l2] -= fdiff / sqrt(gradsq[vector_size + l2]);
			fdiff *= fdiff;
			gradsq[vector_size + l1] += fdiff;
			gradsq[vector_size + l2] += fdiff;
    }
    free(W_updates1);
    free(W_updates2);
    
    fclose(fin);
	return cost;
}

/* Save params to file */
int save_params(int nb_iter) {
    /*
     * nb_iter is the number of iteration (= a full pass through the cooccurrence matrix).
     *   nb_iter > 0 => checkpointing the intermediate parameters, so nb_iter is in the filename of output file.
     *   else        => saving the final paramters, so nb_iter is ignored.
     */

    long long a, b;
    char format[20];
    char output_file[MAX_STRING_LENGTH], output_file_gsq[MAX_STRING_LENGTH];
    char *word = (char*)malloc(sizeof(char) * MAX_STRING_LENGTH + 1);
    FILE *fid, *fout;
    
    if (use_binary > 0) { // Save parameters in binary file
        if (nb_iter <= 0)
            sprintf(output_file,"%s.bin",save_W_file);
        else
            sprintf(output_file,"%s.%03d.bin",save_W_file,nb_iter);

        fout = fopen(output_file,"wb");
        if (fout == NULL) {fprintf(stderr, "Unable to open file %s.\n",save_W_file); return 1;}
        for (a = 0; a < 2 * (long long)vocab_size * (vector_size + 1); a++) fwrite(&W[a], sizeof(real), 1,fout);
        fclose(fout);
    }
    if (use_binary != 1) { // Save parameters in text file
        if (nb_iter <= 0)
            sprintf(output_file,"%s.txt",save_W_file);
        else
            sprintf(output_file,"%s.%03d.txt",save_W_file,nb_iter);
       
        //fout = fopen(output_file,"a+");
		fout = fopen(output_file, "wb");
        if (fout == NULL) {fprintf(stderr, "Unable to open file %s.\n",save_W_file); return 1;}

		char * currIdx = strrchr(input_file, '/');
		if (currIdx) fprintf(fout, ">%s\n", currIdx + 1);
		else fprintf(fout, ">%s\n", input_file);
		
        for (a = 0; a < vocab_size; a++) {
			//fprintf(fout, "%lld",a);
            if (model == 0) { // Save all parameters (including bias)
				for (b = 0; b < (vector_size + 1); b++)
				{
					if (0 == b) fprintf(fout, "%lf", W[a * (vector_size + 1) + b]);
					else fprintf(fout, ",%lf", W[a * (vector_size + 1) + b]);
				}
				for (b = 0; b < (vector_size + 1); b++) fprintf(fout, ",%lf", W[(vocab_size + a) * (vector_size + 1) + b]);
            }
			if (model == 1) // Save only "word" vectors (without bias)
			{
				for (b = 0; b < vector_size; b++)
				{
					if (0 == b) fprintf(fout, "%lf", W[a * (vector_size + 1) + b]);
					else fprintf(fout, ",%lf", W[a * (vector_size + 1) + b]);
				}
			}
			if (model == 2) // Save "word + context word" vectors (without bias)
			{
				for (b = 0; b < vector_size; b++)
				{
					if (0 == b) fprintf(fout, "%lf", W[a * (vector_size + 1) + b] + W[(vocab_size + a) * (vector_size + 1) + b]);
					else fprintf(fout, ",%lf", W[a * (vector_size + 1) + b] + W[(vocab_size + a) * (vector_size + 1) + b]);
				}
			}
			fprintf(fout,"\n");
        }

        fclose(fout);
    }
    return 0;
}

/* Train model */
int train_glove() {
    long long a;
    int save_params_return_code;
    int b;
	real total_cost = 0, prev_cost = -1; 
	long long num_lines = 1e5;

    fprintf(stderr, "TRAINING MODEL\n");

    initialize_parameters();
    
    time_t rawtime;
    struct tm *info;
    char time_buffer[80];
    // Lock-free asynchronous SGD
    for (b = 0; b < num_iter; b++) {
		total_cost = glove_thread();

        time(&rawtime);
        info = localtime(&rawtime);
        strftime(time_buffer,80,"%x - %I:%M.%S%p", info);
        fprintf(stderr, "%s, iter: %03d, cost: %lf\n", time_buffer,  b+1, total_cost/num_lines);
		
		if (stop)
		{
			if (prev_cost < 0 || prev_cost >(total_cost / num_lines)) prev_cost = total_cost / num_lines;
			else break;
		}

    }
    return save_params(0);
}

int find_arg(char *str, int argc, char **argv) {
    int i;
    for (i = 1; i < argc; i++) {
        if (!scmp(str, argv[i])) {
            if (i == argc - 1) {
                printf("No argument given for %s\n", str);
                exit(1);
            }
            return i;
        }
    }
    return -1;
}

int main(int argc, char **argv) {
    int i;
    input_file = (char*)malloc(sizeof(char) * MAX_STRING_LENGTH);
    save_W_file =(char*) malloc(sizeof(char) * MAX_STRING_LENGTH);
	name_file = (char*)malloc(sizeof(char) * MAX_STRING_LENGTH);
    int result = 0;
    
    if (argc == 1) {
        printf("GloVe: Global Vectors for Word Representation\n");
        printf("Original Author: Jeffrey Pennington (jpennin@stanford.edu). Modified by Yang Lu (ylu465@usc.edu)\n\n");
        printf("Usage options:\n");
		printf("\t-k <int>\n");
		printf("\t\tKmer length; default 4\n");
        printf("\t-vector-size <int>\n");
        printf("\t\tDimension of word vector representations (excluding bias term); default 50\n");
        printf("\t-iter <int>\n");
        printf("\t\tNumber of training iterations; default 25\n");
        printf("\t-eta <float>\n");
        printf("\t\tInitial learning rate; default 0.05\n");
        printf("\t-binary <int>\n");
        printf("\t\tSave output in binary format (0: text, 1: binary, 2: both); default 0\n");
        printf("\t-model <int>\n");
        printf("\t\tModel for word vector output (for text output only); default 2\n");
        printf("\t\t   0: output all data, for both word and context word vectors, including bias terms\n");
        printf("\t\t   1: output word vectors, excluding bias terms\n");
        printf("\t\t   2: output word vectors + context word vectors, excluding bias terms\n");
        printf("\t-input-file <file>\n");
        printf("\t\tBinary input file of shuffled cooccurrence data (produced by 'cooccur' and 'shuffle'); \n");
        printf("\t-save-file <file>\n");
        printf("\t\tFilename, excluding extension, for word vector output; default vectors\n");
        result = 0;
    } else {
		if ((i = find_arg((char *)"-k", argc, argv)) > 0) k = atoi(argv[i + 1]);
		if ((i = find_arg((char *)"-stop", argc, argv)) > 0) stop = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-vector-size", argc, argv)) > 0) vector_size = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-iter", argc, argv)) > 0) num_iter = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-eta", argc, argv)) > 0) eta = atof(argv[i + 1]);
        if ((i = find_arg((char *)"-binary", argc, argv)) > 0) use_binary = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-model", argc, argv)) > 0) model = atoi(argv[i + 1]);
        if (model != 0 && model != 1) model = 2;
        if ((i = find_arg((char *)"-save-file", argc, argv)) > 0) strcpy(save_W_file, argv[i + 1]);
        else strcpy(save_W_file, (char *)"vectors");
        if ((i = find_arg((char *)"-input-file", argc, argv)) > 0) strcpy(input_file, argv[i + 1]);
		vocab_size = (long long)pow(4, k);

        result = train_glove();
    }
    free(input_file);
    free(save_W_file);
	free(name_file);
    return result;
}
