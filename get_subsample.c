/*
The MIT License (MIT)

Copyright (c) 2024 Hanfc <h2624366594@gmail.com>

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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <time.h>
#include <sys/stat.h>

#include "log.h"
#include "misc.h"
#include "seqtools.h"


// compare function for qsort
static int compare(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

// calculate the number of sequences in a fasta file
static void seq_counts_gz(gzFile gz_fp, int* seq_num) {
    *seq_num = 0;
    size_t bufsize = 1024;
    char* line = (char*)malloc(bufsize * sizeof(char));

    while (gzgets(gz_fp, line, bufsize) != NULL) {
        while (strlen(line) == bufsize - 1 && line[bufsize - 2] != '\n') {
            bufsize *= 2;
            line = (char*)realloc(line, bufsize * sizeof(char));
            gzgets(gz_fp, line + strlen(line), bufsize - strlen(line));
        }

        if (line[0] == '>') {
            (*seq_num)++;
        }
    }

    free(line);
}

static void seq_counts_ungz(FILE* fp, int* seq_num) {
    *seq_num = 0;
    char* line = NULL;
    size_t len = 0;

    while (getline(&line, &len, fp) != -1) {
        if (line[0] == '>') {
            (*seq_num)++;
        }
    }

    free(line);
}

// subsample a fasta file
void subsample(const char* output, const char* corrected_seq, double factor, int seed) {

    log_message(INFO, "Random select sequence start ...");

    uint64_t i;
    // Check if file exists
    FILE* check_file = fopen(output, "r");
    if (check_file != NULL) {
        log_message(WARNING, "'%s' already exists and rewrites ", output);
        fclose(check_file);
    }

    // log_message(INFO, "The assembly_seq %s", output);
    FILE* fp_result = fopen(output, "w");
    if (fp_result == NULL) {
        log_message(ERROR, "Failed to create output file");
        exit(EXIT_FAILURE);
    }

    gzFile gz_in = NULL;
    FILE* ungz_in = NULL;
    // Calculate the number of sequences in the fasta
    int seq_num = 0;

    if (validate_fasta_file(corrected_seq) == 1) {
        gz_in = gzopen(corrected_seq, "rt");
        if (gz_in == NULL) {
            log_message(ERROR, "Error opening gzipped file");
            fclose(fp_result);
            exit(EXIT_FAILURE);
        }
        seq_counts_gz(gz_in, &seq_num);
        gzclose(gz_in);
        gz_in = gzopen(corrected_seq, "rt");

    } else if (validate_fasta_file(corrected_seq) == 2) {
        ungz_in = fopen(corrected_seq, "r");
        if (ungz_in == NULL) {
            log_message(ERROR, "Error opening file");
            fclose(fp_result);
            exit(EXIT_FAILURE);
        }
        seq_counts_ungz(ungz_in, &seq_num);
        fclose(ungz_in);
        ungz_in = fopen(corrected_seq, "r");
    }  else {
        log_message(ERROR, "Invalid fasta file");
        fclose(fp_result);
        exit(EXIT_FAILURE);
    }

    if (factor < 1) {
        // log_message(INFO, "The factor is %g", factor);
        clock_t start = clock();

        int subset_size = (int)(seq_num * factor);
        if (seed) {
            srand(seed);
        } else {
            srand(time(NULL));
        }

        int* random_index = (int*)malloc(seq_num * sizeof(int));
        
        
        for (i = 0; i < seq_num; i++) {
            random_index[i] = i;
        }

        // Fisher-Yates shuffle algorithm

        for (i = 0; i < subset_size; i++) {
            int j = i + rand() % (seq_num - i);
            int temp = random_index[i];
            random_index[i] = random_index[j];
            random_index[j] = temp;
        }

        qsort(random_index, subset_size, sizeof(int), compare);

        int seq_count = 0;
        int index = 0;
        int flag = 0;
        if (gz_in) {
            size_t bufsize = 1024;
            char* line = (char*)malloc(bufsize * sizeof(char)); 
            while (gzgets(gz_in, line, bufsize) != NULL) {

                while (strlen(line) == bufsize - 1 && line[bufsize - 2] != '\n') {
                    bufsize *= 2;
                    line = (char*)realloc(line, bufsize * sizeof(char));
                    gzgets(gz_in, line + strlen(line), bufsize - strlen(line));
                }

                if (index > subset_size) {
                    break;
                } else {
                    if (line[0] == '>') {
                        if (seq_count == random_index[index]) {
                            fprintf(fp_result, "%s", line);
                            flag = 1;
                            index++;
                        } else {
                            flag = 0;
                        }
                        seq_count++;
                    } else if (flag == 1) {
                        fprintf(fp_result, "%s", line);
                    }
                }
            }

            free(line);
            gzclose(gz_in);
        } else {
            ssize_t seqline;
            char* line = NULL;
            size_t len = 0;
            while ((seqline = getline(&line, &len, ungz_in)) != -1) {
                if (index > subset_size) {
                    break;
                } else {
                    if (line[0] == '>') {
                        if (seq_count == random_index[index]) {
                            fprintf(fp_result, "%s", line);
                            flag = 1;
                            index++;
                        } else {
                            flag = 0;
                        }
                        seq_count++;
                    } else if (flag == 1) {
                        fprintf(fp_result, "%s", line);
                    }
                }
            }
            free(line);
            fclose(ungz_in);
        }

        free(random_index);
        fclose(fp_result);

        clock_t end = clock();
        double time_used = (double)(end - start) / CLOCKS_PER_SEC;
        // log_message(INFO, "Random select sequence end.");
        log_message(INFO, "Time used: %f s", time_used);
    }
}


#ifndef SUBSAMPLE_MAIN
// Function prototypes
void print_usage(const char *prog_name) {
    fprintf(stdout, "Usage:%s -i <fasta_file> -o <output_file> -f <factor> -s <seed>\n", prog_name);
    fprintf(stdout, "Required options:\n");
    fprintf(stdout, "   -i, --input     Input fasta file\n");
    fprintf(stdout, "   -f, --factor    Subsample factor (0~1)\n");
    fprintf(stdout, "   -o, --output    Output directory\n");

    fprintf(stdout, "Optional options:\n");
    fprintf(stdout, "   -s, --seed      Seed for random number generator, default is current time\n");
    fprintf(stdout, "   -h, --help      Display this help message\n");
}


void parse_arguments(int argc, char *argv[], char* input_file, char* output_file, double* factor, int* seed) {
    int i;
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            if (i + 1 < argc) {
                char *input_path = argv[++i];
                strcpy(input_file, input_path);
            }
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (i + 1 < argc) {
                // *prefix = argv[++i];
                char *output_path = argv[++i];
                strcpy(output_file, output_path);
            }
        } else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--factor") == 0) {
                *factor = atof(argv[++i]);
                if (*factor <= 0 || *factor >= 1) {
                    log_message(ERROR, "Invalid factor value");
                    exit(EXIT_FAILURE);
                }
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--seed") == 0) {
                *seed = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            exit(EXIT_SUCCESS);
        } else {
            log_message(ERROR, "Invalid option '%s'", argv[i]);
            print_usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (strlen(input_file) == 0 || strlen(output_file) == 0 || *factor == 0) {
        log_message(ERROR, "Please provide all required arguments");
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
}


int main(int argc, char* argv[]) {
    char input_file[256];
    char output_file[256];
    double factor;
    int seed = 0;

    parse_arguments(argc, argv, input_file, output_file, &factor, &seed);
    subsample(output_file, input_file, factor, seed);
    return 0;
}

#endif