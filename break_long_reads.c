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
#include <unistd.h>
#include <sys/types.h>
#include <zlib.h>

#include "log.h"
#include "misc.h"
#include "seqtools.h"


// Break long reads into shorter reads with uniform truncation
void BreakLongReads(const char *input_seq, const char *output_seq, int break_length) {
    log_message(INFO, "Reads breaking started...");

    uint64_t i;
    if (access(input_seq, F_OK) != 0) {
        log_message(ERROR, "File not found: %s", input_seq);
        exit(EXIT_FAILURE);
    }
    
    if (validate_fasta_file(input_seq) == 2) {
        FILE *inseq = fopen(input_seq, "r");
        FILE *outseq = fopen(output_seq, "w");

        if (inseq == NULL || outseq == NULL) {
            log_message(ERROR, "Failed to open file: %s", input_seq);
            exit(EXIT_FAILURE);
        }

        size_t len = 0;
        char *line = NULL;
        ssize_t read;

        int seq_count = 0;
        while ((read = getline(&line, &len, inseq)) != -1) {
            if (line[0] != '>') {
                line[strcspn(line, "\n")] = '\0';
                int read_length = strlen(line);
                
                if (read_length <= break_length) {
                    // If sequence length is less than or equal to break_length, write directly
                    seq_count++;
                    fprintf(outseq, ">%d\n", seq_count);
                    fprintf(outseq, "%s\n", line);
                } else {
                    // If sequence length exceeds break_length, perform uniform truncation
                    int min_segments = (read_length + break_length - 1) / break_length;
                    int longer_segments = read_length % min_segments;
                    int shorter_segments = min_segments - longer_segments;
                    int shorter_length = read_length / min_segments;
                    int longer_length = shorter_length + 1;
                    
                    int pos = 0;
                    for (i = 0; i < min_segments; i++) {
                        int current_length = (i < longer_segments) ? longer_length : shorter_length;
                        
                        seq_count++;
                        fprintf(outseq, ">%d\n", seq_count);
                        fprintf(outseq, "%.*s\n", current_length, line + pos);
                        pos += current_length;
                    }
                }
            }
        }
        free(line);
        fclose(inseq);
        fclose(outseq);
    }
    log_message(INFO, "Reads breaking finished.");
}


// int main(int argc, char *argv[]) {
//     if (argc != 4) {
//         log_message(ERROR, "Usage: %s <input_seq> <output_seq> <break_length>", argv[0]);
//         exit(EXIT_FAILURE);
//     }

//     const char *input_seq = argv[1];
//     const char *output_seq = argv[2];
//     int break_length = atoi(argv[3]);

//     BreakLongReads(input_seq, output_seq, break_length);

//     return 0;
// }

