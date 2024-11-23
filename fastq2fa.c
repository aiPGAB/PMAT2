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
#include <zlib.h>
#include <time.h>

#include "kseq.h"
#include "log.h"

KSEQ_INIT(gzFile, gzread)

void fq2fa(const char* filename, const char* outfilename) {
    gzFile fp;
    kseq_t* seq;
    FILE* out;
    clock_t start, end;
    double cpu_time_used;
    int seq_count = 0;

    fp = gzopen(filename, "r");
    if (!fp) {
        log_message(ERROR, "Error opening input file");
        exit(EXIT_FAILURE);
    }

    seq = kseq_init(fp);

    out = fopen(outfilename, "w");
    if (!out) {
        log_message(ERROR, "Error opening output file");
        kseq_destroy(seq);
        gzclose(fp);
        exit(EXIT_FAILURE);
    }

    log_message(INFO, "Format conversion in progress...");

    start = clock();

    while (kseq_read(seq) >= 0) {
        seq_count++;
        fprintf(out, ">%d\n%s\n", seq_count, seq->seq.s);
    }

    fclose(out);
    kseq_destroy(seq);
    gzclose(fp);

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    log_message(INFO, "Conversion complete.");
    char time_msg[100];
    snprintf(time_msg, sizeof(time_msg), "Conversion time: %.2f s", cpu_time_used);
    log_message(INFO, time_msg);
}


#ifndef FQ2FA_MAIN
int main(int argc, char* argv[]) {
    
    if (argc != 2) {
        log_message(ERROR, "Usage: %s <input_file>", argv[0]);
        return 1;
    }

    const char* outputfile = replace_extension(argv[1], ".format.fa");
    fq2fa(argv[1], outputfile);
    
    return 0;
}
#endif
