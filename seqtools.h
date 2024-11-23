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

#ifndef SEQTOOLS_H
#define SEQTOOLS_H


/* get_subsample.c */
void subsample(const char* output, const char* corrected_seq, double factor, int seed);

/* fastq2fa.c */
void fq2fa(const char* filename, const char* outfilename);

/* break_long_reads.c */
void BreakLongReads(const char *input_seq, const char *output_seq, int break_length);


/* correct_sequences.c */
void canu_correct(const char* canu_path, const char* input_seq, 
                    int genomeSize, const char* output_dir, 
                    const char* readstype, int cpu);
    
void nextdenovo_correct(const char* nextdenovo_path, const char* canu_path, 
                        const char* input_seq, const char* cfg, int cfg_flag, 
                        const char* output_dir, char* readstype, char* seqtype, 
                        int cpu, int genomeSize);

/* runassembly.c */
void run_Assembly(const char *sif_path, int cpu, const char *assembly_seq, 
                    const char *output_path, int mi, int ml, int mem, float genomesize_bp);


#endif