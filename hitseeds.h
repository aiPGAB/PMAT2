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


#ifndef HITSEEDS_H
#define HITSEEDS_H


typedef struct {
    char gene[100];
    int length;
} PCGs;


typedef struct {
    char *gene;
    char *ctg;
    float identity;
    int length;
} BlastInfo;

typedef struct {
    char *gene;
    char **ctg;
    int num_ctg;
    float *score;
    int *ctglen;
    float *ctgdep;
} PcgCtgs;


typedef struct {
    char *ctg;
    float score;
    int ctglen;
    float ctgdep;
} SortPcgCtgs;


typedef struct {
    int ctgsmp;
    char *ctg;
    int len;
    float depth;
    float score;
} CtgDepth;


/* find the candidate seeds for the given contigs */
void MtHitseeds(const char* exe_path, const char* organelles_type, const char* all_contigs, 
             const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, int** candidate_seeds, int* ctg_threshold, float filter_depth);

void PtHitseeds(const char* exe_path, const char* organelles_type, const char* all_contigs, 
             const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, int* candidate_seeds, int ctg_threshold, float filter_depth); 
             
void AnHitseeds(const char* exe_path, const char* organelles_type, const char* all_contigs, 
             const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, int** candidate_seeds, int* ctg_threshold, float filter_depth); 

#endif