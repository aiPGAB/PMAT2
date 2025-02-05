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

#include "khash.h"


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


typedef struct {
    int num_genes;
    char **gene_ids;
    float *identity;
    int *align_len;
} ContigGenes;


typedef struct {
    char* query_id;
    char* gene_id;
    float identity;
    int align_len;
    int mismatch;
    int gap;
    int qstart;
    int qend;
    int sstart;
    int send;
    float evalue;
    float score;
    int direction;
} BlastDircMatch;

// Conserved genes
static const char *plt_mtpcg[] = {
    "nad1",  "nad2",  "nad3",  "nad4",  "nad5",  "nad6",  "nad7",  "nad9",
    "nad4L", "cob",   "cox1",  "cox2",  "cox3",  "atp1",  "atp4",  "atp6",
    "atp8",  "atp9",  "ccmB", "ccmC", "ccmFc", "ccmFn", "mttB", "matR"
};
static const int plt_mtpcg_len[] = {950, 1467, 357, 1488, 2010, 618, 1185, 573, 
                                    303, 1182, 1584, 783, 798, 1524, 594, 768, 
                                    470, 225, 621, 720, 1320, 1785, 816, 2000};
static const int plt_mtpcg_num = 24;

// Conserved genes
static const char *fug_mtpcg[] = {
    "atp6", "atp8", "atp9", "cob", "cox1", "cox2", "cox3", "nad1", 
    "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6"
};
static const int fug_mtpcg_len[] = {770, 150, 225, 1150, 1520, 750, 810, 1030, 1600, 430, 880, 290, 1950, 630};
static const int fug_mtpcg_num = 14;

// Conserved genes
static const char *anl_mtpcg[] = {
    "atp6", "atp8", "CytB", "cox1", "cox2", "cox3", "nad1", 
    "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6"
};
static const int anl_mtpcg_len[] = {650, 165, 1100, 1500, 1520, 680, 900, 1000, 350, 1350, 290, 1800, 510};
static const int anl_mtpcg_num = 13;

KHASH_MAP_INIT_INT(ctg_genes, ContigGenes*)

/* find the candidate seeds for the given contigs */
void mrun_blastn(const char *all_contigs, const char *db_path, char *final_out, int num_threads, int *num_hits);

void PtHitseeds(const char* exe_path, const char* organelles_type, const char* all_contigs, 
             const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, 
             int* candidate_seeds, int ctg_threshold, float filter_depth, int verbose); 

void HitSeeds(const char* exe_path, const char* organelles_type, const char* all_contigs, 
            const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, 
            int** candidate_seeds, int* ctg_threshold, float filter_depth, int taxo, int verbose);

#endif