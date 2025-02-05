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

#ifndef ORGASS_H
#define ORGASS_H

#include "khash.h"
#include "hitseeds.h"

KHASH_MAP_INIT_INT(hash_count, int)

typedef struct {
    int contig_id;
    int num_genes;
    char **gene_list;
} ContigGeneStats;

typedef struct {
    int total_contigs;
    int total_genes;
    int unique_genes;
    int duplicate_contigs;
    int* duplicate_ids;
    int num_duplicates;
    ContigGeneStats *contig_stats;
} AssessResult;

typedef struct {
    int *contig_array;
    int arr_size;
    int repeat_contig;
} MtStructure;


void orgAss(const char* exe_path, const char* all_contigs, CtgDepth *ctg_depth, 
            const char* output_path, const int* contig_ids, 
            int num_contigs, const char* organelle_type, int taxo);

#endif // ORGASS_H