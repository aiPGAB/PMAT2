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

#ifndef GRAPHTOOLS_H
#define GRAPHTOOLS_H

#include <stdint.h>
#include "hitseeds.h"
#include "BFSseed.h"
#include "khash.h"


typedef struct {
    int pt_nodenum;
    int mt_nodenum;
    int uniq_pt_nodenum;
    int uniq_mt_nodenum;
    int* path_node;
    int* path_utr;
    uint32_t node_num;
    uint64_t path_len;
    uint64_t uniq_mt_pathlen;
    uint32_t inval_num;
    int type;
} pathScore;

typedef struct {
    BFSlinks* links;
    int num_links;
    int* node;
    int num_nodes;
} BFSstructure;

KHASH_MAP_INIT_INT(Ha_structures, BFSstructure*)
KHASH_MAP_INIT_INT(Ha_nodeseq, char*)

typedef struct {
    int ctg;
    char *seq;
} fnainfo;

/* addseq: add sequence to the fna */
void addseq(const char* allgraph, const char* all_fna, CtgDepth* ctgdepth);

/* raw gfa && main gfa */
void optgfa(const char* exe_path, int num_dynseeds, int** dynseeds, BFSlinks** bfslinks, int* num_bfslinks, 
            CtgDepth* ctgdepth, const char* output, const char* all_fna, const char* allgraph, 
            const char* organelles_type, int* mainseeds_num, int** mainseeds, int interfering_ctg_num, 
            int* interfering_ctg, int taxo, float filter_depth, char* cutseq);

/* findSpath: find the shortest path between two contigs */
void findSpath(int node1, int node1utr, int node2, int node2utr, 
              int main_num, BFSlinks* mainlinks, CtgDepth *ctg_depth);

/* BFS_structure: find the BFS structure of the graph */
uint32_t bfs_structure(int node_num, int link_num, BFSlinks* links, int* node_arry, khash_t(Ha_structures)* h_structures);

/* findMpath: find the most likely path between two contigs */
void findMpath(int node1, int node1utr, int node2, int node2utr, int main_num, BFSlinks* mainlinks, CtgDepth *ctg_depth, 
    int* mt_contigs, int mt_num, int* pt_contigs, int pt_num, int* flag_err, float* mt_ratio, int taxo, pathScore *struc_path);

/* copy BFSlinks */
void copy_BFSlinks(BFSlinks* dest, const BFSlinks* src);
// void free_BFSlinks(BFSlinks *links, int num_links);

/* convert path to fasta */
void path2fa(pathScore *path, int ps_num, khash_t(Ha_nodeseq)* node_seq, const char *output);

#endif