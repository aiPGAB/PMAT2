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

#ifndef PMAT_H
#define PMAT_H



/* PMAT function: autoMito */

typedef struct {
    char *input_file;
    char *output_file;
    char *seqtype;
    char *runassembly;
    char *genomesize;
    char *organelles;    // Organelles type (mt/pt)
    int8_t task;
    char *correct_software;
    char *canu_path;
    char *nextdenovo_path;
    char *cfg_file;
    int8_t cfg_flag;
    double factor;
    int32_t seed;
    int32_t breaknum;
    int8_t mi;
    int8_t ml;
    int32_t cpu;
    int8_t taxo;
    int8_t mem;
    int8_t kmersize;
} autoMitoArgs;


typedef struct {
    char *graphinfo;     // Graph info file
    char *assembly_graph;  // Contig graph file
    char *assembly_fna;   // All contigs file
    char *subsample;      // Subsample file
    char *cutseq;
    char *output_file;   // Output file
    char *organelles;    // Organelles type (mt/pt)
    int8_t taxo;
    float depth;         // Depth of sequencing
    int *seeds;          // Array of seed values
    int seedCount;       // Number of seeds
    int cpu;             // Number of CPUs
} graphBuildArgs;


void autoMito(const char* exe_path, autoMitoArgs* opts);


/* PMAT function: graphBuild */
void graphBuild(const char* exe_path, graphBuildArgs* opts);


#endif /* PMAT_H */