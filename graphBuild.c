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
#include <libgen.h> /* dirname */
#include <math.h>

#include "log.h"
#include "version.h"
#include "misc.h"
#include "seqtools.h"
#include "hitseeds.h"
#include "BFSseed.h"
#include "graphtools.h"
#include "gkmer.h"
#include "pmat.h"



void graphBuild(const char* exe_path, graphBuildArgs* opts) {

    log_message(INFO, "graphBuild in progress...");


    /* mkdir output directory */
    char* gfa_dir;
    mkdirfiles(opts->output_file);
    gfa_dir = (char*) malloc(sizeof(*gfa_dir) * (snprintf(NULL, 0, "%s/gfa_result", opts->output_file) + 1));
    sprintf(gfa_dir, "%s/gfa_result", opts->output_file);
    mkdirfiles(gfa_dir);
    free(gfa_dir);

    /* find seeds */
    FILE *graph_file = fopen(opts->assembly_graph, "r");
    if (graph_file == NULL) {
        log_message(ERROR, "Failed to open file: %s", opts->assembly_graph);
        exit(EXIT_FAILURE);
    }

    uint64_t i;
    int num_links = 0;
    size_t len = 0;
    char *line = NULL;
    int flag_C = 0;
    int num_ctg = 0;
    while (getline(&line, &len, graph_file) != -1) {
        if (line[0] == 'C') {
            num_links++;
            flag_C = 1;
        } else if (flag_C) {
            break;
        } else {
            num_ctg++;
        }
    }
    rewind(graph_file);


    Ctglinks* ctglinks = calloc(num_links, sizeof(Ctglinks));
    CtgDepth* ctgdepth = calloc(num_ctg, sizeof(CtgDepth));
    int ctg_arr[200];
    int ctg_arr_idx = 0;
    int ctglink_idx = 0;
    // int ctg_idx = 0;
    int log_idx = 0;
    int log_len = 0;
    flag_C = 0;
    while (getline(&line, &len, graph_file) != -1) {
        if (line[0] == 'C') {
            char* token = strtok(line, "\t");
            token = strtok(NULL, "\t");
            ctglinks[ctglink_idx].lctg = atoi(token);
            token = strtok(NULL, "\t");
            ctglinks[ctglink_idx].lutr = strdup(token);
            token = strtok(NULL, "\t");
            ctglinks[ctglink_idx].rctg = atoi(token);
            token = strtok(NULL, "\t");
            ctglinks[ctglink_idx].rutr = strdup(token);
            token = strtok(NULL, "\t");
            ctglinks[ctglink_idx].linkdepth = atof(token);
            ctglink_idx++;
            flag_C = 1;
        } else if (flag_C) {
            break;
        } else {
            char *token = strtok(line, "\t");
            int ctg_id = atoi(token);
            ctgdepth[ctg_id - 1].ctgsmp = ctg_id;
            token = strtok(NULL, "\t");
            ctgdepth[ctg_id - 1].ctg = strdup(token);
            token = strtok(NULL, "\t");
            ctgdepth[ctg_id - 1].len = atoi(token);
            token = strtok(NULL, "\t");
            ctgdepth[ctg_id - 1].depth = atof(token);
            ctgdepth[ctg_id - 1].score = sqrt(sqrt(ctgdepth[ctg_id - 1].depth) * ctgdepth[ctg_id - 1].len);
            if (ctgdepth[ctg_id - 1].len > log_len) {
                log_len = ctgdepth[ctg_id - 1].len;
                log_idx = ctg_id - 1;
            }
            if (ctgdepth[ctg_id - 1].len > 5000 && ctg_arr_idx < 200) {
                ctg_arr[ctg_arr_idx] = ctgdepth[ctg_id - 1].depth;
                ctg_arr_idx++;
            }
        }
    }
    fclose(graph_file);
    free(line);
    
    // uint64_t longassembly_bp = getFileSize(cut_seq);
    // float seq_depth = (float)longassembly_bp / (float)genomesize_bp;

    float seq_depth = findMedian(ctg_arr, 200);
    float filter_depth;
    if (opts->depth != -1) {
        filter_depth = opts->depth;
    }
    log_message(INFO, "Number of contigs: %d", num_ctg);
    log_message(INFO, "Longest contig: %s %dbp", ctgdepth[log_idx].ctg, ctgdepth[log_idx].len);
    log_message(INFO, "Sequence depth: %.2f", seq_depth);
    // log_message(INFO, "Contig filter depth: %.2f", filter_depth);

    int num_dynseeds = 0;
    int* dynseeds = NULL;
    if (opts->seedCount == 0) {

        if (opts->taxo == 0) {

            if (strcmp(opts->organelles, "mt") == 0) {

                /* pt */
                int pt_num_dynseeds = 0;
                int* pt_dynseeds = NULL;
                int pt_ctg_threshold;
                pt_ctg_threshold = 1;
                pt_dynseeds = calloc(pt_ctg_threshold, sizeof(int));
                PtHitseeds(exe_path, "pt", opts->assembly_fna, opts->output_file, opts->cpu, num_ctg, ctgdepth, pt_dynseeds, pt_ctg_threshold, 2*seq_depth, 1);
                for (i = 0; i < pt_ctg_threshold; i++) {
                    if (pt_dynseeds[i] != 0) {
                        pt_num_dynseeds++;
                    }
                }
                
                int pt_mainseeds_num = 0;
                int* pt_mainseeds = NULL;
                if (pt_num_dynseeds > 0) {
                    filter_depth = 0.3 * ctgdepth[pt_dynseeds[0] - 1].depth;
                    BFSlinks* pt_bfslinks = NULL;
                    int pt_num_BFSlinks = 0;
                    BFSseeds("pt", num_links, num_ctg, ctglinks, ctgdepth, &pt_num_dynseeds, &pt_dynseeds, seq_depth, filter_depth, &pt_bfslinks, &pt_num_BFSlinks);
                    pt_mainseeds = (int*) malloc(sizeof(int) * pt_num_BFSlinks * 2);
                    optgfa(exe_path, pt_num_dynseeds, &pt_dynseeds, &pt_bfslinks, &pt_num_BFSlinks, ctgdepth, opts->output_file, 
                            opts->assembly_fna, opts->assembly_graph, "pt", &pt_mainseeds_num, &pt_mainseeds, 0, NULL, opts->taxo, filter_depth, opts->cutseq);
                    
                    free(pt_bfslinks); 
                }
                free(pt_dynseeds);

                int ctg_threshold = 6;
                dynseeds = calloc(ctg_threshold, sizeof(int));
                if (opts->depth == -1) {
                    if ((2 * seq_depth) > 1.5) {
                        filter_depth = 2 * seq_depth;
                    } else {
                        filter_depth = 1.5;
                    }
                } else {
                    filter_depth = opts->depth;
                }
                HitSeeds(exe_path, "mt", opts->assembly_fna, opts->output_file, opts->cpu, num_ctg, ctgdepth, &dynseeds, &ctg_threshold, filter_depth, opts->taxo, 1);

                for (i = 0; i < ctg_threshold; i++) {
                    if (dynseeds[i] != 0) {
                        num_dynseeds++;
                    }
                }

                if (num_dynseeds > 0) {
                    BFSlinks* bfslinks = NULL;
                    int num_BFSlinks = 0;
                    BFSseeds("mt", num_links, num_ctg, ctglinks, ctgdepth, &num_dynseeds, &dynseeds, seq_depth, filter_depth, &bfslinks, &num_BFSlinks);
                    int mt_mainseeds_num = 0;
                    int* mainseeds = (int*) malloc(sizeof(int) * num_BFSlinks);
                    optgfa(exe_path, num_dynseeds, &dynseeds, &bfslinks, &num_BFSlinks, ctgdepth, opts->output_file, opts->assembly_fna, opts->assembly_graph, "mt", &mt_mainseeds_num, &mainseeds, pt_mainseeds_num, pt_mainseeds, opts->taxo, filter_depth, opts->cutseq);
                    
                    free(bfslinks);
                }
                free(dynseeds);
                free(pt_mainseeds);
                
            } else if (strcmp(opts->organelles, "pt") == 0) {
                if (opts->depth == -1) {
                    filter_depth = 2*seq_depth;
                } else {
                    filter_depth = opts->depth;
                }

                int ctg_threshold = 1;
                dynseeds = calloc(ctg_threshold, sizeof(int));
                PtHitseeds(exe_path, "pt", opts->assembly_fna, opts->output_file, opts->cpu, num_ctg, ctgdepth, dynseeds, ctg_threshold, filter_depth, 1);
                for (i = 0; i < ctg_threshold; i++) {
                    if (dynseeds[i] != 0) {
                        num_dynseeds++;
                    }
                }
                if (num_dynseeds > 0) {
                    if (opts->depth == -1) {
                        filter_depth = 0.3 * ctgdepth[dynseeds[0] - 1].depth;
                    } else {
                        filter_depth = opts->depth;
                    }

                    BFSlinks* bfslinks = NULL;
                    int num_BFSlinks = 0;
                    BFSseeds("pt", num_links, num_ctg, ctglinks, ctgdepth, &num_dynseeds, &dynseeds, seq_depth, filter_depth, &bfslinks, &num_BFSlinks);
                    int pt_mainseeds_num = 0;
                    int* mainseeds = (int*) malloc(sizeof(int) * num_BFSlinks);
                    optgfa(exe_path, num_dynseeds, &dynseeds, &bfslinks, &num_BFSlinks, ctgdepth, opts->output_file, opts->assembly_fna, opts->assembly_graph, "pt", &pt_mainseeds_num, &mainseeds, 0 ,NULL, opts->taxo, filter_depth, opts->cutseq);
                    
                    free(bfslinks);
                    free(mainseeds);
                }
                free(dynseeds);
            }
        } else if (opts->taxo == 1) {
                /* mt */
                if (opts->depth == -1) {
                    filter_depth = 2*seq_depth;
                } else {
                    filter_depth = opts->depth;
                }
                int ctg_threshold = 1;
                dynseeds = calloc(ctg_threshold, sizeof(int));
                HitSeeds(exe_path, "mt", opts->assembly_fna, opts->output_file, opts->cpu, num_ctg, ctgdepth, &dynseeds, &ctg_threshold, filter_depth, opts->taxo, 1);
                for (i = 0; i < ctg_threshold; i++) {
                    if (dynseeds[i] != 0) {
                        num_dynseeds++;
                    }
                }
                if (num_dynseeds > 0) {

                    if (opts->depth == -1) {
                        filter_depth = 0.3 * ctgdepth[dynseeds[0] - 1].depth;
                    } else {
                        filter_depth = opts->depth;
                    }

                    BFSlinks* bfslinks = NULL;
                    int num_BFSlinks = 0;
                    BFSseeds("mt", num_links, num_ctg, ctglinks, ctgdepth, &num_dynseeds, &dynseeds, seq_depth, filter_depth, &bfslinks, &num_BFSlinks);
                    int mt_mainseeds_num = 0;
                    int* mainseeds = (int*) malloc(sizeof(int) * num_BFSlinks);
                    optgfa(exe_path, num_dynseeds, &dynseeds, &bfslinks, &num_BFSlinks, ctgdepth, opts->output_file, opts->assembly_fna, opts->assembly_graph, "mt", &mt_mainseeds_num, &mainseeds, 0, NULL, opts->taxo, filter_depth, opts->cutseq);
                    
                    free(bfslinks);
                    free(mainseeds);
                }
                free(dynseeds);

        } else if (opts->taxo == 2) {
                /* mt */
                if (opts->depth == -1) {
                    filter_depth = 2*seq_depth;
                } else {
                    filter_depth = opts->depth;
                }                
                int ctg_threshold = 1;
                dynseeds = calloc(ctg_threshold, sizeof(int));
                HitSeeds(exe_path, "mt", opts->assembly_fna, opts->output_file, opts->cpu, num_ctg, ctgdepth, &dynseeds, &ctg_threshold, filter_depth, opts->taxo, 1);
                for (i = 0; i < ctg_threshold; i++) {
                    if (dynseeds[i] != 0) {
                        num_dynseeds++;
                    }
                }
                if (num_dynseeds > 0) {

                    if (opts->depth == -1) {
                        filter_depth = 0.3 * ctgdepth[dynseeds[0] - 1].depth;
                    } else {
                        filter_depth = opts->depth;
                    }

                    BFSlinks* bfslinks = NULL;
                    int num_BFSlinks = 0;
                    BFSseeds("mt", num_links, num_ctg, ctglinks, ctgdepth, &num_dynseeds, &dynseeds, seq_depth, filter_depth, &bfslinks, &num_BFSlinks);
                    int mt_mainseeds_num = 0;
                    int* mainseeds = (int*) malloc(sizeof(int) * num_BFSlinks);
                    optgfa(exe_path, num_dynseeds, &dynseeds, &bfslinks, &num_BFSlinks, ctgdepth, opts->output_file, opts->assembly_fna, opts->assembly_graph, "mt", &mt_mainseeds_num, &mainseeds, 0, NULL, opts->taxo, filter_depth, opts->cutseq);
                    
                    free(bfslinks);
                    free(mainseeds);
                }
                free(dynseeds);
            }

    } else {
        num_dynseeds = opts->seedCount;
        dynseeds = calloc(num_dynseeds, sizeof(int));
        for (i = 0; i < num_dynseeds; i++) 
        {
            if (opts->seeds[i] < 1 || opts->seeds[i] > num_ctg) {
                log_message(ERROR, "Seed %d is out of range, please check the input seeds.", opts->seeds[i]);
                exit(EXIT_FAILURE);
            }
            dynseeds[i] = opts->seeds[i];
        }
        if (opts->taxo == 0) {
            if (strcmp(opts->organelles, "mt") == 0) {
                /* pt */
                int pt_num_dynseeds = 0;
                int* pt_dynseeds = NULL;
                int pt_ctg_threshold;
                pt_ctg_threshold = 1;
                pt_dynseeds = calloc(pt_ctg_threshold, sizeof(int));
                PtHitseeds(exe_path, "pt", opts->assembly_fna, opts->output_file, opts->cpu, num_ctg, ctgdepth, pt_dynseeds, pt_ctg_threshold, 2*seq_depth, 1);
                for (i = 0; i < pt_ctg_threshold; i++) {
                    if (pt_dynseeds[i] != 0) {
                        pt_num_dynseeds++;
                    }
                }
                
                int pt_mainseeds_num = 0;
                int* pt_mainseeds = NULL;
                if (pt_num_dynseeds > 0) {
                    filter_depth = (ctgdepth[pt_dynseeds[0] - 1].depth > 4) ? 0.3 * ctgdepth[pt_dynseeds[0] - 1].depth : 2;
                    BFSlinks* pt_bfslinks = NULL;
                    int pt_num_BFSlinks = 0;
                    BFSseeds("pt", num_links, num_ctg, ctglinks, ctgdepth, &pt_num_dynseeds, &pt_dynseeds, seq_depth, filter_depth, &pt_bfslinks, &pt_num_BFSlinks);
                    pt_mainseeds = (int*) malloc(sizeof(int) * pt_num_BFSlinks * 2);
                    optgfa(exe_path, pt_num_dynseeds, &pt_dynseeds, &pt_bfslinks, &pt_num_BFSlinks, ctgdepth, opts->output_file, 
                            opts->assembly_fna, opts->assembly_graph, "pt", &pt_mainseeds_num, &pt_mainseeds, 0, NULL, opts->taxo, filter_depth, opts->cutseq);
                    
                    free(pt_bfslinks); 
                }
                free(pt_dynseeds);
                
                /* mt */
                if (opts->depth == -1) {
                    if ((2 * seq_depth) > 1.5) {
                        filter_depth = 2 * seq_depth;
                    } else {
                        filter_depth = 1.5;
                    }
                } else {
                    filter_depth = opts->depth;
                }

                BFSlinks* bfslinks = NULL;
                int num_BFSlinks = 0;
                BFSseeds("mt", num_links, num_ctg, ctglinks, ctgdepth, &num_dynseeds, &dynseeds, seq_depth, filter_depth, &bfslinks, &num_BFSlinks);
                int mt_mainseeds_num = 0;
                int* mainseeds = (int*) malloc(sizeof(int) * num_BFSlinks);
                optgfa(exe_path, num_dynseeds, &dynseeds, &bfslinks, &num_BFSlinks, ctgdepth, opts->output_file, opts->assembly_fna, opts->assembly_graph, "mt", &mt_mainseeds_num, &mainseeds, pt_mainseeds_num, pt_mainseeds, opts->taxo, filter_depth, opts->cutseq);
                
                free(bfslinks);
                free(mainseeds);
                free(pt_mainseeds);
            } else if (strcmp(opts->organelles, "pt") == 0) {
                if (opts->depth == -1) {
                    filter_depth = (ctgdepth[dynseeds[0] - 1].depth > 4) ? 0.3 * ctgdepth[dynseeds[0] - 1].depth : 2;
                    if (filter_depth < 3*seq_depth) filter_depth = 3*seq_depth;
                } else {
                    filter_depth = opts->depth;
                }

                BFSlinks* bfslinks = NULL;
                int num_BFSlinks = 0;
                BFSseeds("pt", num_links, num_ctg, ctglinks, ctgdepth, &num_dynseeds, &dynseeds, seq_depth, filter_depth, &bfslinks, &num_BFSlinks);
                int pt_mainseeds_num = 0;
                int* mainseeds = (int*) malloc(sizeof(int) * num_BFSlinks);
                optgfa(exe_path, num_dynseeds, &dynseeds, &bfslinks, &num_BFSlinks, ctgdepth, opts->output_file, opts->assembly_fna, opts->assembly_graph, "pt", &pt_mainseeds_num, &mainseeds, 0 ,NULL, opts->taxo, filter_depth, opts->cutseq);
                
                free(bfslinks);
                free(mainseeds);
            }
        } else if (opts->taxo == 1 || opts->taxo == 2) {
            if (strcmp(opts->organelles, "mt") == 0) {
                if (opts->depth == -1) {
                    filter_depth = (ctgdepth[dynseeds[0] - 1].depth > 4) ? 0.3 * ctgdepth[dynseeds[0] - 1].depth : 2;
                    if (filter_depth < 2*seq_depth) filter_depth = 2*seq_depth;
                } else {
                    filter_depth = opts->depth;
                }

                BFSlinks* bfslinks = NULL;
                int num_BFSlinks = 0;
                BFSseeds("mt", num_links, num_ctg, ctglinks, ctgdepth, &num_dynseeds, &dynseeds, seq_depth, filter_depth, &bfslinks, &num_BFSlinks);
                int mt_mainseeds_num = 0;
                int* mainseeds = (int*) malloc(sizeof(int) * num_BFSlinks);
                optgfa(exe_path, num_dynseeds, &dynseeds, &bfslinks, &num_BFSlinks, ctgdepth, opts->output_file, opts->assembly_fna, opts->assembly_graph, "mt", &mt_mainseeds_num, &mainseeds, 0, NULL, opts->taxo, filter_depth, opts->cutseq);
                
                free(bfslinks);
                free(mainseeds);
            }
        }
    }

    for (i = 0; i < num_ctg; i++) 
    {
        free(ctgdepth[i].ctg);
    }
    free(ctgdepth);

    for (i = 0; i < num_links; i++) 
    {
        free(ctglinks[i].lutr);
        free(ctglinks[i].rutr);
    }
    
    free(ctglinks);
}