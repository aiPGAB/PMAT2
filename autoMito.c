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
#include <stdint.h>
#include <libgen.h> /* dirname */
#include <math.h>

#include "log.h"
#include "kalloc.h"
#include "version.h"
#include "misc.h"
#include "seqtools.h"
#include "hitseeds.h"
#include "BFSseed.h"
#include "graphtools.h"
#include "gkmer.h"
#include "pmat.h"



void autoMito(const char* exe_path, autoMitoArgs* opts) {

    log_message(INFO, "autoMito in progress...");

    // get sequence type
    char* readstype;
    if (strcmp(opts->seqtype, "hifi") == 0) {
        readstype = "hifi";
    } else if (strcmp(opts->seqtype, "ont") == 0) {
        readstype = "nanopore";
    } else if (strcmp(opts->seqtype, "clr") == 0) {
        readstype = "pacbio";
    } else {
        log_message(ERROR, "Invalid sequence type: %s", opts->seqtype);
        exit(EXIT_FAILURE);
    }

    mkdirfiles(opts->output_file);
    uint64_t genomesize_bp = 0;
    uint64_t i;
    if (strcmp(readstype, "hifi") != 0) {
        if (opts->genomesize != NULL) {
            to_lower(opts->genomesize);

            int glen = strlen(opts->genomesize);
            if (opts->genomesize[strlen(opts->genomesize) - 1] == 'g') {
                opts->genomesize[glen - 1] = '\0';
                genomesize_bp = (uint64_t)(atof(opts->genomesize) * 1000000000);

            } else if (opts->genomesize[strlen(opts->genomesize) - 1] == 'm') {
                opts->genomesize[glen - 1] = '\0';
                genomesize_bp = (uint64_t)(atof(opts->genomesize) * 1000000);

            } else if (opts->genomesize[strlen(opts->genomesize) - 1] == 'k') {
                opts->genomesize[glen - 1] = '\0';
                genomesize_bp = (uint64_t)(atof(opts->genomesize) * 1000);

            } else if (is_digits(opts->genomesize)) {
                genomesize_bp = (uint64_t)atof(opts->genomesize);
            } else {
                log_message(ERROR, "Invalid genome size: %s", opts->genomesize);
                exit(EXIT_FAILURE);
            }

        } else {

            log_message(INFO, "Kmer frequency counting...");
            char* gkmer_dir;
            char* gkmer_histo;

            char* dir = dirname(strdup(exe_path));
            char* genomescope = (char*)malloc(sizeof(*genomescope) * (snprintf(NULL, 0, "%s/lib/genomescope.R", dir) + 1));
            sprintf(genomescope, "%s/lib/genomescope.R", dir);
            if (which_executable(genomescope) == 0) {
                log_message(ERROR, "Failed to find genomescope.R in %s/lib", dir);
                exit(EXIT_FAILURE);
            }
            free(dir); 

            gkmer_dir = (char*)malloc(sizeof(*gkmer_dir) * (snprintf(NULL, 0, "%s/gkmer", opts->output_file) + 1));
            sprintf(gkmer_dir, "%s/gkmer", opts->output_file);
            mkdirfiles(gkmer_dir);

            gkmer_histo = (char*)malloc(sizeof(*gkmer_histo) * (snprintf(NULL, 0, "%s/gkmer/gkmer_histo.txt", opts->output_file) + 1));
            sprintf(gkmer_histo, "%s/gkmer/gkmer_histo.txt", opts->output_file);

            yak_copt_t opt;
            yak_copt_init(&opt);
            opt.k = opts->kmersize;
            opt.n_thread = opts->cpu;
            
            gkmerAPI(&opt, opts->input_file, gkmer_histo); //*

            char* command = NULL;
            size_t cmd_len = snprintf(NULL, 0, "%s %s %d 15000 %s 1000 0", 
                genomescope, gkmer_histo, opts->kmersize, gkmer_dir) + 1;
            command = (char*) malloc(cmd_len);
            snprintf(command, cmd_len, "%s %s %d 15000 %s 1000 0", 
                genomescope, gkmer_histo, opts->kmersize, gkmer_dir);
            execute_command(command, 0, 1);
            free(genomescope); free(command);
            char* summ = NULL;
            size_t summ_len = snprintf(NULL, 0, "%s/summary.txt", gkmer_dir) + 1;
            summ = (char*) malloc(summ_len);
            snprintf(summ, summ_len, "%s/summary.txt", gkmer_dir);
            checkfile(summ);
            free(gkmer_dir); free(gkmer_histo);

            FILE* fp2 = fopen(summ, "r");
            if (fp2 == NULL) {
                fprintf(stderr, "Failed to open file %s\n", summ);
                exit(EXIT_FAILURE);
            }
            size_t line_len = 0;
            char* line = NULL;
            while (getline(&line, &line_len, fp2)!= -1) {
                if (strstr(line, "Genome Haploid Length") != NULL) {
                    char *max_part = strstr(line, "bp");
                    if (max_part != NULL) {
                        max_part += 2;
                        while (*max_part == ' ' || *max_part == '\t') max_part++;
                        opts->genomesize = (char*) malloc(strlen(max_part) + 1);
                        sscanf(max_part, "%s", opts->genomesize);
                        remove_commas(opts->genomesize);
                        genomesize_bp = (uint64_t)atof(opts->genomesize);
                        if (genomesize_bp < 1) {
                            log_message(ERROR, "Invalid genome size: %s bp (failed to converge)", opts->genomesize);
                            free(opts->genomesize);
                            exit(EXIT_FAILURE);
                        }
                    }
                } else if (strstr(line, "Model Fit") != NULL) {
                    double fit_rate = 0;
                    sscanf(line, "%*s %*s %lf", &fit_rate); 
                    if (fit_rate < 90) {
                        log_message(ERROR, "Invalid genome size: %s bp (failed to converge)", opts->genomesize);
                        free(opts->genomesize);
                        exit(EXIT_FAILURE);
                    }
                }
            }
            log_message(INFO, "Kmer size: %d; Estimated genome size: %llu bp", opts->kmersize, genomesize_bp);

            if (opts->genomesize != NULL) free(opts->genomesize);
            fclose(fp2); 
            free(line); free(summ); 
        }
    }

    /* mkdir output directory */
    char* subsample_dir;
    char* assembly_dir;
    char* gfa_dir;

    size_t len_subsample = snprintf(NULL, 0, "%s/subsample", opts->output_file) + 1;
    subsample_dir = (char*)malloc(len_subsample);
    sprintf(subsample_dir, "%s/subsample", opts->output_file);
    mkdirfiles(subsample_dir);

    size_t len_assembly = snprintf(NULL, 0, "%s/assembly_result", opts->output_file) + 1;
    assembly_dir = (char*)malloc(len_assembly);
    sprintf(assembly_dir, "%s/assembly_result", opts->output_file);

    size_t len_gfa = snprintf(NULL, 0, "%s/gfa_result", opts->output_file) + 1;
    gfa_dir = (char*)malloc(len_gfa);
    sprintf(gfa_dir, "%s/gfa_result", opts->output_file);
    mkdirfiles(gfa_dir);

    char* high_quality_seq = (char*)malloc(sizeof(*high_quality_seq) * (snprintf(NULL, 0, "%s/subsample/PMAT_assembly_seq.fa", opts->output_file) + 1));
    sprintf(high_quality_seq, "%s/subsample/PMAT_assembly_seq.fa", opts->output_file);

    if (strcmp(opts->seqtype, "hifi") == 0) {
        fq2fa(opts->input_file, high_quality_seq); //*
    } else if (strcmp(opts->seqtype, "ont") == 0 || strcmp(opts->seqtype, "clr") == 0) {
        if (opts->task == 1) {
            char* correct_seq = (char*)malloc(sizeof(*correct_seq) * (snprintf(NULL, 0, "%s/correct_out/PMAT.correctedReads.fasta%s", 
                opts->output_file, (strcmp(opts->correct_software, "canu") == 0 ? ".gz" : "")) + 1));
            sprintf(correct_seq, "%s/correct_out/PMAT.correctedReads.fasta%s", 
                opts->output_file, (strcmp(opts->correct_software, "canu") == 0 ? ".gz" : ""));

            to_lower(opts->correct_software);
            if (strcmp(opts->correct_software, "canu") == 0) {
                canu_correct(opts->canu_path, opts->input_file, 
                            genomesize_bp, opts->output_file, 
                            readstype, opts->cpu);
            } else if (strcmp(opts->correct_software, "nextdenovo") == 0) {
                nextdenovo_correct(opts->nextdenovo_path, opts->canu_path, 
                                    opts->input_file, opts->cfg_file, opts->cfg_flag, 
                                    opts->output_file, opts->seqtype, readstype, opts->cpu, genomesize_bp);
            }
            checkfile(correct_seq);
            fq2fa(correct_seq, high_quality_seq);
            free(correct_seq);

        } else if (opts->task == 0) {
            fq2fa(opts->input_file, high_quality_seq);
        } else {
            log_message(ERROR, "Invalid task type: %d", opts->task);
            exit(EXIT_FAILURE);
        }
    } else {
        log_message(ERROR, "Invalid sequence type (hifi/ont/clr): %s ", opts->seqtype);
        exit(EXIT_FAILURE);
    }

    // select subsample
    char* subsample_seq = (char*)malloc(sizeof(*subsample_seq) * (snprintf(NULL, 0, "%s/subsample/PMAT_subsample_seq.fa", opts->output_file) + 1));
    sprintf(subsample_seq, "%s/subsample/PMAT_subsample_seq.fa", opts->output_file);

    char* cut_seq = (char*)malloc(sizeof(*cut_seq) * (snprintf(NULL, 0, "%s/subsample/PMAT_cut_seq.fa", opts->output_file) + 1));
    sprintf(cut_seq, "%s/subsample/PMAT_cut_seq.fa", opts->output_file);

    if (opts->factor == 1) {
        BreakLongReads(high_quality_seq, cut_seq, opts->breaknum); //*
        remove_file(high_quality_seq);
    } else {
        subsample(subsample_seq, high_quality_seq, opts->factor, opts->seed); //*
        BreakLongReads(subsample_seq, cut_seq, opts->breaknum); //*
        remove_file(high_quality_seq);
        remove_file(subsample_seq);
    }


    /* run assembly */
    char* dir_pmat = dirname(strdup(exe_path));
    char sif_path[4096];
    snprintf(sif_path, sizeof(sif_path), "%s/container/runAssembly.sif", dir_pmat);
    free(dir_pmat);
    if (is_file(sif_path) == 0) {
        log_message(ERROR, "Failed to find container: %s", sif_path);
        exit(EXIT_FAILURE);
    }
    run_Assembly(sif_path, opts->cpu, cut_seq, opts->output_file, opts->mi, opts->ml, opts->mem, genomesize_bp); //*

    char* assembly_fna = (char*)malloc(sizeof(*assembly_fna) * (snprintf(NULL, 0, "%s/assembly_result/PMATAllContigs.fna", opts->output_file) + 1));
    sprintf(assembly_fna, "%s/assembly_result/PMATAllContigs.fna", opts->output_file);

    char* assembly_graph = (char*)malloc(sizeof(*assembly_graph) * (snprintf(NULL, 0, "%s/assembly_result/PMATContigGraph.txt", opts->output_file) + 1));
    sprintf(assembly_graph, "%s/assembly_result/PMATContigGraph.txt", opts->output_file);


    /* find seeds */

    FILE *graph_file = fopen(assembly_graph, "r");
    if (graph_file == NULL) {
        log_message(ERROR, "Failed to open file: %s", assembly_graph);
        exit(EXIT_FAILURE);
    }

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
    int num_taxa = 200;
    if (opts->taxo == 2) num_taxa = 100;
    int ctg_arr[num_taxa];
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
            if (ctgdepth[ctg_id - 1].len > 5000 && ctg_arr_idx < num_taxa) {
                ctg_arr[ctg_arr_idx] = ctgdepth[ctg_id - 1].depth;
                ctg_arr_idx++;
            }
        }
    }
    fclose(graph_file);
    free(line);
    
    /* addseq */
    addseq(assembly_graph, assembly_fna, ctgdepth);

    FILE *fin = fopen(assembly_graph, "r");
    if (!fin) {
        log_message(ERROR, "Failed to open file: %s", assembly_graph);
        exit(EXIT_FAILURE);
    }

    char* assembly_info = (char*)malloc(sizeof(*assembly_info) * (snprintf(NULL, 0, "%s/assembly_result/PMATContiginfo.txt", opts->output_file) + 1));
    sprintf(assembly_info, "%s/assembly_result/PMATContiginfo.txt", opts->output_file);

    FILE *fout = fopen(assembly_info, "w");
    if (!fout) {
        log_message(ERROR, "Failed to open file: %s", assembly_info);
        fclose(fin);
        exit(EXIT_FAILURE);
    }

    size_t line_len = 0;

    while (getline(&line, &line_len, fin) != -1) {
        if (line[0] == 'S') {
            break;
        }
        fputs(line, fout);
    }

    free(line);
    fclose(fin);
    fclose(fout);
    remove_file(assembly_graph);
    rename_file(assembly_info, assembly_graph);
    free(assembly_info);

    // uint64_t longassembly_bp = getFileSize(cut_seq);
    // float seq_depth = longassembly_bp / genomesize_bp;
    float seq_depth = findMedian(ctg_arr, num_taxa);
    float filter_depth;

    log_message(INFO, "Number of contigs: %d", num_ctg);
    log_message(INFO, "Longest contig: %s %dbp", ctgdepth[log_idx].ctg, ctgdepth[log_idx].len);
    log_message(INFO, "Sequence depth: %.2f", seq_depth);

    if (opts->taxo == 0) {

        /* plants */
        /* pt */
        int pt_num_dynseeds = 0;
        int* pt_dynseeds = NULL;
        int pt_ctg_threshold;
        pt_ctg_threshold = 1;
        pt_dynseeds = calloc(pt_ctg_threshold, sizeof(int));
        PtHitseeds(exe_path, "pt", assembly_fna, opts->output_file, opts->cpu, num_ctg, ctgdepth, pt_dynseeds, pt_ctg_threshold, 2*seq_depth, 0);
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
                    assembly_fna, assembly_graph, "pt", &pt_mainseeds_num, &pt_mainseeds, 0, NULL, 0, filter_depth, cut_seq);
            
            // for (int i = 0; i < pt_num_BFSlinks; i++) {
            //     free(pt_bfslinks[i].lctg); free(pt_bfslinks[i].lutr); free(pt_bfslinks[i].rctg); free(pt_bfslinks[i].rutr);
            // }
            free(pt_bfslinks); 
        }
        free(pt_dynseeds);

        if (strcmp(opts->organelles, "mt") == 0) {
            /* mt */

            int mt_ctg_threshold = 6;
            int mt_num_dynseeds = 0;
            int* mt_dynseeds = NULL;
            if ((2 * seq_depth) > 1.5) {
                filter_depth = 2 * seq_depth;
            } else {
                filter_depth = 1.5;
            }
            mt_dynseeds = calloc(mt_ctg_threshold, sizeof(int));
            HitSeeds(exe_path, "mt", assembly_fna, opts->output_file, opts->cpu, num_ctg, ctgdepth, &mt_dynseeds, &mt_ctg_threshold, 1.5*seq_depth, 0, 0);
            for (i = 0; i < mt_ctg_threshold; i++) {
                if (mt_dynseeds[i] != 0) {
                    mt_num_dynseeds++;
                }
            }

            if (mt_num_dynseeds > 0) {
                BFSlinks* mt_bfslinks = NULL;
                int mt_num_BFSlinks = 0;
                BFSseeds("mt", num_links, num_ctg, ctglinks, ctgdepth, &mt_num_dynseeds, &mt_dynseeds, seq_depth, filter_depth, &mt_bfslinks, &mt_num_BFSlinks);
                int mt_mainseeds_num = 0;
                int* mt_mainseeds = (int*) malloc(sizeof(int) * mt_num_BFSlinks * 2);
                optgfa(exe_path, mt_num_dynseeds, &mt_dynseeds, &mt_bfslinks, &mt_num_BFSlinks, ctgdepth, opts->output_file, 
                        assembly_fna, assembly_graph, "mt", &mt_mainseeds_num, &mt_mainseeds, pt_mainseeds_num, pt_mainseeds, 0, filter_depth, cut_seq);
                
                // for (int i = 0; i < mt_num_BFSlinks; i++) {
                //     free(mt_bfslinks[i].lctg); free(mt_bfslinks[i].lutr); free(mt_bfslinks[i].rctg); free(mt_bfslinks[i].rutr);
                // }
                free(mt_bfslinks); 
                // free(mt_mainseeds);
            }
            free(mt_dynseeds);
            free(pt_mainseeds);
        }


    } else if (opts->taxo == 1) {
        /* animals */
        /* mt */
        int mt_ctg_threshold = 1;
        int mt_num_dynseeds = 0;
        int* mt_dynseeds = NULL;
        if ((2 * seq_depth) > 1.5) {
            filter_depth = 2 * seq_depth;
        } else {
            filter_depth = 1.5;
        }
        mt_dynseeds = calloc(mt_ctg_threshold, sizeof(int));
        HitSeeds(exe_path, "mt", assembly_fna, opts->output_file, opts->cpu, num_ctg, ctgdepth, &mt_dynseeds, &mt_ctg_threshold, 2*seq_depth, 1, 0);
        for (i = 0; i < mt_ctg_threshold; i++) {
            if (mt_dynseeds[i] != 0) {
                mt_num_dynseeds++;
            }
        }
        if (mt_num_dynseeds > 0) {
            BFSlinks* mt_bfslinks = NULL;
            int mt_num_BFSlinks = 0;
            filter_depth = (ctgdepth[mt_dynseeds[0] - 1].depth > 4) ? 0.3 * ctgdepth[mt_dynseeds[0] - 1].depth : 2;
            BFSseeds("mt", num_links, num_ctg, ctglinks, ctgdepth, &mt_num_dynseeds, &mt_dynseeds, seq_depth, filter_depth, &mt_bfslinks, &mt_num_BFSlinks);
            int mt_mainseeds_num = 0;
            int* mt_mainseeds = (int*) malloc(sizeof(int) * mt_num_BFSlinks * 2);
            optgfa(exe_path, mt_num_dynseeds, &mt_dynseeds, &mt_bfslinks, &mt_num_BFSlinks, ctgdepth, opts->output_file, 
                    assembly_fna, assembly_graph, "mt", &mt_mainseeds_num, &mt_mainseeds, 0, NULL, 1, filter_depth, cut_seq);
            
            // for (int i = 0; i < mt_num_BFSlinks; i++) {
            //     free(mt_bfslinks[i].lctg); free(mt_bfslinks[i].lutr); free(mt_bfslinks[i].rctg); free(mt_bfslinks[i].rutr);
            // }
            free(mt_bfslinks); 
            free(mt_mainseeds);
        }
        free(mt_dynseeds);        

    } else if (opts->taxo == 2) {
        /* Fungi */
        /* mt */
        int mt_ctg_threshold = 1;
        int mt_num_dynseeds = 0;
        int* mt_dynseeds = NULL;
        if ((2 * seq_depth) > 1.5) {
            filter_depth = 2 * seq_depth;
        } else {
            filter_depth = 1.5;
        }
        mt_dynseeds = calloc(mt_ctg_threshold, sizeof(int));
        HitSeeds(exe_path, "mt", assembly_fna, opts->output_file, opts->cpu, num_ctg, ctgdepth, &mt_dynseeds, &mt_ctg_threshold, 2*seq_depth, 2, 0);
        for (i = 0; i < mt_ctg_threshold; i++) {
            if (mt_dynseeds[i] != 0) {
                mt_num_dynseeds++;
            }
        }
        if (mt_num_dynseeds > 0) {
            BFSlinks* mt_bfslinks = NULL;
            int mt_num_BFSlinks = 0;
            filter_depth = (ctgdepth[mt_dynseeds[0] - 1].depth > 4) ? 0.3 * ctgdepth[mt_dynseeds[0] - 1].depth : 2;
            BFSseeds("mt", num_links, num_ctg, ctglinks, ctgdepth, &mt_num_dynseeds, &mt_dynseeds, seq_depth, filter_depth, &mt_bfslinks, &mt_num_BFSlinks);
            int mt_mainseeds_num = 0;
            int* mt_mainseeds = (int*) malloc(sizeof(int) * mt_num_BFSlinks * 2);
            optgfa(exe_path, mt_num_dynseeds, &mt_dynseeds, &mt_bfslinks, &mt_num_BFSlinks, ctgdepth, opts->output_file, 
                    assembly_fna, assembly_graph, "mt", &mt_mainseeds_num, &mt_mainseeds, 0, NULL, 2, filter_depth, cut_seq);
            
            // for (int i = 0; i < mt_num_BFSlinks; i++) {
            //     free(mt_bfslinks[i].lctg); free(mt_bfslinks[i].lutr); free(mt_bfslinks[i].rctg); free(mt_bfslinks[i].rutr);
            // }
            free(mt_bfslinks); 
            free(mt_mainseeds);
        }
        free(mt_dynseeds);        
    }

    for (i = 0; i < num_ctg; i++) {
        free(ctgdepth[i].ctg);
    }
    free(ctgdepth);

    for (i = 0; i < num_links; i++) {
        free(ctglinks[i].lutr);  free(ctglinks[i].rutr);
    }
    
    free(ctglinks);

    /* free memory */
    if(high_quality_seq!= NULL) free(high_quality_seq);
    if(subsample_seq!= NULL) free(subsample_seq);
    if(cut_seq!= NULL) free(cut_seq);
    if(assembly_fna!= NULL) free(assembly_fna);
    if(assembly_graph!= NULL) free(assembly_graph);
}