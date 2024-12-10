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
#include <libgen.h> // dirname
#include <math.h> // sqrt



#include "log.h"
#include "misc.h"
#include "hitseeds.h"


static void run_blastn(const char *all_contigs, const char *db_path, char *blastn_out, int num_threads, int *num_hits);
int compare_ctg_scores(const void *a, const void *b);

// Conserved genes
const char *conserved_genes[] = {
    "atp1", "atp4", "atp6", "atp8", "atp9", "cob", "cox1", "cox2", "cox3", "mttB", "matR", 
    "rpl2", "rpl10", "rpl16", "rps1", "rps3", "rps4", "rps7", "rps12", "nad1", "nad2", 
    "nad3", "nad4", "nad4L", "nad5", "nad6", "nad7", "nad9", "ccmB", "ccmC", "ccmFn", 
    "ccmFc", "sdh3", "sdh4"
};
size_t gene_count = sizeof(conserved_genes) / sizeof(conserved_genes[0]);

// Exon genes
const char *exon_genes[] = {"rpl2", "rps3", "nad1", "nad2", "nad4", "nad5", "nad7", "ccmFc"};
size_t exon_gene_count = sizeof(exon_genes) / sizeof(exon_genes[0]);

// Conserved genes
const char *fu_conserved_genes[] = {
    "atp6", "atp8", "atp9", "cob", "cox1", "cox2", "cox3", "nad1", 
    "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6", "rps3", "rrnL", "rrnS"
};
const int fu_gene_len[] = {770, 150, 225, 1150, 1520, 750, 810, 1030, 1600, 430, 880, 290, 1950, 630, 1378, 3340, 1720};
size_t fu_gene_count = sizeof(fu_conserved_genes) / sizeof(fu_conserved_genes[0]);

// Exon genes
const char *fu_exon_genes[] = {"cox1", "rrnL", "rrnS"};
size_t fu_exon_gene_count = sizeof(fu_exon_genes) / sizeof(fu_exon_genes[0]);


void MtHitseeds(const char* exe_path, const char* organelles_type, const char* all_contigs, 
             const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, int** candidate_seeds, int* ctg_threshold, float filter_depth) {

    log_message(INFO, "Finding Mt seeds...");
    char *dir = dirname(strdup(exe_path));  // strdup to avoid modifying exe_path
    size_t dir_len = strlen(dir);
    size_t db_path_len = dir_len + strlen("/Conserved_PCGs_db/Plant_conserved_mtgene_nt.fa") + 1;
    char db_path[db_path_len];
    
    if (strcmp(organelles_type, "mt") == 0) {
        snprintf(db_path, db_path_len, "%s/Conserved_PCGs_db/Plant_conserved_mtgene_nt.fa", dir);
    } else {
        log_message(ERROR, "Invalid organelles type %s", organelles_type);
    }

    // log_message(INFO, "The path of runAssembler is %s", sif_path);
    if (access(db_path, F_OK) != 0) {
        log_message(ERROR, "Failed to find the database file %s", db_path);
        free(dir);
        exit(EXIT_FAILURE);
    }
    free(dir);
    FILE *fasta = fopen(db_path, "r");
    if (!fasta) {
        log_message(ERROR, "Failed to open file %s", db_path);
        exit(EXIT_FAILURE);
    }
    
    size_t line_len;
    char *line = NULL;
    int num_pcgs = 100;
    PCGs *pcgs = malloc(sizeof(PCGs) * num_pcgs);
    int idx = 0;
    
    while ((line_len = getline(&line, &line_len, fasta)) != -1) {
        if (line[0] == '>') {
            memset(pcgs[idx].gene, 0, sizeof(pcgs[idx].gene));
            int parsed = sscanf(line, ">%*[^_]_%*[^_]_%s", pcgs[idx].gene);
            if (parsed != 1) {
                log_message(ERROR, "Failed to parse gene name from line %s", line);
                exit(EXIT_FAILURE);
            }
        } else {
            pcgs[idx].length = line_len - 1;
            
            idx += 1;
            if (idx == num_pcgs) {
                num_pcgs *= 2;
                pcgs = realloc(pcgs, sizeof(PCGs) * num_pcgs);
                if (!pcgs) {
                    log_message(ERROR, "Failed to allocate memory for PCGs");
                    exit(EXIT_FAILURE);
                }
            }
        }

    }

    free(line);
    fclose(fasta);

    PCGs *conv_pcgs = malloc(sizeof(PCGs) * gene_count);
    for (size_t i = 0; i < gene_count; i++) {
        const char *gene = conserved_genes[i];
        conv_pcgs[i].length = 0;
        strncpy(conv_pcgs[i].gene, gene, sizeof(conv_pcgs[i].gene));
        int num = 0;

        for (int j = 0; j < idx; j++) {
            if (strcmp(pcgs[j].gene, gene) == 0) {
                num += 1;
                conv_pcgs[i].length += pcgs[j].length;
            }
        }
        
        if (num > 1) {
            /* caculate the average */
            conv_pcgs[i].length /= num;
        }
    }
    free(pcgs);

    // Allocate memory for the output path
    mkdirfiles(output_path);

    size_t blastn_out_len = strlen(output_path) + strlen("/PMAT_mt_blastn.txt") + 1;
    char blastn_out[blastn_out_len];
    snprintf(blastn_out, blastn_out_len, "%s/PMAT_mt_blastn.txt", output_path);
    
    /* Run blastn */
    int num_hits = 0;
    run_blastn(all_contigs, db_path, blastn_out, num_threads, &num_hits); //*
    FILE *blastn_file = fopen(blastn_out, "r");
    if (!blastn_file) {
        log_message(ERROR, "Failed to open file %s", blastn_out);
        exit(EXIT_FAILURE);
    }
    /* PCG; Contig; identity; length */
    BlastInfo *blast_info = calloc(num_hits, sizeof(BlastInfo));
    int blast_idx = 0;

    if (!blast_info) {
        log_message(ERROR, "Failed to allocate memory for BlastInfo");
        exit(EXIT_FAILURE);
    }

    /* PCGs; Contigs; Score */
    PcgCtgs *pcg_ctgs = calloc(gene_count, sizeof(PcgCtgs));
    size_t blt_len;
    char *blt_line = NULL;
    
    int seeds_ctg = 0;
    /* Find the best hit for each PCG */
    while ((blt_len = getline(&blt_line, &blt_len, blastn_file)) != -1) {
        char* tempgene;
        char* organgene;
        char* tempctg;
        float tempiden;
        int templen;
        if (blt_line[0] == '#') {
            continue;
        } else {
            char* token = strtok(blt_line, "\t");
            tempctg = strdup(token);
            token = strtok(NULL, "\t");
            organgene = strdup(token);
            token = strtok(NULL, "\t");
            tempiden = atof(token);
            token = strtok(NULL, "\t");
            templen = atoi(token);
        }

        char *ptr = strchr(organgene, '-');
        if (ptr) {
            *ptr = '\0';
        }

        char* loc = strrchr(organgene, '_');
        tempgene = strdup(loc + 1);

        int tem_flag = 0;
        int tempctg_int = rm_contig(tempctg);
        if (ctg_depth[tempctg_int - 1].depth > filter_depth) {
            for (size_t i = 0; i < gene_count; i++) {
                if (strcmp(tempgene, conv_pcgs[i].gene) == 0) {
                    if (findstr(exon_genes, exon_gene_count, tempgene) == 1 && tempiden > 90 && templen > 500) {

                        if (pcg_ctgs[i].ctg == NULL) {
                            pcg_ctgs[i].ctg = malloc(sizeof(char*)); // malloc for the first contig
                            pcg_ctgs[i].ctg[0] = strdup(tempctg);

                            pcg_ctgs[i].gene = strdup(tempgene);
                            pcg_ctgs[i].num_ctg = 1;

                            pcg_ctgs[i].score = malloc(sizeof(float)); // malloc for the first score
                            pcg_ctgs[i].score[0] = sqrt(sqrt(tempiden * templen) * ctg_depth[tempctg_int - 1].score);

                            pcg_ctgs[i].ctglen = malloc(sizeof(int));
                            pcg_ctgs[i].ctglen[0] = ctg_depth[tempctg_int - 1].len;

                            pcg_ctgs[i].ctgdep = malloc(sizeof(float));
                            pcg_ctgs[i].ctgdep[0] = ctg_depth[tempctg_int - 1].depth;

                        } else {
                            pcg_ctgs[i].ctg = realloc(pcg_ctgs[i].ctg, (pcg_ctgs[i].num_ctg + 1) * sizeof(char*));
                            pcg_ctgs[i].ctg[pcg_ctgs[i].num_ctg] = strdup(tempctg);

                            pcg_ctgs[i].score = realloc(pcg_ctgs[i].score, (pcg_ctgs[i].num_ctg + 1) * sizeof(float));
                            pcg_ctgs[i].score[pcg_ctgs[i].num_ctg] = sqrt(sqrt(tempiden * templen) * ctg_depth[tempctg_int - 1].score);

                            pcg_ctgs[i].ctglen = realloc(pcg_ctgs[i].ctglen, (pcg_ctgs[i].num_ctg + 1) * sizeof(int));
                            pcg_ctgs[i].ctglen[pcg_ctgs[i].num_ctg] = ctg_depth[tempctg_int - 1].len;

                            pcg_ctgs[i].ctgdep = realloc(pcg_ctgs[i].ctgdep, (pcg_ctgs[i].num_ctg + 1) * sizeof(float));                     
                            pcg_ctgs[i].ctgdep[pcg_ctgs[i].num_ctg] = ctg_depth[tempctg_int - 1].depth;

                            pcg_ctgs[i].num_ctg += 1;
                        }
                        seeds_ctg += 1;
                        blast_idx += 1;
                        break;

                    } else if (findstr(exon_genes, exon_gene_count, tempgene) == 0 && tempiden > 90 && templen > 0.85*conv_pcgs[i].length) {

                        if (pcg_ctgs[i].ctg == NULL) {
                            pcg_ctgs[i].ctg = malloc(sizeof(char*)); // malloc for the first contig
                            pcg_ctgs[i].ctg[0] = strdup(tempctg);

                            pcg_ctgs[i].gene = strdup(tempgene);
                            pcg_ctgs[i].num_ctg = 1;
                            
                            pcg_ctgs[i].score = malloc(sizeof(float)); // malloc for the first score
                            pcg_ctgs[i].score[0] = sqrt(sqrt(tempiden * templen) * ctg_depth[tempctg_int - 1].score);
                            
                            pcg_ctgs[i].ctglen = malloc(sizeof(int));
                            pcg_ctgs[i].ctglen[0] = ctg_depth[tempctg_int - 1].len;
                            
                            pcg_ctgs[i].ctgdep = malloc(sizeof(float));
                            pcg_ctgs[i].ctgdep[0] = ctg_depth[tempctg_int - 1].depth;
                            seeds_ctg += 1;
                        } else if (pcg_ctgs[i].score[0] < sqrt(tempiden * templen)) {
                            free(pcg_ctgs[i].ctg[0]);
                            pcg_ctgs[i].ctg[0] = strdup(tempctg);
                            pcg_ctgs[i].score[0] = sqrt(sqrt(tempiden * templen) * ctg_depth[tempctg_int - 1].score);
                            pcg_ctgs[i].ctglen[0] = ctg_depth[tempctg_int - 1].len;
                            pcg_ctgs[i].ctgdep[0] = ctg_depth[tempctg_int - 1].depth;
                            seeds_ctg += 1;
                        }
                        blast_idx += 1;
                        break;
                    }
                }
            }
        }
    }
    
    /* Sort the scores in descending order */
    SortPcgCtgs *sort_pcg_ctgs = malloc(sizeof(SortPcgCtgs) * 30);
    size_t sort_idx = 0;
    for (size_t i = 0; i < gene_count; i++) {
        if (sort_idx >= 30) break;
        if (pcg_ctgs[i].ctg != NULL) {
            for (int j = 0; j < pcg_ctgs[i].num_ctg; j++) {
                sort_pcg_ctgs[sort_idx].ctg = pcg_ctgs[i].ctg[j];
                sort_pcg_ctgs[sort_idx].score = pcg_ctgs[i].score[j];
                sort_pcg_ctgs[sort_idx].ctglen = pcg_ctgs[i].ctglen[j];
                sort_pcg_ctgs[sort_idx].ctgdep = pcg_ctgs[i].ctgdep[j];
                sort_idx += 1;
                if (sort_idx >= 30) break;
            }
        }
    }


    if (sort_idx == 0) {
        log_message(WARNING, "No seed contigs found (mt), please use GraphBuild command.");
    } else {
        qsort(sort_pcg_ctgs, sort_idx, sizeof(SortPcgCtgs), compare_ctg_scores);
        *ctg_threshold = sort_idx;
        *candidate_seeds = calloc(sort_idx, sizeof(int));

        log_message(INFO, "Seed finding process is complete.");
        log_info(" _______________________________________________________\n");
        log_info(" Contig Name    Length (bp)   Depth (x)     Score    \n");
        log_info(" -------------  ------------  ------------  ------------\n");
        for (size_t i = 0; i < sort_idx; i++) {
            log_info(" %-12s   %-12d  %-12g  %-10.2f\n", 
                    sort_pcg_ctgs[i].ctg, 
                    sort_pcg_ctgs[i].ctglen, 
                    sort_pcg_ctgs[i].ctgdep, 
                    sort_pcg_ctgs[i].score
                    );
                    (*candidate_seeds)[i] = rm_contig(sort_pcg_ctgs[i].ctg);
        }
        log_info(" _______________________________________________________\n");
        log_info("\n");
    }

    free(sort_pcg_ctgs); free(blt_line); free(pcg_ctgs); fclose(blastn_file); free(blast_info); free(conv_pcgs);
}



void PtHitseeds(const char* exe_path, const char* organelles_type, const char* all_contigs, 
             const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, int* candidate_seeds, int ctg_threshold, float filter_depth) {

    log_message(INFO, "Finding Pt seeds...");
    char *dir = dirname(strdup(exe_path));  // strdup to avoid modifying exe_path
    size_t dir_len = strlen(dir);
    size_t db_path_len = dir_len + strlen("/Conserved_PCGs_db/Plant_conserved_cpgene_nt.fa") + 1;
    char db_path[db_path_len];
    
    if (strcmp(organelles_type, "pt") == 0) {
        snprintf(db_path, db_path_len, "%s/Conserved_PCGs_db/Plant_conserved_cpgene_nt.fa", dir);
    } else  {
        log_message(ERROR, "Invalid organelles type: %s", organelles_type);
    }

    free(dir);


    // Allocate memory for the output path
    mkdirfiles(output_path);

    size_t blastn_out_len = strlen(output_path) + strlen("/PMAT_pt_blastn.txt") + 1;
    char blastn_out[blastn_out_len];
    snprintf(blastn_out, blastn_out_len, "%s/PMAT_pt_blastn.txt", output_path);
    
    /* Run blastn */
    int num_hits = 0;
    run_blastn(all_contigs, db_path, blastn_out, num_threads, &num_hits); //*
    FILE *blastn_file = fopen(blastn_out, "r");
    if (!blastn_file) {
        log_message(ERROR, "Failed to open file %s", blastn_out);
        exit(EXIT_FAILURE);
    }
    /* PCG; Contig; identity; length */
    BlastInfo *blast_info = calloc(num_hits, sizeof(BlastInfo));
    int blast_idx = 0;

    if (!blast_info) {
        log_message(ERROR, "Failed to allocate memory for BlastInfo");
        exit(EXIT_FAILURE);
    }

    /* PCGs; Contigs; Score */
    SortPcgCtgs *ptpcg_ctgs = NULL;
    size_t blt_len;
    char *blt_line = NULL;
    
    int seeds_ctg = 0;
    /* Find the best hit for each PCG */
    while ((blt_len = getline(&blt_line, &blt_len, blastn_file)) != -1) {
        char* organgene;
        char* tempctg;
        int inttempctg;
        float tempiden;
        int templen;
        if (blt_line[0] == '#') {
            continue;
        } else {
            char* token = strtok(blt_line, "\t");
            tempctg = strdup(token);
            inttempctg = rm_contig(tempctg);
            token = strtok(NULL, "\t");
            organgene = strdup(token);
            token = strtok(NULL, "\t");
            tempiden = atof(token);
            token = strtok(NULL, "\t");
            templen = atoi(token);
        }

        int tem_flag = 0;
        int tempctg_int = rm_contig(tempctg);
        if (ctg_depth[tempctg_int - 1].depth > filter_depth) {
            if (tempiden > 70 && templen > 500) {
                ptpcg_ctgs = realloc(ptpcg_ctgs, (seeds_ctg + 1) * sizeof(PcgCtgs));
                ptpcg_ctgs[seeds_ctg].ctg = strdup(tempctg);
                ptpcg_ctgs[seeds_ctg].score = sqrt((ctg_depth[inttempctg - 1].depth)*sqrt(sqrt(ctg_depth[inttempctg - 1].len * tempiden * templen)));
                ptpcg_ctgs[seeds_ctg].ctglen = ctg_depth[inttempctg - 1].len;
                ptpcg_ctgs[seeds_ctg].ctgdep = ctg_depth[inttempctg - 1].depth;
                seeds_ctg += 1;
            }
        }
    }
    

    if (seeds_ctg == 0) {
        log_message(WARNING, "No seed contigs found (pt), please use GraphBuild command.");
    } else {
        qsort(ptpcg_ctgs, seeds_ctg, sizeof(SortPcgCtgs), compare_ctg_scores);


        log_message(INFO, "Seed finding process is complete.");
        log_info(" _______________________________________________________\n");
        log_info(" Contig Name    Length (bp)   Depth (x)     Score    \n");
        log_info(" -------------  ------------  ------------  ------------\n");
        for (size_t i = 0; i < seeds_ctg; i++) {
            log_info(" %-12s   %-12d  %-12g  %-10.2f\n", 
                    ptpcg_ctgs[i].ctg, 
                    ptpcg_ctgs[i].ctglen, 
                    ptpcg_ctgs[i].ctgdep, 
                    ptpcg_ctgs[i].score
                    );
            if (i < ctg_threshold) {
                candidate_seeds[i] = rm_contig(ptpcg_ctgs[i].ctg);
            }
        }
        log_info(" _______________________________________________________\n");
        log_info("\n");
    }

    free(ptpcg_ctgs);
    free(blt_line);
    fclose(blastn_file);
    free(blast_info);
}

void AnHitseeds(const char* exe_path, const char* organelles_type, const char* all_contigs, 
             const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, int** candidate_seeds, int* ctg_threshold, float filter_depth) {

    log_message(INFO, "Finding Mt seeds...");
    char *dir = dirname(strdup(exe_path));  // strdup to avoid modifying exe_path
    size_t dir_len = strlen(dir);
    size_t db_path_len = dir_len + strlen("/Conserved_PCGs_db/Animal_conserved_cpgene_nt.fa") + 1;
    char db_path[db_path_len];
    
    if (strcmp(organelles_type, "mt") == 0) {
        snprintf(db_path, db_path_len, "%s/Conserved_PCGs_db/Animal_conserved_mtgene_nt.fa", dir);
    } else  {
        log_message(ERROR, "Invalid organelles type: %s", organelles_type);
    }

    free(dir);

    // Allocate memory for the output path
    mkdirfiles(output_path);

    size_t blastn_out_len = strlen(output_path) + strlen("/PMAT_mt_blastn.txt") + 1;
    char blastn_out[blastn_out_len];
    snprintf(blastn_out, blastn_out_len, "%s/PMAT_mt_blastn.txt", output_path);
    
    /* Run blastn */
    int num_hits = 0;
    run_blastn(all_contigs, db_path, blastn_out, num_threads, &num_hits); //*
    FILE *blastn_file = fopen(blastn_out, "r");
    if (!blastn_file) {
        log_message(ERROR, "Failed to open file %s", blastn_out);
        exit(EXIT_FAILURE);
    }

    BlastInfo *blast_info = calloc(num_hits, sizeof(BlastInfo));
    int blast_idx = 0;

    if (!blast_info) {
        log_message(ERROR, "Failed to allocate memory for BlastInfo");
        exit(EXIT_FAILURE);
    }

    SortPcgCtgs *mtpcg_ctgs = NULL;
    size_t blt_len;
    char *blt_line = NULL;
    
    int seeds_ctg = 0;
    /* Find the best hit for each PCG */
    while ((blt_len = getline(&blt_line, &blt_len, blastn_file)) != -1) {
        char* organgene;
        char* tempctg;
        int inttempctg;
        float tempiden;
        int templen;
        if (blt_line[0] == '#') {
            continue;
        } else {
            char* token = strtok(blt_line, "\t");
            tempctg = strdup(token);
            inttempctg = rm_contig(tempctg);
            token = strtok(NULL, "\t");
            organgene = strdup(token);
            token = strtok(NULL, "\t");
            tempiden = atof(token);
            token = strtok(NULL, "\t");
            templen = atoi(token);
        }

        int tem_flag = 0;
        int tempctg_int = rm_contig(tempctg);
        if (ctg_depth[tempctg_int - 1].depth > filter_depth) {
            if (tempiden > 70 && templen > 500) {
                mtpcg_ctgs = realloc(mtpcg_ctgs, (seeds_ctg + 1) * sizeof(PcgCtgs));
                mtpcg_ctgs[seeds_ctg].ctg = strdup(tempctg);
                mtpcg_ctgs[seeds_ctg].score = sqrt((ctg_depth[inttempctg - 1].depth)*sqrt(sqrt(ctg_depth[inttempctg - 1].len * tempiden * templen)));
                mtpcg_ctgs[seeds_ctg].ctglen = ctg_depth[inttempctg - 1].len;
                mtpcg_ctgs[seeds_ctg].ctgdep = ctg_depth[inttempctg - 1].depth;
                seeds_ctg += 1;
            }
        }
    }
    

    if (seeds_ctg == 0) {
        log_message(WARNING, "No seed contigs found (mt), please use GraphBuild command.");
    } else {
        qsort(mtpcg_ctgs, seeds_ctg, sizeof(SortPcgCtgs), compare_ctg_scores);
        *candidate_seeds = calloc(seeds_ctg, sizeof(int));
        int hit_depth = ctg_depth[rm_contig(mtpcg_ctgs[0].ctg) - 1].depth;

        log_message(INFO, "Seed finding process is complete.");
        log_info(" _______________________________________________________\n");
        log_info(" Contig Name    Length (bp)   Depth (x)     Score    \n");
        log_info(" -------------  ------------  ------------  ------------\n");
        (*ctg_threshold) = 0;
        for (size_t i = 0; i < seeds_ctg; i++) {
            log_info(" %-12s   %-12d  %-12g  %-10.2f\n", 
                    mtpcg_ctgs[i].ctg, 
                    mtpcg_ctgs[i].ctglen, 
                    mtpcg_ctgs[i].ctgdep, 
                    mtpcg_ctgs[i].score
                    );
            if (ctg_depth[rm_contig(mtpcg_ctgs[i].ctg) - 1].depth > 0.3*hit_depth && ctg_depth[rm_contig(mtpcg_ctgs[i].ctg) - 1].depth < 3*hit_depth) {
                (*candidate_seeds)[i] = rm_contig(mtpcg_ctgs[i].ctg);
                (*ctg_threshold)++;
            }
        }
        log_info(" _______________________________________________________\n");
        log_info("\n");
    }

    free(mtpcg_ctgs);
    free(blt_line);
    fclose(blastn_file);
    free(blast_info);
}


// void FuHitseeds(const char* exe_path, const char* organelles_type, const char* all_contigs, 
//              const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, int** candidate_seeds, int* ctg_threshold, float filter_depth) {

//     log_message(INFO, "Finding Mt seeds...");
//     char *dir = dirname(strdup(exe_path));  // strdup to avoid modifying exe_path
//     size_t dir_len = strlen(dir);
//     size_t db_path_len = dir_len + strlen("/Conserved_PCGs_db/Fungi_conserved_cpgene_nt.fa") + 1;
//     char db_path[db_path_len];
    
//     if (strcmp(organelles_type, "mt") == 0) {
//         snprintf(db_path, db_path_len, "%s/Conserved_PCGs_db/Fungi_conserved_mtgene_nt.fa", dir);
//     } else  {
//         log_message(ERROR, "Invalid organelles type: %s", organelles_type);
//     }

//     free(dir);

//     // Allocate memory for the output path
//     mkdirfiles(output_path);

//     size_t blastn_out_len = strlen(output_path) + strlen("/PMAT_mt_blastn.txt") + 1;
//     char blastn_out[blastn_out_len];
//     snprintf(blastn_out, blastn_out_len, "%s/PMAT_mt_blastn.txt", output_path);
    
//     /* Run blastn */
//     int num_hits = 0;
//     run_blastn(all_contigs, db_path, blastn_out, num_threads, &num_hits); //*
//     FILE *blastn_file = fopen(blastn_out, "r");
//     if (!blastn_file) {
//         log_message(ERROR, "Failed to open file %s", blastn_out);
//         exit(EXIT_FAILURE);
//     }

//     BlastInfo *blast_info = calloc(num_hits, sizeof(BlastInfo));
//     int blast_idx = 0;

//     if (!blast_info) {
//         log_message(ERROR, "Failed to allocate memory for BlastInfo");
//         exit(EXIT_FAILURE);
//     }

//     SortPcgCtgs *mtpcg_ctgs = NULL;
//     size_t blt_len;
//     char *blt_line = NULL;
    
//     int seeds_ctg = 0;
//     /* Find the best hit for each PCG */
//     while ((blt_len = getline(&blt_line, &blt_len, blastn_file)) != -1) {
//         char* organgene;
//         char* tempctg;
//         int inttempctg;
//         float tempiden;
//         int templen;
//         if (blt_line[0] == '#') {
//             continue;
//         } else {
//             char* token = strtok(blt_line, "\t");
//             tempctg = strdup(token);
//             inttempctg = rm_contig(tempctg);
//             token = strtok(NULL, "\t");
//             organgene = strdup(token);
//             token = strtok(NULL, "\t");
//             tempiden = atof(token);
//             token = strtok(NULL, "\t");
//             templen = atoi(token);
//         }

//         int tem_flag = 0;
//         int tempctg_int = rm_contig(tempctg);
//         if (ctg_depth[tempctg_int - 1].depth > filter_depth) {
//             if (tempiden > 80 && templen > 300) {
//                 mtpcg_ctgs = realloc(mtpcg_ctgs, (seeds_ctg + 1) * sizeof(PcgCtgs));
//                 mtpcg_ctgs[seeds_ctg].ctg = strdup(tempctg);
//                 mtpcg_ctgs[seeds_ctg].score = sqrt((ctg_depth[inttempctg - 1].depth)*sqrt(sqrt(ctg_depth[inttempctg - 1].len * tempiden * templen)));
//                 mtpcg_ctgs[seeds_ctg].ctglen = ctg_depth[inttempctg - 1].len;
//                 mtpcg_ctgs[seeds_ctg].ctgdep = ctg_depth[inttempctg - 1].depth;
//                 seeds_ctg += 1;
//             }
//         }
//     }
    

//     if (seeds_ctg == 0) {
//         log_message(WARNING, "No seed contigs found (mt), please use GraphBuild command.");
//     } else {
//         qsort(mtpcg_ctgs, seeds_ctg, sizeof(SortPcgCtgs), compare_ctg_scores);
//         *candidate_seeds = calloc(seeds_ctg, sizeof(int));
//         int hit_depth = ctg_depth[rm_contig(mtpcg_ctgs[0].ctg) - 1].depth;

//         log_message(INFO, "Seed finding process is complete.");
//         log_info(" _______________________________________________________\n");
//         log_info(" Contig Name    Length (bp)   Depth (x)     Score    \n");
//         log_info(" -------------  ------------  ------------  ------------\n");
//         (*ctg_threshold) = 0;
//         for (size_t i = 0; i < seeds_ctg; i++) {
//             log_info(" %-12s   %-12d  %-12g  %-10.2f\n", 
//                     mtpcg_ctgs[i].ctg, 
//                     mtpcg_ctgs[i].ctglen, 
//                     mtpcg_ctgs[i].ctgdep, 
//                     mtpcg_ctgs[i].score
//                     );
//             if (ctg_depth[rm_contig(mtpcg_ctgs[i].ctg) - 1].depth > 0.3*hit_depth && ctg_depth[rm_contig(mtpcg_ctgs[i].ctg) - 1].depth < 3*hit_depth) {
//                 (*candidate_seeds)[i] = rm_contig(mtpcg_ctgs[i].ctg);
//                 (*ctg_threshold)++;
//             }
//         }
//         log_info(" _______________________________________________________\n");
//         log_info("\n");
//     }

//     free(mtpcg_ctgs);
//     free(blt_line);
//     fclose(blastn_file);
//     free(blast_info);
// }


void FuHitseeds(const char* exe_path, const char* organelles_type, const char* all_contigs, 
             const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, int** candidate_seeds, int* ctg_threshold, float filter_depth) {

    log_message(INFO, "Finding Mt seeds...");
    char *dir = dirname(strdup(exe_path));  // strdup to avoid modifying exe_path
    size_t dir_len = strlen(dir);
    size_t db_path_len = dir_len + strlen("/Conserved_PCGs_db/Fungi_conserved_mtgene_nt.fa") + 1;
    char db_path[db_path_len];
    
    if (strcmp(organelles_type, "mt") == 0) {
        snprintf(db_path, db_path_len, "%s/Conserved_PCGs_db/Fungi_conserved_mtgene_nt.fa", dir);
    } else {
        log_message(ERROR, "Invalid organelles type %s", organelles_type);
    }

    PCGs *conv_pcgs = malloc(sizeof(PCGs) * fu_gene_count);
    for (size_t i = 0; i < fu_gene_count; i++) {
        strncpy(conv_pcgs[i].gene, fu_conserved_genes[i], sizeof(conv_pcgs[i].gene));
        conv_pcgs[i].length = fu_gene_len[i];
    }

    // Allocate memory for the output path
    mkdirfiles(output_path);

    size_t blastn_out_len = strlen(output_path) + strlen("/PMAT_mt_blastn.txt") + 1;
    char blastn_out[blastn_out_len];
    snprintf(blastn_out, blastn_out_len, "%s/PMAT_mt_blastn.txt", output_path);
    
    /* Run blastn */
    int num_hits = 0;
    run_blastn(all_contigs, db_path, blastn_out, num_threads, &num_hits); //*
    FILE *blastn_file = fopen(blastn_out, "r");
    if (!blastn_file) {
        log_message(ERROR, "Failed to open file %s", blastn_out);
        exit(EXIT_FAILURE);
    }
    /* PCG; Contig; identity; length */
    BlastInfo *blast_info = calloc(num_hits, sizeof(BlastInfo));
    int blast_idx = 0;

    if (!blast_info) {
        log_message(ERROR, "Failed to allocate memory for BlastInfo");
        exit(EXIT_FAILURE);
    }

    /* PCGs; Contigs; Score */
    PcgCtgs *pcg_ctgs = calloc(fu_gene_count, sizeof(PcgCtgs));
    size_t blt_len;
    char *blt_line = NULL;
    
    int seeds_ctg = 0;
    /* Find the best hit for each PCG */
    while ((blt_len = getline(&blt_line, &blt_len, blastn_file)) != -1) {
        char* tempgene;
        char* organgene;
        char* tempctg;
        float tempiden;
        int templen;
        if (blt_line[0] == '#') {
            continue;
        } else {
            char* token = strtok(blt_line, "\t");
            tempctg = strdup(token);
            token = strtok(NULL, "\t");
            organgene = strdup(token);
            token = strtok(NULL, "\t");
            tempiden = atof(token);
            token = strtok(NULL, "\t");
            templen = atoi(token);
        }

        char *ptr = strchr(organgene, '-');
        if (ptr) {
            *ptr = '\0';
        }

        char* loc = strrchr(organgene, '_');
        tempgene = strdup(loc + 1);

        int tem_flag = 0;
        int tempctg_int = rm_contig(tempctg);
        if (ctg_depth[tempctg_int - 1].depth > filter_depth) {
            for (size_t i = 0; i < fu_gene_count; i++) {
                if (strcmp(tempgene, conv_pcgs[i].gene) == 0) {
                    if (findstr(fu_exon_genes, fu_exon_gene_count, tempgene) == 1 && tempiden > 70 && templen > 300) {

                        if (pcg_ctgs[i].ctg == NULL) {
                            pcg_ctgs[i].ctg = malloc(sizeof(char*)); // malloc for the first contig
                            pcg_ctgs[i].ctg[0] = strdup(tempctg);

                            pcg_ctgs[i].gene = strdup(tempgene);
                            pcg_ctgs[i].num_ctg = 1;

                            pcg_ctgs[i].score = malloc(sizeof(float)); // malloc for the first score
                            pcg_ctgs[i].score[0] = sqrt(sqrt(tempiden * templen) * ctg_depth[tempctg_int - 1].score);

                            pcg_ctgs[i].ctglen = malloc(sizeof(int));
                            pcg_ctgs[i].ctglen[0] = ctg_depth[tempctg_int - 1].len;

                            pcg_ctgs[i].ctgdep = malloc(sizeof(float));
                            pcg_ctgs[i].ctgdep[0] = ctg_depth[tempctg_int - 1].depth;

                        } else {
                            pcg_ctgs[i].ctg = realloc(pcg_ctgs[i].ctg, (pcg_ctgs[i].num_ctg + 1) * sizeof(char*));
                            pcg_ctgs[i].ctg[pcg_ctgs[i].num_ctg] = strdup(tempctg);

                            pcg_ctgs[i].score = realloc(pcg_ctgs[i].score, (pcg_ctgs[i].num_ctg + 1) * sizeof(float));
                            pcg_ctgs[i].score[pcg_ctgs[i].num_ctg] = sqrt(sqrt(tempiden * templen) * ctg_depth[tempctg_int - 1].score);

                            pcg_ctgs[i].ctglen = realloc(pcg_ctgs[i].ctglen, (pcg_ctgs[i].num_ctg + 1) * sizeof(int));
                            pcg_ctgs[i].ctglen[pcg_ctgs[i].num_ctg] = ctg_depth[tempctg_int - 1].len;

                            pcg_ctgs[i].ctgdep = realloc(pcg_ctgs[i].ctgdep, (pcg_ctgs[i].num_ctg + 1) * sizeof(float));                     
                            pcg_ctgs[i].ctgdep[pcg_ctgs[i].num_ctg] = ctg_depth[tempctg_int - 1].depth;

                            pcg_ctgs[i].num_ctg += 1;
                        }
                        seeds_ctg += 1;
                        blast_idx += 1;
                        break;

                    } else if (findstr(fu_exon_genes, fu_exon_gene_count, tempgene) == 0 && tempiden > 70 && templen > 0.5*conv_pcgs[i].length) {

                        if (pcg_ctgs[i].ctg == NULL) {
                            pcg_ctgs[i].ctg = malloc(sizeof(char*)); // malloc for the first contig
                            pcg_ctgs[i].ctg[0] = strdup(tempctg);

                            pcg_ctgs[i].gene = strdup(tempgene);
                            pcg_ctgs[i].num_ctg = 1;
                            
                            pcg_ctgs[i].score = malloc(sizeof(float)); // malloc for the first score
                            pcg_ctgs[i].score[0] = sqrt(sqrt(tempiden * templen) * ctg_depth[tempctg_int - 1].score);
                            
                            pcg_ctgs[i].ctglen = malloc(sizeof(int));
                            pcg_ctgs[i].ctglen[0] = ctg_depth[tempctg_int - 1].len;
                            
                            pcg_ctgs[i].ctgdep = malloc(sizeof(float));
                            pcg_ctgs[i].ctgdep[0] = ctg_depth[tempctg_int - 1].depth;
                            seeds_ctg += 1;
                        } else if (pcg_ctgs[i].score[0] < sqrt(tempiden * templen)) {
                            free(pcg_ctgs[i].ctg[0]);
                            pcg_ctgs[i].ctg[0] = strdup(tempctg);
                            pcg_ctgs[i].score[0] = sqrt(sqrt(tempiden * templen) * ctg_depth[tempctg_int - 1].score);
                            pcg_ctgs[i].ctglen[0] = ctg_depth[tempctg_int - 1].len;
                            pcg_ctgs[i].ctgdep[0] = ctg_depth[tempctg_int - 1].depth;
                            seeds_ctg += 1;
                        }
                        blast_idx += 1;
                        break;
                    }
                }
            }
        }
    }
    
    /* Sort the scores in descending order */
    SortPcgCtgs *sort_pcg_ctgs = malloc(sizeof(SortPcgCtgs) * 30);
    size_t sort_idx = 0;
    for (size_t i = 0; i < fu_gene_count; i++) {
        if (sort_idx >= 30) break;
        if (pcg_ctgs[i].ctg != NULL) {
            for (int j = 0; j < pcg_ctgs[i].num_ctg; j++) {
                sort_pcg_ctgs[sort_idx].ctg = pcg_ctgs[i].ctg[j];
                sort_pcg_ctgs[sort_idx].score = pcg_ctgs[i].score[j];
                sort_pcg_ctgs[sort_idx].ctglen = pcg_ctgs[i].ctglen[j];
                sort_pcg_ctgs[sort_idx].ctgdep = pcg_ctgs[i].ctgdep[j];
                sort_idx += 1;
                if (sort_idx >= 30) break;
            }
        }
    }


    if (sort_idx == 0) {
        log_message(WARNING, "No seed contigs found (mt), please use GraphBuild command.");
    } else {
        qsort(sort_pcg_ctgs, sort_idx, sizeof(SortPcgCtgs), compare_ctg_scores);
        *ctg_threshold = sort_idx;
        *candidate_seeds = calloc(sort_idx, sizeof(int));

        log_message(INFO, "Seed finding process is complete.");
        log_info(" _______________________________________________________\n");
        log_info(" Contig Name    Length (bp)   Depth (x)     Score    \n");
        log_info(" -------------  ------------  ------------  ------------\n");
        for (size_t i = 0; i < sort_idx; i++) {
            log_info(" %-12s   %-12d  %-12g  %-10.2f\n", 
                    sort_pcg_ctgs[i].ctg, 
                    sort_pcg_ctgs[i].ctglen, 
                    sort_pcg_ctgs[i].ctgdep, 
                    sort_pcg_ctgs[i].score
                    );
                    (*candidate_seeds)[i] = rm_contig(sort_pcg_ctgs[i].ctg);
        }
        log_info(" _______________________________________________________\n");
        log_info("\n");
    }

    free(sort_pcg_ctgs); free(blt_line); free(pcg_ctgs); fclose(blastn_file); free(blast_info); free(conv_pcgs);
}

/* sort the scores in descending order */
int compare_ctg_scores(const void *a, const void *b) {
    float score_a = ((SortPcgCtgs *)a)->score;
    float score_b = ((SortPcgCtgs *)b)->score;
    
    if (score_a > score_b) return -1;
    if (score_a < score_b) return 1;
    return 0;
}

static void run_blastn(const char *all_contigs, const char *db_path, char *blastn_out, int num_threads, int *num_hits) {
    size_t comond_len = snprintf(NULL, 0, "blastn -db %s -query %s -outfmt 6 -num_threads %d -num_alignments 1 -max_hsps 1 > %s", 
                                db_path, all_contigs, num_threads, blastn_out) + 1;
    char* command = malloc(comond_len);
    snprintf(command, comond_len, "blastn -db %s -query %s -outfmt 6 -num_threads %d -num_alignments 1 -max_hsps 1 > %s", db_path, all_contigs, num_threads, blastn_out);

    execute_command(command, 0, 1);

    if (access(blastn_out, F_OK) != 0) {
        log_message(ERROR, "Failed to run blastn");
        free(command);
        return;
    }

    free(command);

    size_t line_len;
    char *line = NULL;
    FILE *blastn_file = fopen(blastn_out, "r");
    if (!blastn_file) {
        log_message(ERROR, "Failed to open file %s", blastn_out);
        return;
    }

    while ((line_len = getline(&line, &line_len, blastn_file)) != -1) {
        if (line[0] == '#') {
            continue;
        } else {
            (*num_hits) += 1;
        }
    }

    fclose(blastn_file);
    free(line);
}

#ifndef HITSEEDS_MAIN
int main(int argc, char **argv) {
    if (argc != 6) {
        log_message(ERROR, "Usage: %s <organelles_type> <all_contigs> <all_graph> <num_threads> <output_path>", argv[0]);
        return -1;
    }

    const char* all_graph = argv[3];
    if (access(all_graph, F_OK) != 0) {
        log_message(ERROR, "Failed to open file %s", all_graph);
        return -1;    
    }


    FILE *graph_file = fopen(all_graph, "r");

    size_t line_len;
    char *line = NULL;
    int num_ctgs = 0;

    while ((line_len = getline(&line, &line_len, graph_file)) != -1) {
        if (line[0] != 'C') {
            num_ctgs += 1;
        } else {
            break;
        }
    }

    CtgDepth *ctg_depth = malloc(sizeof(CtgDepth) * num_ctgs);
    if (!ctg_depth) {
        log_message(ERROR, "Failed to allocate memory for CtgDepth");
        return -1;
    }

    rewind(graph_file);
    int idx = 0;
    while ((line_len = getline(&line, &line_len, graph_file)) != -1) {
        if (line[0] != 'C') {
            char *token = strtok(line, "\t");
            int ctg_id = atoi(token);
            ctg_depth[ctg_id - 1].ctgsmp = ctg_id;
            token = strtok(NULL, "\t");
            ctg_depth[ctg_id - 1].ctg = strdup(token);
            token = strtok(NULL, "\t");
            ctg_depth[ctg_id - 1].len = atoi(token);
            token = strtok(NULL, "\t");
            ctg_depth[ctg_id - 1].depth = atof(token);
            ctg_depth[ctg_id - 1].score = sqrt(sqrt(ctg_depth[idx].depth) * ctg_depth[idx].len);
            idx += 1;
        } else {
            break;
        }
    }
    
    // for (int i = 0; i < idx; i++) {
    //     log_info("Contig: %s, Length: %d, Depth: %g, Score: %g\n", ctg_depth[i].ctg, ctg_depth[i].len, ctg_depth[i].depth, ctg_depth[i].score);
    // }

    free(line);
    fclose(graph_file);

    char *exe_path = realpath(argv[0], NULL);
    char **candidate_seeds = calloc(6, sizeof(char*));
    hitseeds(exe_path, argv[1], argv[2], argv[5], atoi(argv[4]), num_ctgs, ctg_depth, &candidate_seeds, 6);
    for (int i = 0; i < 6; i++) {
        if (candidate_seeds[i] != NULL) {
            printf("%s\n", candidate_seeds[i]);
        }
    }

    /* Free memory */
    for (int i = 0; i < 6; i++) {
        if (candidate_seeds[i] != NULL) {
            free(candidate_seeds[i]);
        }
    }
    for (int i = 0; i < num_ctgs; i++) {
        free(ctg_depth[i].ctg);
    }
    
    free(exe_path);
    free(ctg_depth);
    free(candidate_seeds);
    return 0;
}
#endif