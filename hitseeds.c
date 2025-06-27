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
#include <sys/stat.h> // stat

#include "log.h"
#include "misc.h"
#include "hitseeds.h"
#include "khash.h"


static void run_blastn(const char *all_contigs, const char *db_path, char *blastn_out, int num_threads, int *num_hits);
void mrun_blastn(const char *all_contigs, const char *db_path, char *final_out, int num_threads, int *num_hits);
int compare_ctg_scores(const void *a, const void *b);


static float calc_score_plant(float depth, int len, int num_genes);
static float calc_score_animal(float depth, int len, int num_genes);
static float calc_score_fungi(float depth, int len, int num_genes);


static float sigmoid_score(float score) {
    float sigmoid = 1.0/(1.0 + exp(-score/100));
    return sigmoid * 10;
}

static float calc_score_plant(float depth, int len, int num_genes) {
    float base_score = sqrt(depth * sqrt(len)); 
    float norm_score = sigmoid_score(base_score);
    return norm_score * num_genes * num_genes;
}

static float calc_score_animal(float depth, int len, int num_genes) {
    float base_score = sqrt(depth * sqrt(len)); 
    float norm_score = sigmoid_score(base_score);
    return norm_score * num_genes * 2;
}

static float calc_score_fungi(float depth, int len, int num_genes) {
    float base_score = sqrt(depth * sqrt(len));
    float norm_score = sigmoid_score(base_score);
    return norm_score * num_genes * num_genes;
}

void HitSeeds(const char* exe_path, const char* organelles_type, const char* all_contigs,
             const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, 
             int** candidate_seeds, int* ctg_threshold, float filter_depth,
             int taxo, int verbose) {

    const char *db_suffix;
    const char **mtpcg;
    const int *mtpcg_len;
    int mtpcg_num;
    
    switch(taxo) {
        case 0: // Plant
            db_suffix = "/Conserved_PCGs_db/Plant_conserved_mtgene_nt.fa";
            mtpcg = plt_mtpcg;
            mtpcg_len = plt_mtpcg_len;
            mtpcg_num = plt_mtpcg_num;
            break;
        case 1: // Animal  
            db_suffix = "/Conserved_PCGs_db/Animal_conserved_mtgene_nt.fa";
            mtpcg = anl_mtpcg;
            mtpcg_len = anl_mtpcg_len;
            mtpcg_num = anl_mtpcg_num;
            break;
        case 2: // Fungi
            db_suffix = "/Conserved_PCGs_db/Fungi_conserved_mtgene_nt.fa";
            mtpcg = fug_mtpcg;
            mtpcg_len = fug_mtpcg_len; 
            mtpcg_num = fug_mtpcg_num;
            break;
        default:
            log_message(ERROR, "Invalid taxo type: %d", taxo);
            return;
    }

    log_message(INFO, "Finding Mt seeds...");
    char *dir = dirname(strdup(exe_path));
    char* db_path = (char*)malloc(sizeof(*db_path) * (snprintf(NULL, 0, "%s%s", dir, db_suffix) + 1));
    sprintf(db_path, "%s%s", dir, db_suffix);
    free(dir);

    uint64_t i;
    khint_t k;

    mkdirfiles(output_path);

    char* blastn_out = (char*)malloc(sizeof(*blastn_out) * (snprintf(NULL, 0, "%s/PMAT_mt_blastn.txt", output_path) + 1));
    sprintf(blastn_out, "%s/PMAT_mt_blastn.txt", output_path);
    
    /* Run blastn */
    int num_hits = 0;
    mrun_blastn(all_contigs, db_path, blastn_out, num_threads, &num_hits); //* 
    
    khash_t(ctg_genes) *h_ctg_genes = kh_init(ctg_genes);
    
    FILE *blastn_file = fopen(blastn_out, "r");
    if (!blastn_file) {
        log_message(ERROR, "Failed to open file %s", blastn_out);
        exit(EXIT_FAILURE);
    }

    char *line = NULL;
    size_t len = 0;
    int best_ctg = 0;
    int best_gnum = 0;

    while (getline(&line, &len, blastn_file) != -1) {
        char query[256];
        char gene[256];
        float identity;
        int align_len;
        int gene_idx = -1;
        
        sscanf(line, "%s\t%s\t%f\t%d", query, gene, &identity, &align_len);
        int ctg_id = rm_contig(query);

        for (i = 0; i < mtpcg_num; i++) {
            if (strcmp(gene, mtpcg[i]) == 0) {
                gene_idx = i;
                break;
            }
        }
        
        if (gene_idx >= 0 && identity > 70 && align_len > 0.4 * mtpcg_len[gene_idx]) {
            k = kh_get(ctg_genes, h_ctg_genes, ctg_id);
            if (k == kh_end(h_ctg_genes)) {
                int ret;
                k = kh_put(ctg_genes, h_ctg_genes, ctg_id, &ret);
                ContigGenes *genes = malloc(sizeof(ContigGenes));
                genes->num_genes = 0;
                genes->gene_ids = NULL;
                genes->identity = NULL;
                genes->align_len = NULL;
                kh_val(h_ctg_genes, k) = genes;
            }
            
            ContigGenes *genes = kh_val(h_ctg_genes, k);
            
            genes->num_genes++;
            if (best_gnum < genes->num_genes) {
                best_gnum = genes->num_genes;
                best_ctg = ctg_id;
            } else if (best_gnum == genes->num_genes && ctg_depth[ctg_id - 1].len > ctg_depth[best_ctg - 1].len) {
                best_ctg = ctg_id;
            }

            genes->gene_ids = realloc(genes->gene_ids, genes->num_genes * sizeof(char*));
            genes->identity = realloc(genes->identity, genes->num_genes * sizeof(float));
            genes->align_len = realloc(genes->align_len, genes->num_genes * sizeof(int));
            
            genes->gene_ids[genes->num_genes - 1] = strdup(gene);
            genes->identity[genes->num_genes - 1] = identity;
            genes->align_len[genes->num_genes - 1] = align_len;
        }
    }

    float hit_depth_filter = 0.0;
    if (best_gnum > 1) {
        hit_depth_filter = ctg_depth[best_ctg - 1].depth;
    } else {
        hit_depth_filter = (10 * filter_depth) / 3;
    }
    float (*calc_score)(float depth, int len, int num_genes);
    switch(taxo) {
        case 0: // Plant
            calc_score = calc_score_plant;
            break;
        case 1: // Animal
            calc_score = calc_score_animal; 
            break;
        case 2: // Fungi
            calc_score = calc_score_fungi;
            break;
        default:
            log_message(ERROR, "Invalid taxo type: %d", taxo);
            kh_destroy(ctg_genes, h_ctg_genes);
            if(line) free(line);
            if(blastn_file) fclose(blastn_file);
            return;
    }

    size_t sort_idx = 0;
    SortPcgCtgs *sort_pcg_ctgs = malloc(sizeof(SortPcgCtgs) * 30);
    
    for (k = kh_begin(h_ctg_genes); k != kh_end(h_ctg_genes); ++k) {
        if (!kh_exist(h_ctg_genes, k)) continue;
        
        int ctg_id = kh_key(h_ctg_genes, k);
        ContigGenes *genes = kh_val(h_ctg_genes, k);
        
        if (genes->num_genes > 0 && ctg_depth[ctg_id - 1].depth > 0.3 * hit_depth_filter) {
            if (sort_idx < 30) {
                sort_pcg_ctgs[sort_idx].ctg = strdup(ctg_depth[ctg_id - 1].ctg);
                sort_pcg_ctgs[sort_idx].score = calc_score(
                    ctg_depth[ctg_id - 1].depth,
                    ctg_depth[ctg_id - 1].len,
                    genes->num_genes
                );
                sort_pcg_ctgs[sort_idx].ctglen = ctg_depth[ctg_id - 1].len;
                sort_pcg_ctgs[sort_idx].ctgdep = ctg_depth[ctg_id - 1].depth;
                sort_idx++;
            }
        }
    }


    if (sort_idx == 0 && verbose == 0) {
        log_message(WARNING, "No seed contigs found (mt), please use GraphBuild command.");
    } else if (sort_idx == 0 && verbose == 1) {
        log_message(WARNING, "No seed contigs found (mt).");
    } else {
        qsort(sort_pcg_ctgs, sort_idx, sizeof(SortPcgCtgs), compare_ctg_scores);
        if (taxo == 0) {
            *ctg_threshold = sort_idx;
        } else {
            *ctg_threshold = 1;
        }
        *candidate_seeds = calloc(*ctg_threshold, sizeof(int));

        log_message(INFO, "Seed finding process is complete.");
        log_info(" _______________________________________________________\n");
        log_info(" Contig Name    Length (bp)   Depth (x)     Score    \n");
        log_info(" -------------  ------------  ------------  ------------\n");
        for (i = 0; i < sort_idx; i++) {
            log_info(" %-12s   %-12d  %-12g  %-10.2f\n", 
                    sort_pcg_ctgs[i].ctg, 
                    sort_pcg_ctgs[i].ctglen, 
                    sort_pcg_ctgs[i].ctgdep, 
                    sort_pcg_ctgs[i].score
                    );
                    if (i < *ctg_threshold) {
                        (*candidate_seeds)[i] = rm_contig(sort_pcg_ctgs[i].ctg);
                    }
                    if (i == 49) {
                        break;
                    }
        }
        log_info(" _______________________________________________________\n");
        log_info("\n");
    }

    for (k = kh_begin(h_ctg_genes); k != kh_end(h_ctg_genes); ++k) {
        if (!kh_exist(h_ctg_genes, k)) continue;
        ContigGenes *genes = kh_val(h_ctg_genes, k);
        for (i = 0; i < genes->num_genes; i++) {
            free(genes->gene_ids[i]);
        }
        free(genes->gene_ids);
        free(genes->identity);
        free(genes->align_len);
        free(genes);
    }
    kh_destroy(ctg_genes, h_ctg_genes);
    
    free(sort_pcg_ctgs); free(line); fclose(blastn_file);
}


void PtHitseeds(const char* exe_path, const char* organelles_type, const char* all_contigs, 
             const char* output_path, int num_threads, int num_ctgs, CtgDepth *ctg_depth, int* candidate_seeds, int ctg_threshold, float filter_depth, int verbose) {

    log_message(INFO, "Finding Pt seeds...");
    char *dir = dirname(strdup(exe_path));
    char* db_path = (char*)malloc(sizeof(*db_path) * (snprintf(NULL, 0, "%s/Conserved_PCGs_db/Plant_conserved_cpgene_nt.fa", dir) + 1));
    sprintf(db_path, "%s/Conserved_PCGs_db/Plant_conserved_cpgene_nt.fa", dir);
    free(dir);

    mkdirfiles(output_path);

    char* blastn_out = (char*)malloc(sizeof(*blastn_out) * (snprintf(NULL, 0, "%s/PMAT_pt_blastn.txt", output_path) + 1));
    sprintf(blastn_out, "%s/PMAT_pt_blastn.txt", output_path);
    
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
    uint64_t i;
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

        // int tem_flag = 0;
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
    

    if (seeds_ctg == 0 && verbose == 0) {
        log_message(WARNING, "No seed contigs found (pt), please use GraphBuild command.");
    } else if (seeds_ctg == 0 && verbose == 1) {
        log_message(WARNING, "No seed contigs found (pt).");
    } else {
        qsort(ptpcg_ctgs, seeds_ctg, sizeof(SortPcgCtgs), compare_ctg_scores);


        log_message(INFO, "Seed finding process is complete.");
        log_info(" _______________________________________________________\n");
        log_info(" Contig Name    Length (bp)   Depth (x)     Score    \n");
        log_info(" -------------  ------------  ------------  ------------\n");
        for (i = 0; i < seeds_ctg; i++) {
            log_info(" %-12s   %-12d  %-12g  %-10.2f\n", 
                    ptpcg_ctgs[i].ctg, 
                    ptpcg_ctgs[i].ctglen, 
                    ptpcg_ctgs[i].ctgdep, 
                    ptpcg_ctgs[i].score
                    );
            if (i < ctg_threshold) {
                candidate_seeds[i] = rm_contig(ptpcg_ctgs[i].ctg);
            }
            if (i == 49) {
                break;
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

    execute_command(command, 0, 0);

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


static int compare_blast_dirc_matches(const void* a, const void* b) {
    const BlastDircMatch* ma = (const BlastDircMatch*)a;
    const BlastDircMatch* mb = (const BlastDircMatch*)b;
    
    int id_a = atoi(ma->query_id + 6);
    int id_b = atoi(mb->query_id + 6);
    if (id_a != id_b) {
        return id_a - id_b;
    }
    
    int cmp = strcmp(ma->gene_id, mb->gene_id);
    if (cmp != 0) return cmp;
    
    return ma->sstart - mb->sstart;
}

static int check_overlap(int start1, int end1, int start2, int end2) {
    return !(end1 < start2 || end2 < start1);
}

static int in_keep_indices(int idx, int* indices, int count) {
    int i;
    for (i = 0; i < count; i++) {
        if (indices[i] == idx) return 1;
    }
    return 0;
}

static void filter_same_gene_matches(BlastDircMatch* matches, int* match_count) {
    if (*match_count <= 1) return;

    int i;
    int* keep_indices = calloc(*match_count, sizeof(int));
    int keep_count = 0;

    char current_ctg[256] = "";
    char current_gene[256] = "";
    int last_keep_idx = -1;

    for (i = 0; i < *match_count; i++) {
        if (strcmp(current_ctg, matches[i].query_id) != 0 || 
           strcmp(current_gene, matches[i].gene_id) != 0) {
            keep_indices[keep_count++] = i;
            last_keep_idx = i;
            strcpy(current_ctg, matches[i].query_id);
            strcpy(current_gene, matches[i].gene_id);
            continue;
        }

        if (check_overlap(matches[last_keep_idx].sstart, matches[last_keep_idx].send,
                        matches[i].sstart, matches[i].send)) {
            if (matches[i].score > matches[last_keep_idx].score || 
               (matches[i].score == matches[last_keep_idx].score && 
                matches[i].align_len > matches[last_keep_idx].align_len)) {
                keep_indices[keep_count-1] = i;
                last_keep_idx = i;
            }
        } else {
            keep_indices[keep_count++] = i;
            last_keep_idx = i;
        }
    }

    BlastDircMatch* temp = malloc(sizeof(BlastDircMatch) * keep_count);
    for (i = 0; i < keep_count; i++) {
        temp[i] = matches[keep_indices[i]];
    }

    for (i = 0; i < *match_count; i++) {
        if(!in_keep_indices(i, keep_indices, keep_count)) {
            free(matches[i].query_id);
            free(matches[i].gene_id);
        }
    }
    memcpy(matches, temp, sizeof(BlastDircMatch) * keep_count);
    *match_count = keep_count;

    free(keep_indices);
    free(temp);
}



void mrun_blastn(const char *all_contigs, const char *db_path, char *final_out, int num_threads, int *num_hits) {
    char* temp_out = (char*)malloc(sizeof(*temp_out) * (snprintf(NULL, 0, "%s.temp", final_out) + 1));
    sprintf(temp_out, "%s.temp", final_out);

    uint64_t i;
    char* command = (char*)malloc(sizeof(*command) * (snprintf(NULL, 0, "blastn -db %s -query %s -outfmt 6 -num_threads %d > %s",
                                db_path, all_contigs, num_threads, temp_out) + 1));
    sprintf(command, "blastn -db %s -query %s -outfmt 6 -num_threads %d > %s",
            db_path, all_contigs, num_threads, temp_out);
    execute_command(command, 0, 0);
    free(command);

    struct stat st;
    stat(temp_out, &st);
    if (st.st_size == 0) {
        FILE *final_fp = fopen(final_out, "w");
        if(final_fp) fclose(final_fp);
        remove(temp_out);
        *num_hits = 0;
        // log_message(WARNING, "No BLAST hits found");
        return;
    }

    FILE *temp_fp = fopen(temp_out, "r");
    if (!temp_fp) {
        log_message(ERROR, "Failed to open temp blast output file");
        remove(temp_out);
        return;
    }

    int total_matches = 0;
    char line[1024];
    while(fgets(line, sizeof(line), temp_fp)) {
        total_matches++;
    }

    rewind(temp_fp);

    BlastDircMatch* matches = (BlastDircMatch*)malloc(sizeof(BlastDircMatch) * total_matches);
    int match_count = 0;

    while (fgets(line, sizeof(line), temp_fp)) {
        if(strlen(line) < 10) continue;

        BlastDircMatch *match = &matches[match_count];
        char gene_full[256];
        
        if(sscanf(line, "%ms\t%[^\t]\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f",
            &match->query_id,
            gene_full,
            &match->identity,
            &match->align_len,
            &match->mismatch,
            &match->gap,
            &match->qstart,
            &match->qend,
            &match->sstart,
            &match->send,
            &match->evalue,
            &match->score) != 12) {
            continue;
        }

        char *last_underscore = strrchr(gene_full, '_');
        if (last_underscore != NULL) {
            match->gene_id = strdup(last_underscore + 1);
        } else {
            match->gene_id = strdup(gene_full);
        }

        if (match->sstart > match->send) {
            int temp = match->sstart;
            match->sstart = match->send;
            match->send = temp;
            match->direction = 2;
        } else {
            match->direction = 1;
        }
        
        match_count++;
    }

    /* Sort the matches by query id, gene id, and start position */
    qsort(matches, match_count, sizeof(BlastDircMatch), compare_blast_dirc_matches);
    
    /* Filter out overlapping matches */
    filter_same_gene_matches(matches, &match_count);

    FILE *dirc_fp = fopen(final_out, "w");
    if (dirc_fp) {
        for (i = 0; i < match_count; i++) {
            BlastDircMatch *m = &matches[i];
            fprintf(dirc_fp, "%s\t%s\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2g\t%.1f\t%d\n",
                m->query_id,
                m->gene_id,
                m->identity,
                m->align_len,
                m->mismatch,
                m->gap,
                m->qstart,
                m->qend,
                m->sstart,
                m->send,
                m->evalue,
                m->score,
                m->direction
            );
        }
        fclose(dirc_fp);
    }

    for (i = 0; i < match_count; i++) {
        free(matches[i].query_id);
        free(matches[i].gene_id);
    }
    free(matches);
    
    fclose(temp_fp);
    remove(temp_out);

    *num_hits = match_count;
}