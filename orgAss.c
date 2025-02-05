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
#include <unistd.h> // access
#include <sys/types.h> // stat
#include <libgen.h> // dirname

#include "misc.h"
#include "khash.h"
#include "orgAss.h"
#include "log.h"
#include "khash.h"
#include "hitseeds.h"


static int get_gene_length(const char* gene, int taxo);
static void print_report_file(FILE* fp, AssessResult* result, int gene_count, CtgDepth* ctg_depth, int num_contigs, const int *contig_ids);
static void print_report_console(AssessResult* result, int gene_count, CtgDepth* ctg_depth, int num_contigs, const int *contig_ids);

typedef struct {
    char* gene_name;
    float identity;
    int align_len;
    int gene_len;
} GeneStats;

typedef struct {
    char** genes;
    GeneStats* stats;
    int count;
} GeneMap;

static void update_gene_stats(GeneMap* map, const char* gene, float identity, int align_len, int gene_len) {
    int i;

    for (i = 0; i < map->count; i++) {
        if (strcmp(map->genes[i], gene) == 0) {
            int old_len = map->stats[i].align_len;
            float old_identity = map->stats[i].identity;
            map->stats[i].align_len += align_len;
            map->stats[i].identity = (old_identity * old_len + identity * align_len) / 
                                   (old_len + align_len);
            return;
        }
    }
    
    map->count++;
    map->genes = realloc(map->genes, map->count * sizeof(char*));
    map->stats = realloc(map->stats, map->count * sizeof(GeneStats));
    map->genes[map->count - 1] = strdup(gene);
    map->stats[map->count - 1].gene_name = strdup(gene);
    map->stats[map->count - 1].identity = identity;
    map->stats[map->count - 1].align_len = align_len;
    map->stats[map->count - 1].gene_len = gene_len;
}

void orgAss(const char* exe_path, const char* all_contigs, CtgDepth *ctg_depth, const char* output_path, const int* contig_ids, int num_contigs, const char* organelle_type, int taxo) {

    AssessResult* result = calloc(1, sizeof(AssessResult));
    
    int i, j;
    khint_t k;

    const char **target_genes;
    int gene_count;
    switch(taxo) {
        case 0:
            target_genes = plt_mtpcg;
            gene_count = plt_mtpcg_num;
            break;
        case 1:
            target_genes = anl_mtpcg;
            gene_count = anl_mtpcg_num;
            break;
        case 2:
            target_genes = fug_mtpcg;
            gene_count = fug_mtpcg_num;
            break;
        default:
            log_message(ERROR, "Invalid taxo type");
            free(result);
            return;
    }


    char blast_file[4096];
    snprintf(blast_file, sizeof(blast_file), "%s/PMAT_mt_blastn.txt", output_path);

    if (access(blast_file, F_OK) == -1) {
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
        char *dir = dirname(strdup(exe_path));
        size_t dir_len = strlen(dir);
        uint64_t i;
        size_t db_path_len = dir_len + strlen(db_suffix) + 1;
        char db_path[db_path_len];
        snprintf(db_path, db_path_len, "%s%s", dir, db_suffix);
        int num_hits = 0;
        mrun_blastn(all_contigs, db_path, blast_file, 6, &num_hits);
        free(dir);
    }

    FILE *fp = fopen(blast_file, "r");
    if (!fp) {
        log_message(ERROR, "Cannot open : %s", blast_file);
        free(result);
        return;
    }

    khash_t(hash_count) *h_counts = kh_init(hash_count);
    result->duplicate_ids = malloc(num_contigs * sizeof(int));
    result->num_duplicates = 0;

    for (i = 0; i < num_contigs; i++) {
        int ret;
        k = kh_get(hash_count, h_counts, contig_ids[i]);
        if (k == kh_end(h_counts)) {
            k = kh_put(hash_count, h_counts, contig_ids[i], &ret);
            kh_value(h_counts, k) = 1;
        } else {
            kh_value(h_counts, k)++;
        }
    }

    for (k = kh_begin(h_counts); k != kh_end(h_counts); ++k) {
        if (kh_exist(h_counts, k)) {
            int id = kh_key(h_counts, k);
            int count = kh_value(h_counts, k);
            if (count > 1) {
                result->duplicate_ids[result->num_duplicates++] = id;
                result->duplicate_contigs++;
            }
        }
    }
    
    result->total_contigs = num_contigs;
    result->contig_stats = calloc(result->total_contigs, sizeof(ContigGeneStats));
    
    GeneMap* gene_maps = calloc(num_contigs, sizeof(GeneMap));
    
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        char query[256], gene[256];
        float identity;
        int align_len, gene_len;
        
        sscanf(line, "%s\t%s\t%f\t%d", query, gene, &identity, &align_len);
        int ctg_id = rm_contig(query);
        
        gene_len = get_gene_length(gene, taxo);
        if (gene_len == 0) continue;
        
        if (identity > 50) {
            for (i = 0; i < num_contigs; i++) {
                if (contig_ids[i] == ctg_id) {
                    ContigGeneStats *stat = &result->contig_stats[i];
                    if (stat->contig_id == 0) {
                        stat->contig_id = ctg_id;
                        stat->num_genes = 0;
                        stat->gene_list = NULL;
                    }
                    update_gene_stats(&gene_maps[i], gene, identity, align_len, gene_len);
                }
            }
        }
    }
    
    char **found_genes = NULL;
    int found_count = 0;
    
    for (i = 0; i < num_contigs; i++) {
        ContigGeneStats *stat = &result->contig_stats[i];
        if(stat->contig_id != 0) {
            int valid_genes = 0;
            stat->gene_list = malloc(gene_maps[i].count * sizeof(char*));
            
            for (j = 0; j < gene_maps[i].count; j++) {
                if (gene_maps[i].stats[j].align_len >= 0.4 * gene_maps[i].stats[j].gene_len) {
                    stat->gene_list[valid_genes] = strdup(gene_maps[i].genes[j]);
                    result->total_genes++;

                    int is_new = 1;
                    int m;
                    for (m = 0; m < found_count; m++) {
                        if(strcmp(found_genes[m], gene_maps[i].genes[j]) == 0) {
                            is_new = 0;
                            break;
                        }
                    }
                    
                    if(is_new) {
                        found_count++;
                        found_genes = realloc(found_genes, found_count * sizeof(char*));
                        found_genes[found_count - 1] = strdup(gene_maps[i].genes[j]);
                    }
                    
                    valid_genes++;
                }
            }
            
            stat->num_genes = valid_genes;
            if (valid_genes < gene_maps[i].count) {
                stat->gene_list = realloc(stat->gene_list, valid_genes * sizeof(char*));
            }
        }
    }
    
    result->unique_genes = found_count;

    for (i = 0; i < found_count; i++) {
        free(found_genes[i]);
    }
    free(found_genes);

    for (i = 0; i < num_contigs; i++) {
        for (j = 0; j < gene_maps[i].count; j++) {
            free(gene_maps[i].genes[j]);
            free(gene_maps[i].stats[j].gene_name);
        }
        free(gene_maps[i].genes);
        free(gene_maps[i].stats);
    }
    free(gene_maps);
    
    fclose(fp);

    char report_file[4096];
    snprintf(report_file, sizeof(report_file), "%s/PMAT_orgAss.txt", output_path);
    FILE *report = fopen(report_file, "w");
    
    print_report_file(report, result, gene_count, ctg_depth, num_contigs, contig_ids);
    fclose(report);

    print_report_console(result, gene_count, ctg_depth, num_contigs, contig_ids);

    if (result->duplicate_ids) {
        free(result->duplicate_ids);
    }

    if(result->contig_stats) {
        for (i = 0; i < result->total_contigs; i++) {
            if (result->contig_stats[i].gene_list) {
                for (j = 0; j < result->contig_stats[i].num_genes; j++) {
                    free(result->contig_stats[i].gene_list[j]);
                }
                free(result->contig_stats[i].gene_list);
            }
        }
        free(result->contig_stats);
    }

    free(result);
    
    kh_destroy(hash_count, h_counts);
    
}


static void print_report_file(FILE* fp, AssessResult* result, int gene_count, CtgDepth *ctg_depth, int num_contigs, const int *contig_ids) {
    int i, j;
    fprintf(fp, "\n");
    fprintf(fp, " ==========================================================\n");
    fprintf(fp, "             Mitochondrial Assembly Assessment             \n");
    fprintf(fp, " ==========================================================\n");
    fprintf(fp, "\n");
    
    int total_length = 0;
    float total_depth = 0;
    int num_gene_contigs = 0;
    for (i = 0; i < num_contigs; i++) {
        total_length += ctg_depth[contig_ids[i] - 1].len;
        total_depth += ctg_depth[contig_ids[i] - 1].depth;
    }
    float avg_depth = total_depth / result->total_contigs;
    
    fprintf(fp, " Basic Statistics:\n");
    fprintf(fp, " ----------------------------------------------------------\n");
    fprintf(fp, " Total contigs:          %-4d\n", result->total_contigs);
    fprintf(fp, " Total length:           %.1f kb\n", total_length/1000.0);
    fprintf(fp, " Average depth:          %.1f x\n", avg_depth);
    fprintf(fp, " Total genes found:      %d/%-2d (%.1f%%)\n", 
             result->unique_genes, gene_count,
             (float)result->unique_genes/gene_count * 100);
    fprintf(fp, " Duplicated contigs:     %-4d\n", result->duplicate_contigs);
    fprintf(fp, "\n");

    fprintf(fp, " Per-contig Details:\n");
    fprintf(fp, " ----------------------------------------------------------\n");
    fprintf(fp, " %-10s  %-8s  %-20s\n", 
           "Contig ID", "Genes", "Gene List");
    fprintf(fp, " ----------------------------------------------------------\n");
    
    int found_genes = 0;
    for (i = 0; i < result->total_contigs; i++) {
        ContigGeneStats* stat = &result->contig_stats[i];
        
        if(stat->num_genes == 0) continue;
        
        found_genes = 1;
        char gene_list[256] = {0};
        
        int offset = 0;
        for (j = 0; j < stat->num_genes && offset < 255; j++) {
            offset += snprintf(gene_list + offset, 255 - offset,
                           "%s%s", stat->gene_list[j],
                           (j == stat->num_genes-1) ? "" : ",");
        }
        
        
        fprintf(fp, " %-10d  %-8d  %-20s\n",
              stat->contig_id, stat->num_genes, gene_list);
    }
    
    if (!found_genes) {
        fprintf(fp, " No contigs contain genes.\n");
    }
    
    fprintf(fp, " ----------------------------------------------------------\n");
    fprintf(fp, "\n");
}


static void print_report_console(AssessResult* result, int gene_count, CtgDepth *ctg_depth, int num_contigs, const int *contig_ids) {
    int i, j;
    log_info("\n");
    log_info(" ==========================================================\n");
    log_info("             Mitochondrial Assembly Assessment             \n");
    log_info(" ==========================================================\n");
    log_info("\n");
    
    int total_length = 0;
    float total_depth = 0;
    int num_gene_contigs = 0;
    for (i = 0; i < num_contigs; i++) {
        total_length += ctg_depth[contig_ids[i] - 1].len;
        total_depth += ctg_depth[contig_ids[i] - 1].depth;
    }
    float avg_depth = total_depth / result->total_contigs;
    
    log_info(" Basic Statistics:\n");
    log_info(" ----------------------------------------------------------\n");
    log_info(" Total contigs:          %-4d\n", result->total_contigs);
    log_info(" Total length:           %.1f kb\n", total_length/1000.0);
    log_info(" Average depth:          %.1f x\n", avg_depth);
    log_info(" Total genes found:      %d/%-2d (%.1f%%)\n", 
             result->unique_genes, gene_count,
             (float)result->unique_genes/gene_count * 100);
    log_info(" Duplicated contigs:     %-4d\n", result->duplicate_contigs);
    log_info("\n");

    log_info(" Per-contig Details:\n");
    log_info(" ----------------------------------------------------------\n");
    log_info(" %-10s  %-8s  %-20s\n", 
             "Contig ID", "Genes", "Gene List");
    log_info(" ----------------------------------------------------------\n");
    
    int found_genes = 0;
    for (i = 0; i < result->total_contigs; i++) {
        ContigGeneStats* stat = &result->contig_stats[i];
        
        if(stat->num_genes == 0) continue;
        
        found_genes = 1;
        char gene_list[256] = {0};
        
        int offset = 0;
        for (j = 0; j < stat->num_genes && offset < 255; j++) {
            offset += snprintf(gene_list + offset, 255 - offset,
                           "%s%s", stat->gene_list[j],
                           (j == stat->num_genes-1) ? "" : ",");
        }
        
        if(strlen(gene_list) > 20) {
            gene_list[17] = '.';
            gene_list[18] = '.';
            gene_list[19] = '.';
            gene_list[20] = '\0';
        }
        
        log_info(" %-10d  %-8d  %-20s\n",
                stat->contig_id, stat->num_genes, gene_list);
    }
    
    if (!found_genes) {
        log_info(" No contigs contain genes.\n");
    }
    
    log_info(" ----------------------------------------------------------\n");
    log_info("\n");
}

static int get_gene_length(const char* gene, int taxo) {
    const char **gene_list;
    const int *len_list;
    int count;
    int i;
    switch(taxo) {
        case 0:
            gene_list = plt_mtpcg;
            len_list = plt_mtpcg_len;
            count = plt_mtpcg_num;
            break;
        case 1:
            gene_list = anl_mtpcg;
            len_list = anl_mtpcg_len;
            count = anl_mtpcg_num;
            break;
        case 2:
            gene_list = fug_mtpcg;
            len_list = fug_mtpcg_len;
            count = fug_mtpcg_num;
            break;
        default:
            return 0;
    }
    
    for (i = 0; i < count; i++) {
        if(strcmp(gene, gene_list[i]) == 0) {
            return len_list[i];
        }
    }
    return 0;
}

