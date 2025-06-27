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
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> // dirname
#include <unistd.h>

#include "khash.h"
#include "graphtools.h"
#include "misc.h"
#include "log.h"
#include "BFSseed.h"
#include "hitseeds.h"
#include "orgAss.h"

typedef struct {
    int* node;
    int num;
} nodeNum;

typedef struct {
    char seq[1001];
    int type;
} nodeKmer;

KHASH_MAP_INIT_INT(Ha_nodelink, nodeNum)
KHASH_MAP_INIT_INT(Ha_nodekmer, int)
KHASH_MAP_INIT_INT(Ha_direct, int)



void addseq(const char* allgraph, const char* all_fna, CtgDepth* ctgdepth) {
    
    char* line = NULL;
    size_t len = 0;

    FILE* temp_fpall = fopen(all_fna, "r");
    int fnactg_num = 0;
    while (getline(&line, &len, temp_fpall) != -1) {
        if (line[0] == '>') {
            fnactg_num++;
        }
    }
    rewind(temp_fpall);
    int* fnactg = malloc(fnactg_num * sizeof(int));
    int fnactg_idx = 0;
    while (getline(&line, &len, temp_fpall) != -1) {
        if (line[0] == '>') {
            char *token = strtok(line, " ");
            char *tempctg = token + 1;
            int intctg = rm_contig(tempctg);
            fnactg[fnactg_idx] = intctg;
            fnactg_idx++;
        }
    }
    fclose(temp_fpall);

    FILE* fpgraph = fopen(allgraph, "r");
    FILE* fpout = fopen(all_fna, "a");
    int flag_I = 0;
    while (getline(&line, &len, fpgraph) != -1) {
        if (line[0] == 'I') {
            flag_I = 1;
            char *token = strtok(line, "\t");
            token = strtok(NULL, "\t");
            int tempctg = atoi(token);
            token = strtok(NULL, "\t");
            char *tempseq = strdup(token);
            if (findint(fnactg, fnactg_num, tempctg) == 0) {
                fprintf(fpout, ">%s length=%d numreads=%g\n", ctgdepth[tempctg - 1].ctg, ctgdepth[tempctg - 1].len, ctgdepth[tempctg - 1].depth);
                char* seq = strdup(tempseq);
                to_upper(seq);
                size_t seq_len = strlen(seq);
                size_t pos = 0;
                while (pos < seq_len) {
                    size_t chunk_size = (seq_len - pos >= 60) ? 60 : seq_len - pos;
                    fprintf(fpout, "%.*s\n", (int)chunk_size, seq + pos);
                    pos += chunk_size;
                }
                free(seq);
            }
        } else if (flag_I) {
            break;
        }
    }

    /* free memory */
    free(line);
    free(fnactg);
    fclose(fpout);
    fclose(fpgraph);
}


static void maingraph(BFSlinks* bfslinks, BFSlinks* mainlinks, CtgDepth* ctg_depth, int num_dynseeds, int num_bfslinks, int* mainseeds, 
                     int* main_num, int* mainseeds_num, int* rm_ctg, int rm_num, float filter_depth) {

    uint64_t i;
    BFSlinks* templinks = malloc(num_bfslinks * sizeof(BFSlinks));
    for (i = 0; i < num_bfslinks; i++) {
        copy_BFSlinks(&templinks[i], &bfslinks[i]);
    }
    *main_num = num_bfslinks;
    int flag_num = 0;

    while(1) {
        *mainseeds_num = 0;
        khash_t(Ha_direct) *map_direct = kh_init(Ha_direct);
        for (i = 0; i < *main_num; i++) {
            khint_t k;
            int ret;
            int temp_lctgsmp = templinks[i].lctgsmp;
            int temp_rctgsmp = templinks[i].rctgsmp;
            if (temp_lctgsmp != temp_rctgsmp) {
                k = kh_get(Ha_direct, map_direct, temp_lctgsmp);
                int lutr_value = templinks[i].lutrsmp;
                if (k != kh_end(map_direct) && kh_val(map_direct, k) != lutr_value) {
                    mainseeds[*mainseeds_num] = temp_lctgsmp;
                    (*mainseeds_num)++;
                    kh_val(map_direct, k) = 0;
                } else if (k == kh_end(map_direct)) {
                    mainseeds[*mainseeds_num] = temp_lctgsmp;
                    (*mainseeds_num)++;
                    k = kh_put(Ha_direct, map_direct, temp_lctgsmp, &ret);
                    kh_val(map_direct, k) = templinks[i].lutrsmp;              
                }

                k = kh_get(Ha_direct, map_direct, temp_rctgsmp);
                int rutr_value = templinks[i].rutrsmp;
                if (k != kh_end(map_direct) && kh_val(map_direct, k) != rutr_value) {
                    mainseeds[*mainseeds_num] = temp_rctgsmp;
                    (*mainseeds_num)++;
                    kh_val(map_direct, k) = 0;
                } else if (k == kh_end(map_direct)) {
                    mainseeds[*mainseeds_num] = temp_rctgsmp;
                    (*mainseeds_num)++;
                    k = kh_put(Ha_direct, map_direct, temp_rctgsmp, &ret);
                    kh_val(map_direct, k) = templinks[i].rutrsmp;                    
                }

            } else {
                mainseeds[*mainseeds_num] = temp_lctgsmp;
                (*mainseeds_num)++;
                mainseeds[*mainseeds_num] = temp_rctgsmp;
                (*mainseeds_num)++;
                k = kh_get(Ha_direct, map_direct, temp_lctgsmp);
                if (k == kh_end(map_direct)) {
                    k = kh_put(Ha_direct, map_direct, temp_lctgsmp, &ret);
                    kh_val(map_direct, k) = 0;
                }
                // k = kh_get(Ha_direct, map_direct, temp_rctgsmp);
                // if (k == kh_end(map_direct)) {
                //     mainseeds[*mainseeds_num] = temp_rctgsmp;
                //     (*mainseeds_num)++;
                //     k = kh_put(Ha_direct, map_direct, temp_lctgsmp, &ret);
                //     kh_val(map_direct, k) = 0;
                // }
            }       

        }

        kh_destroy(Ha_direct, map_direct);

        if (*mainseeds_num == 0) {
            break;
        }

        removeUnique(mainseeds, mainseeds_num);

        flag_num = 0;
        for (i = 0; i < *main_num; i++) {
            if (findint(mainseeds, *mainseeds_num, templinks[i].lctgsmp) == 1 &&
                findint(mainseeds, *mainseeds_num, templinks[i].rctgsmp) == 1 &&
                findint(rm_ctg, rm_num, templinks[i].lctgsmp) == 0 &&
                findint(rm_ctg, rm_num, templinks[i].rctgsmp) == 0 &&
                ctg_depth[templinks[i].lctgsmp - 1].depth >= filter_depth &&
                ctg_depth[templinks[i].rctgsmp - 1].depth >= filter_depth) {

                // mainlinks[flag_num] = templinks[i];
                copy_BFSlinks(&mainlinks[flag_num], &templinks[i]);
                flag_num++;
            }
        }

        if (flag_num == *main_num) {
            break;
        }
        *main_num = flag_num;
        for (i = 0; i < *main_num; i++) {
            copy_BFSlinks(&templinks[i], &mainlinks[i]);
        }
        
    }
    // free
    free(templinks);
    removeDup(mainseeds, mainseeds_num);
}

static void run_db(const char *cutseq) {

    log_message(INFO, "Building database...");
    size_t dbcomond_len = snprintf(NULL, 0, "makeblastdb -in %s -dbtype nucl -out %s.db", 
                                    cutseq, cutseq) + 1;
    char* dbcommand = malloc(dbcomond_len);
    snprintf(dbcommand, dbcomond_len, "makeblastdb -in %s -dbtype nucl -out %s.db", 
                                    cutseq, cutseq);

    execute_command(dbcommand, 0, 0);
}

static void run_blastn(const char *cutseq, const char *ctgseq, char *blastn_out, int num_threads, int *num_hits) {

    size_t comond_len = snprintf(NULL, 0, "blastn -query %s -db %s.db -outfmt 6 -num_threads %d > %s", 
                                    ctgseq, cutseq, num_threads, blastn_out) + 1;
    char* command = malloc(comond_len);
    snprintf(command, comond_len, "blastn -query %s -db %s.db -outfmt 6 -num_threads %d > %s", 
                                    ctgseq, cutseq, num_threads, blastn_out);

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


void optgfa(const char* exe_path, int num_dynseeds, int** dynseeds, BFSlinks** bfslinks, int* num_bfslinks, 
            CtgDepth* ctgdepth, const char* output, const char* all_fna, const char* allgraph, 
            const char* organelles_type, int* mainseeds_num, int** mainseeds, int interfering_ctg_num, 
            int* interfering_ctg, int taxo, float filter_depth, char* cutseq) 
{
    // addseq(allgraph, all_fna, ctgdepth);

    FILE* fpall = fopen(all_fna, "r");
    if (fpall == NULL) {
        log_message(ERROR, "Failed to open %s", all_fna);
        exit(EXIT_FAILURE);
    }

    fnainfo* fnainfos = (fnainfo*) calloc(num_dynseeds, sizeof(fnainfo));
    
    char* line = NULL;
    size_t len = 0;
    int tmp_num_dynseeds = 0;
    int flag_seq = 0;
    char* tempseq = NULL;
    size_t tempseq_len = 0;
    int flag_end = 0;
    int* temp_dynseed = calloc(num_dynseeds, sizeof(int));
    khash_t(Ha_nodekmer) *nodeKmer_hash = kh_init(Ha_nodekmer);
    size_t kmer1000_len = snprintf(NULL, 0, "%s/Kmer1000.fa", output) + 1;
    char kmer1000[kmer1000_len];
    snprintf(kmer1000, kmer1000_len, "%s/Kmer1000.fa", output);
    FILE* kout = fopen(kmer1000, "w");
    while (getline(&line, &len, fpall) != -1) 
    {
        if (line[0] == '>') 
        {
            char *token = strtok(line, " ");
            char *tempctg = token + 1;

            int intctg = rm_contig(tempctg);
            if (findint(*dynseeds, num_dynseeds, intctg) == 1) {
                fnainfos[tmp_num_dynseeds].ctg = intctg;
                
                if (tmp_num_dynseeds > 0) {
                    fnainfos[tmp_num_dynseeds - 1].seq = strdup(tempseq);
                    to_upper(fnainfos[tmp_num_dynseeds - 1].seq);

                    if (ctgdepth[fnainfos[tmp_num_dynseeds - 1].ctg - 1].len > 1000) {
                        char seqKmer[1001];

                        char* tempseq_end = tempseq + strlen(tempseq) - 500;
                        char* tempseq_start = tempseq;
                        strncpy(seqKmer, tempseq_end, 500);
                        seqKmer[500] = '\0';
                        strncat(seqKmer, tempseq_start, 500);
                        int ret;
                        khint_t k = kh_put(Ha_nodekmer, nodeKmer_hash, fnainfos[tmp_num_dynseeds - 1].ctg, &ret);
                        kh_value(nodeKmer_hash, k) = 1;
                        fprintf(kout, ">%d\n%s\n", fnainfos[tmp_num_dynseeds - 1].ctg, seqKmer);
                    }

                    free(tempseq);
                    flag_end = 0;
                }
                tempseq = malloc(1);
                tempseq[0] = '\0';
                flag_seq = 1;
                if (findint(temp_dynseed, num_dynseeds, intctg) == 1) {
                    log_message(ERROR, "%s is not a valid contig in %s", ctgdepth[intctg - 1].ctg, all_fna);
                    exit(EXIT_FAILURE);
                }
                temp_dynseed[tmp_num_dynseeds] = intctg;
                tmp_num_dynseeds++;
                tempseq_len = 0;
            } else {
                flag_seq = 0;
            }

        } else if (flag_seq) 
        {
            flag_end = 1;
            line[strcspn(line, "\n")] = 0;
            size_t seqlen = strlen(line);

            tempseq = realloc(tempseq, tempseq_len + seqlen + 1);
            if (tempseq == NULL) {
                log_message(ERROR, "Failed to allocate memory for tempseq");
                exit(EXIT_FAILURE);
            }
            strcat(tempseq, line);
            tempseq_len += seqlen;
        }
    }

    if (flag_end) {
        fnainfos[tmp_num_dynseeds - 1].seq = strdup(tempseq);
        to_upper(fnainfos[tmp_num_dynseeds - 1].seq);

        if (ctgdepth[fnainfos[tmp_num_dynseeds - 1].ctg - 1].len > 1000) {
            char seqKmer[1001];

            char* tempseq_end = tempseq + strlen(tempseq) - 500;
            char* tempseq_start = tempseq;
            strncpy(seqKmer, tempseq_end, 500);
            seqKmer[500] = '\0';
            strncat(seqKmer, tempseq_start, 500);
            int ret;
            khint_t k = kh_put(Ha_nodekmer, nodeKmer_hash, fnainfos[tmp_num_dynseeds - 1].ctg, &ret);
            kh_value(nodeKmer_hash, k) = 1;
            fprintf(kout, ">%d\n%s\n", fnainfos[tmp_num_dynseeds - 1].ctg, seqKmer);
        }

        free(tempseq);
    }
    free(line);
    fclose(kout);
    fclose(fpall);

    /* circular / linear */
    int num_hits = 0;
    // char* cutseq = malloc(strlen(output) + strlen("/subsample/PMAT_cut_seq.fa") + 1);
    // if (cutseq == NULL) {
    //     fprintf(stderr, "Memory allocation failed!\n");
    //     exit(EXIT_FAILURE);
    // }
    // snprintf(cutseq, strlen(output) + strlen("/subsample/PMAT_cut_seq.fa") + 1, 
    //     "%s/subsample/PMAT_cut_seq.fa", output);
    char* cutdb = malloc(sizeof(*cutdb) * (snprintf(NULL, 0, "%s.db.ndb", cutseq) + 1));
    snprintf(cutdb, sizeof(*cutdb) * (snprintf(NULL, 0, "%s.db.ndb", cutseq) + 1), 
        "%s.db.ndb", cutseq);

    char* blast_out = (char*) malloc(sizeof(*blast_out) * (snprintf(NULL, 0, "%s/PMAT_kmer1000.txt", output) + 1));
    snprintf(blast_out, sizeof(*blast_out) * (snprintf(NULL, 0, "%s/PMAT_kmer1000.txt", output) + 1), 
        "%s/PMAT_kmer1000.txt", output);
    
    if (is_file(cutdb) == 0) run_db(cutseq);
    run_blastn(cutseq, kmer1000, blast_out, 8, &num_hits);
    remove_file(kmer1000);
    /* the type of contig structure */
    uint64_t i, j, n;
    if (num_hits > 0) {
        FILE* fpblast = fopen(blast_out, "r");
        if (fpblast == NULL) {
            log_message(ERROR, "Failed to open blast_out file %s", blast_out);
            exit(EXIT_FAILURE);
        }
        
        size_t blt_len;
        char *blt_line = NULL;
        while ((blt_len = getline(&blt_line, &len, fpblast)) != -1) {
            int tempctg;
            float tempiden;
            int temps;
            int tempe;

            if (blt_line[0] == '#') {
                continue;
            } else {
                sscanf(blt_line, "%d\t%*d\t%f\t%*d\t%*d\t%*d\t%d\t%d\t%*f\t%*d", &tempctg, &tempiden, &temps, &tempe);
            }
            if (kh_value(nodeKmer_hash, kh_get(Ha_nodekmer, nodeKmer_hash, tempctg)) != 0){
                if (tempiden > 0.99 && temps <= 450 && tempe >= 550) {
                    int add_flag = 1;
                    for (i = 0; i < *num_bfslinks; i++) {
                        if ((*bfslinks)[i].lctgsmp == tempctg && (*bfslinks)[i].rctgsmp == tempctg) {
                            add_flag = 0;
                            kh_value(nodeKmer_hash, kh_get(Ha_nodekmer, nodeKmer_hash, tempctg)) = 0;
                            break;
                        }
                    }
                    kh_value(nodeKmer_hash, kh_get(Ha_nodekmer, nodeKmer_hash, tempctg)) += 1;

                    if (add_flag && kh_value(nodeKmer_hash, kh_get(Ha_nodekmer, nodeKmer_hash, tempctg)) >= 4) {
                        (*bfslinks) = realloc(*bfslinks, (*num_bfslinks + 1)*sizeof(BFSlinks));
                        if (*bfslinks == NULL) {
                            log_message(ERROR, "Failed to allocate memory for BFSlinks");
                            exit(EXIT_FAILURE);
                        }
                        // (*bfslinks)[*num_bfslinks].lctg = strdup(ctgdepth[tempctg - 1].ctg);
                        (*bfslinks)[*num_bfslinks].lutrsmp = 5;
                        (*bfslinks)[*num_bfslinks].lctgsmp = tempctg;
                        (*bfslinks)[*num_bfslinks].lctglen = ctgdepth[tempctg - 1].len;
                        
                        // (*bfslinks)[*num_bfslinks].rctg = strdup(ctgdepth[tempctg - 1].ctg);
                        (*bfslinks)[*num_bfslinks].rutrsmp = 3;
                        (*bfslinks)[*num_bfslinks].rctgsmp = tempctg;
                        (*bfslinks)[*num_bfslinks].rctglen = ctgdepth[tempctg - 1].len;
                        
                        (*bfslinks)[*num_bfslinks].linkdepth = ctgdepth[tempctg - 1].depth;
                        (*num_bfslinks)++;
                        
                        kh_value(nodeKmer_hash, kh_get(Ha_nodekmer, nodeKmer_hash, tempctg)) = 0;
                    }
                }
            }
        }
    }
    remove_file(blast_out);
    free(cutdb); free(blast_out);

    // create output/gfa_result directory
    char* gfa_output = (char*) malloc(sizeof(*gfa_output) * (snprintf(NULL, 0, "%s/gfa_result", output) + 1));
    sprintf(gfa_output, "%s/gfa_result", output);
    mkdirfiles(gfa_output);

    size_t maingfa_out_len = snprintf(NULL, 0, "%s/PMAT_%s_main.gfa", gfa_output, organelles_type) + 1;
    char maingfa[maingfa_out_len];
    snprintf(maingfa, maingfa_out_len, "%s/PMAT_%s_main.gfa", gfa_output, organelles_type);

    size_t rawgfa_out_len = snprintf(NULL, 0, "%s/PMAT_%s_raw.gfa", gfa_output, organelles_type) + 1;
    char rawgfa[rawgfa_out_len];
    snprintf(rawgfa, rawgfa_out_len, "%s/PMAT_%s_raw.gfa", gfa_output, organelles_type);

    size_t rawfa_len = snprintf(NULL, 0, "%s/gfa_%s.fa", gfa_output, organelles_type) + 1;
    char rawfa[rawfa_len];
    snprintf(rawfa, rawfa_len, "%s/gfa_%s.fa", gfa_output, organelles_type);


    /* raw graph */
    FILE* fprawgfa = fopen(rawgfa, "w");
    FILE* fprawfa = fopen(rawfa, "w");
    if (fprawgfa == NULL) {
        log_message(ERROR, "Failed to open rawgfa file %s", rawgfa);
        exit(EXIT_FAILURE);
    }
    if (fprawfa == NULL) {
        log_message(ERROR, "Failed to open rawfa file %s", rawfa);
        exit(EXIT_FAILURE);
    }
    log_message(INFO, "Raw seeds (%s): %d", organelles_type, num_dynseeds);
    if (num_dynseeds > 0) {
        for (i = 0; i < num_dynseeds; i++) {
            int seed = (*dynseeds)[i];
            int ctg_RC = ctgdepth[seed - 1].len * ctgdepth[seed - 1].depth;

            for (j = 0; j < num_dynseeds; j++) {
                if (seed == fnainfos[j].ctg) {
                    fprintf(fprawgfa, "S\t%d\t%s\tLN:i:%d\tRC:i:%d\n", seed, fnainfos[j].seq, ctgdepth[seed - 1].len, ctg_RC);
                    fprintf(fprawfa, ">%d Len:%d Dep:%.2f\n%s\n", seed, ctgdepth[seed - 1].len, ctgdepth[seed - 1].depth, fnainfos[j].seq);
                }
            }
        }

        for (i = 0; i < *num_bfslinks; i++) {
            fprintf(fprawgfa, "L\t%d\t", (*bfslinks)[i].lctgsmp);
            if ((*bfslinks)[i].lutrsmp == (*bfslinks)[i].rutrsmp) {
                if ((*bfslinks)[i].lutrsmp == 3) {
                    fprintf(fprawgfa, "-\t%d\t+\t0M\n", (*bfslinks)[i].rctgsmp);
                } else if ((*bfslinks)[i].lutrsmp == 5) {
                    fprintf(fprawgfa, "+\t%d\t-\t0M\n", (*bfslinks)[i].rctgsmp);
                }
            } else {
                if ((*bfslinks)[i].lutrsmp == 3) {
                    fprintf(fprawgfa, "-\t%d\t-\t0M\n", (*bfslinks)[i].rctgsmp);
                } else if ((*bfslinks)[i].lutrsmp == 5) {
                    fprintf(fprawgfa, "+\t%d\t+\t0M\n", (*bfslinks)[i].rctgsmp);
                }
            }
        }
    } else {
        log_message(WARNING, "No raw seeds found.");
    }
    fclose(fprawgfa);
    fclose(fprawfa);



    BFSlinks* mainlinks = (BFSlinks*) malloc((*num_bfslinks) * sizeof(BFSlinks));
    int main_num = 0;
    int linear_f = 0;
    
    // // int* mainseeds = (int*) malloc((*num_bfslinks) * 2 * sizeof(int));
    // maingraph(*bfslinks, mainlinks, ctgdepth, num_dynseeds, *num_bfslinks, *mainseeds, &main_num, mainseeds_num, NULL, 0, 0);
    // if (*mainseeds_num == 0) {
    //     linear_f = 1;
    //     for (int i = 0; i < *num_bfslinks; i++) {
    //         copy_BFSlinks(&mainlinks[i], &(*bfslinks)[i]);
    //         main_num++;
    //     }
    //     for (int i = 0; i < num_dynseeds; i++) {
    //         (*mainseeds)[i] = (*dynseeds)[i];
    //         (*mainseeds_num)++;
    //     }
    // } else if (strcmp(organelles_type, "pt") == 0 || taxo == 1) {
    //     if (findint(*mainseeds, *mainseeds_num, (*dynseeds)[0]) == 0) {
    //         linear_f = 1;
    //         (*mainseeds_num) = num_dynseeds;
    //         main_num = *num_bfslinks;
    //         for (int i = 0; i < *num_bfslinks; i++) {
    //             copy_BFSlinks(&mainlinks[i], &(*bfslinks)[i]);
    //         }
    //         for (int i = 0; i < num_dynseeds; i++) {
    //             (*mainseeds)[i] = (*dynseeds)[i];
    //         }
    //     }
    // }
    // if (*mainseeds_num > 0) {

    //     /* main graph */
    //     FILE* fpmaingfa = fopen(maingfa, "w");
    //     if (fpmaingfa == NULL) {
    //         log_message(ERROR, "Failed to open maingfa file %s", maingfa);
    //         exit(EXIT_FAILURE);
    //     }
    //     for (int i = 0; i < *mainseeds_num; i++) {
    //         int seed = (*mainseeds)[i];
    //         int ctg_RC = ctgdepth[seed - 1].len * ctgdepth[seed - 1].depth;

    //         for (int j = 0; j < num_dynseeds; j++) {
    //             if (seed == fnainfos[j].ctg) {
    //                 fprintf(fpmaingfa, "S\t%d\t%s\tLN:i:%d\tRC:i:%d\n", seed, fnainfos[j].seq, ctgdepth[seed - 1].len, ctg_RC);
    //             }
    //         }
    //     }
    //     for (int i = 0; i < main_num; i++) {
    //         fprintf(fpmaingfa, "L\t%d\t", mainlinks[i].lctgsmp);
    //         if (mainlinks[i].lutrsmp == mainlinks[i].rutrsmp) {
    //             if (mainlinks[i].lutrsmp == 3) {
    //                 fprintf(fpmaingfa, "-\t%d\t+\t0M\n", mainlinks[i].rctgsmp);
    //             } else if (mainlinks[i].lutrsmp == 5) {
    //                 fprintf(fpmaingfa, "+\t%d\t-\t0M\n", mainlinks[i].rctgsmp);
    //             }
    //         } else {
    //             if (mainlinks[i].lutrsmp == 3) {
    //                 fprintf(fpmaingfa, "-\t%d\t-\t0M\n", mainlinks[i].rctgsmp);
    //             } else if (mainlinks[i].lutrsmp == 5) {
    //                 fprintf(fpmaingfa, "+\t%d\t+\t0M\n", mainlinks[i].rctgsmp);
    //             }
    //         }
    //     }
    //     fclose(fpmaingfa);
    // }
    int ass_ctg_num = 0;
    int ass_ctg_mall = 100;
    int* ass_ctg_arr = (int*) malloc(ass_ctg_mall * sizeof(int));
    int ps_num = 0;
    pathScore *ps_struct = (pathScore*) malloc(100 * sizeof(pathScore));
    khint_t k, ka, kb;
    int rm_flag = 0;
    if (1) {
        /* Capturing all mitochondrial structures */
        khash_t(Ha_structures)* h_structures = kh_init(Ha_structures);
        uint32_t structure_num = bfs_structure(num_dynseeds, *num_bfslinks, *bfslinks, *dynseeds, h_structures);
        int struc = 0;
        int main_seeds = 0;
        uint64_t max_structure_num = 1;
        
        float temp_filter_depth = 0.0;
        float maxscore = 0.0;
        
        for (i = 0; i < structure_num; i++) 
        {
            k = kh_get(Ha_structures, h_structures, i);
            if (k != kh_end(h_structures)) {
                BFSstructure* temp_struct = kh_value(h_structures, k);
                for (j = 0; j < temp_struct->num_nodes; j++) {
                    if (maxscore < ctgdepth[(temp_struct->node[j]) - 1].len &&
                        findint(interfering_ctg, interfering_ctg_num, temp_struct->node[j]) == 0) {
                        temp_filter_depth = ctgdepth[(temp_struct->node[j]) - 1].depth;
                        maxscore = ctgdepth[(temp_struct->node[j]) - 1].len;
                    }
                }
            }
        }
        if (temp_filter_depth < filter_depth) {
            temp_filter_depth = (10/3) * filter_depth;
        }

        for (i = 0; i < structure_num; i++)
        {
            max_structure_num = 1;
            linear_f = 0;
            struc++;
            k = kh_get(Ha_structures, h_structures, i);
            if (k != kh_end(h_structures)) {
                BFSstructure* temp_struct = kh_value(h_structures, k);

                /* links > 0 and nodes > 0 */
                BFSlinks* temp_mainlinks = (BFSlinks*) malloc((*num_bfslinks) * sizeof(BFSlinks));
                int* temp_mainseeds = (int*) malloc(2 * temp_struct->num_links * sizeof(int));
                int temp_main_num = 0;
                int temp_mainseeds_num = 0;

                BFSlinks* temp_structlinks = (BFSlinks*) malloc((temp_struct->num_links) * sizeof(BFSlinks));
                for (j = 0; j < temp_struct->num_links; j++) {
                    copy_BFSlinks(&temp_structlinks[j], &temp_struct->links[j]);
                }
                maingraph(temp_structlinks, temp_mainlinks, ctgdepth, temp_struct->num_nodes, temp_struct->num_links, temp_mainseeds, 
                            &temp_main_num, &temp_mainseeds_num, NULL, 0, 0.3*temp_filter_depth);
                if (temp_mainseeds_num == 0) {
                    linear_f = 1;
                    for (j = 0; j < temp_struct->num_nodes; j++) {
                        if (ctgdepth[(temp_struct->node[j]) - 1].depth > 0.3*temp_filter_depth) {
                            temp_mainseeds_num++;
                        }
                    }

                    if (temp_mainseeds_num > 0) {
                        temp_mainseeds_num = 0;
                        for (j = 0; j < temp_struct->num_nodes; j++) {
                            temp_mainseeds[j] = temp_struct->node[j];
                            temp_mainseeds_num++;
                        }
                        for (j = 0; j < temp_struct->num_links; j++) {
                            copy_BFSlinks(&temp_mainlinks[j], &temp_struct->links[j]);
                            temp_main_num++;
                        }
                    }

                } else if (strcmp(organelles_type, "pt") == 0 || taxo == 1 || taxo == 2) {
                    float flag_depth = ctgdepth[(*dynseeds)[0] - 1].depth;
                    for (j = 0; j < temp_struct->num_nodes; j++) {
                        float temp_depth = ctgdepth[(temp_struct->node[j]) - 1].depth;
                        int temp_len = ctgdepth[(temp_struct->node[j]) - 1].len;
                        if (temp_depth > 0.4 * flag_depth &&
                            temp_depth < 2 * flag_depth &&
                            temp_len > 100 &&
                            findint(temp_mainseeds, temp_mainseeds_num, temp_struct->node[j]) == 0 &&
                            findint(temp_struct->node, temp_struct->num_nodes, temp_struct->node[j]) == 1) {

                            linear_f = 1;
                            temp_mainseeds_num = temp_struct->num_nodes;
                            temp_main_num = temp_struct->num_links;
                            for (n = 0; n < temp_struct->num_links; n++) {
                                copy_BFSlinks(&temp_mainlinks[n], &temp_struct->links[n]);
                            }
                            for (n = 0; n < temp_struct->num_nodes; n++) {
                                temp_mainseeds[n] = temp_struct->node[n];
                            }
                            break;
                        }
                    }
                    // if (findint(temp_mainseeds, temp_mainseeds_num, (*dynseeds)[0]) == 0 &&
                    //     findint(temp_struct->node, temp_struct->num_nodes, (*dynseeds)[0]) == 1) {
                    //     linear_f = 1;
                    //     temp_mainseeds_num = temp_struct->num_nodes;
                    //     temp_main_num = temp_struct->num_links;
                    //     for (int j = 0; j < temp_struct->num_links; j++) {
                    //         copy_BFSlinks(&temp_mainlinks[j], &temp_struct->links[j]);
                    //     }
                    //     for (int j = 0; j < temp_struct->num_nodes; j++) {
                    //         temp_mainseeds[j] = temp_struct->node[j];
                    //     }
                    // }
                }
                
                /* skip complex structures (node > 500) */
                if (temp_mainseeds_num > 200) {
                    struc--;
                    continue;
                }

                for (j = 0; j < temp_mainseeds_num; j++) 
                {
                    if (findint(*mainseeds, *mainseeds_num, temp_mainseeds[j]) == 0) {
                        (*mainseeds)[(*mainseeds_num)] = temp_mainseeds[j];
                        (*mainseeds_num)++;
                    }
                }


                /* h_links 1:{2,3,4;3}；2:{1,5,6;3}；3:{1,4;2}... */
                khash_t(Ha_nodelink)* h_links = kh_init(Ha_nodelink);

                for (j = 0; j < temp_mainseeds_num; j++) 
                {
                    nodeNum tempNodeNum;
                    tempNodeNum.num = 0;
                    tempNodeNum.node = (int*) malloc((temp_mainseeds_num + 1) * sizeof(int));
                    int ret;

                    ka = kh_put(Ha_nodelink, h_links, temp_mainseeds[j], &ret);
                    for (n = 0; n < temp_main_num; n++) 
                    {
                        if (temp_mainlinks[n].lctgsmp == temp_mainseeds[j]) {
                            tempNodeNum.node[tempNodeNum.num] = temp_mainlinks[n].rctgsmp;
                            tempNodeNum.num++;
                        } else if (temp_mainlinks[n].rctgsmp == temp_mainseeds[j])
                        {
                            tempNodeNum.node[tempNodeNum.num] = temp_mainlinks[n].lctgsmp;
                            tempNodeNum.num++;
                        }
                        
                    }

                    kh_value(h_links, ka) = tempNodeNum;
                    if (tempNodeNum.num > 1) max_structure_num = max_structure_num * (tempNodeNum.num - 1);
                }

                int transfer_num = 0;
                int non_transfer_num = 0;
                
                if (interfering_ctg_num > 0) {
                    int* transfer_ctg = (int*) malloc(interfering_ctg_num * sizeof(int));
                    for (ka = kh_begin(h_links); ka != kh_end(h_links); ++ka) 
                    {
                        if (kh_exist(h_links, ka)) {
                            int key = kh_key(h_links, ka);
                            
                            int flag_transfer = 0;

                            nodeNum tempNodeNum = kh_value(h_links, kh_get(Ha_nodelink, h_links, key));

                            if (findint(interfering_ctg, interfering_ctg_num, key) == 0 && ctgdepth[key - 1].len > 30) {
                                for (j = 0; j < tempNodeNum.num; j++) 
                                {
                                    if (findint(interfering_ctg, interfering_ctg_num, tempNodeNum.node[j]) == 0)
                                        flag_transfer++;
                                }
                            }

                            for (j = 0; j < tempNodeNum.num; j++) 
                            {
                                if (findint(interfering_ctg, interfering_ctg_num, tempNodeNum.node[j]) == 1 &&
                                    findint(interfering_ctg, interfering_ctg_num, key) == 1) {
                                        non_transfer_num++;
                                } else if (findint(interfering_ctg, interfering_ctg_num, tempNodeNum.node[j]) == 1 &&
                                    findint(interfering_ctg, interfering_ctg_num, key) == 0 &&
                                    findint(transfer_ctg, transfer_num, tempNodeNum.node[j]) == 0 &&
                                    ctgdepth[tempNodeNum.node[j] - 1].len > 30) {
                                        if (flag_transfer > 0) {
                                            transfer_ctg[transfer_num] = tempNodeNum.node[j];
                                            transfer_num++;
                                        } else if (ctgdepth[key - 1].len > 2000) {
                                            transfer_ctg[transfer_num] = tempNodeNum.node[j];
                                            transfer_num++;
                                        }
                                } else if (findint(interfering_ctg, interfering_ctg_num, tempNodeNum.node[j]) == 0 &&
                                    findint(interfering_ctg, interfering_ctg_num, key) == 1 &&
                                    findint(transfer_ctg, transfer_num, key) == 0) {
                                        int temp_flag = 0;
                                        nodeNum tempNodeNum2 = kh_value(h_links, kh_get(Ha_nodelink, h_links, tempNodeNum.node[j]));
                                        for (n = 0; n < tempNodeNum2.num; n++) 
                                        {
                                            if (findint(interfering_ctg, interfering_ctg_num, tempNodeNum2.node[n]) == 0)
                                                temp_flag ++;
                                        }
                                        if (temp_flag > 0) {
                                            transfer_ctg[transfer_num] = key;
                                            transfer_num++;
                                        } else if (ctgdepth[tempNodeNum.node[j] - 1].len > 2000) {
                                            transfer_ctg[transfer_num] = key;
                                            transfer_num++;
                                        }
                                }
                            }
                        }
                    }
                }

                for (ka = kh_begin(h_links); ka != kh_end(h_links); ++ka) 
                {
                    if (kh_exist(h_links, ka)) {
                        nodeNum tempNodeNum = kh_value(h_links, ka);
                        free(tempNodeNum.node);
                    }
                }
                kh_destroy(Ha_nodelink, h_links);
                if (transfer_num == 0 && non_transfer_num > 0) {
                    struc--;
                    rm_flag = 1;
                    continue;
                } else if (temp_main_num > 0) {
                    log_info("Structure %d: \n", struc);
                    // if (max_structure_num > 10000000000) {
                    //     log_message(WARNING, "Too many ctg-pairs, Failed to find M-path."); 
                    //     continue;
                    // }
                    int ctg_s = 0;
                    int ctg_len = 0;
                    int utr_s = 5;
                    int utr_e = 3;
                    int temp_linear = 0;

                    for (j = 0; j < temp_mainseeds_num; j++) 
                    {
                        if (ctg_len < ctgdepth[(temp_mainseeds[j]) - 1].len &&
                            findint(interfering_ctg, interfering_ctg_num, temp_mainseeds[j]) == 0) {
                            ctg_s = temp_mainseeds[j];
                            ctg_len = ctgdepth[temp_mainseeds[j] - 1].len;
                        }
                    }
                    if (linear_f == 1) {
                        int temp_len = 0;
                        int temp_ctg = 0;
                        int temp_utr = 0;
                        for (j = 0; j < temp_mainseeds_num; j++) 
                        {
                            int num = 0;
                            int dy_utr = 0;
                            for (n = 0; n < temp_main_num; n++) 
                            {
                                int temp_lutr = temp_mainlinks[n].lutrsmp;
                                int temp_rutr = temp_mainlinks[n].rutrsmp;

                                if (temp_mainlinks[n].lctgsmp == temp_mainseeds[j] && temp_lutr != dy_utr) {
                                    num++;
                                    dy_utr = temp_lutr;
                                } else if (temp_mainlinks[n].rctgsmp == temp_mainseeds[j] && temp_rutr != dy_utr)
                                {
                                    num++;
                                    dy_utr = temp_rutr;
                                }
                                if (num == 2) break;
                            }
                            if (num == 1 && ctgdepth[temp_mainseeds[j] - 1].len > temp_len) {
                                temp_len = ctgdepth[temp_mainseeds[j] - 1].len;
                                temp_ctg = temp_mainseeds[j];
                                temp_utr = dy_utr;
                            }
                        }
                        if (temp_ctg != 0) {ctg_s = temp_ctg; utr_s = (temp_utr == 3) ? 5 : 3; utr_e = (temp_utr == 3) ? 3 : 5; temp_linear = 1;}
                    }
                    int flag_err = 0;
                    float mt_ratio = 0.0;
                    float ctg_s_depth = ctgdepth[ctg_s - 1].depth;
                    /* Filter based on longest mitochondrial contig depth */

                    // maingraph(*bfslinks, mainlinks, ctgdepth, num_dynseeds, *num_bfslinks, *mainseeds, &main_num, mainseeds_num, NULL, 0, 0.5*ctg_s_depth);
                    pathScore struct_path;
                    findMpath(ctg_s, utr_s, ctg_s, utr_e, temp_main_num, temp_mainlinks, ctgdepth, temp_mainseeds, temp_mainseeds_num, interfering_ctg, interfering_ctg_num, &flag_err, &mt_ratio, taxo, &struct_path);
                    if (mt_ratio < 0.7 && temp_linear == 0) {
                        mt_ratio = 0.0;
                        flag_err = 0;
                        findMpath(ctg_s, utr_e, ctg_s, utr_s, temp_main_num, temp_mainlinks, ctgdepth, temp_mainseeds, temp_mainseeds_num, interfering_ctg, interfering_ctg_num, &flag_err, &mt_ratio, taxo, &struct_path);
                    }

                    if (mt_ratio < 0.1 || flag_err == 1) {
                        log_message(WARNING, "Failed to find M-path");
                    } else {
                        float struc_depth = 0.0;
                        for (j = 0; j < struct_path.node_num; j++) {
                            struc_depth += ctgdepth[struct_path.path_node[j] - 1].depth;
                        }
                        struc_depth /= struct_path.node_num;
                        log_info("———————————————————————————————————————\n");
                        log_info(" M-path  Length (bp)  Depth (x)  Score\n");
                        log_info("-------  -----------  ---------  ------\n");
                        log_info(" %-7s %-12ld %-10.2f %-5.2f\n", struct_path.type == 0 ? "C" : "L", struct_path.path_len, struc_depth, mt_ratio * 100);
                        // printf("Graph path %s %.2f%% %ld bp:\n", path_score.type == 0 ? "C" : "L", ratio * 100, path_score.path_len);
                        log_info("\n");
                        log_info("** ");
                        for (j = 0; j < (struct_path.node_num - 1); j++) {
                            // log_info("%d (%c) -> ", struct_path.path_node[j], (struct_path.path_utr[j] == 3 ? '-' : '+'));
                            log_info("%d -> ", struct_path.path_node[j]);
                            if (ass_ctg_num > ass_ctg_mall) {
                                ass_ctg_mall += 100;
                                ass_ctg_arr = realloc(ass_ctg_arr, ass_ctg_mall * sizeof(int));
                            }

                            ass_ctg_arr[ass_ctg_num] = struct_path.path_node[j];
                            ass_ctg_num++;
                        }
                        if (struct_path.type == 0) {
                            log_info("%d\n", struct_path.path_node[struct_path.node_num - 1]);
                        } else {
                            // log_info("%d (%c)\n", struct_path.path_node[struct_path.node_num - 1], (struct_path.path_utr[struct_path.node_num - 1] == 3 ? '-' : '+'));
                            log_info("%d\n", struct_path.path_node[struct_path.node_num - 1]);
                            if (ass_ctg_num > ass_ctg_mall) {
                                ass_ctg_mall += 100;
                                ass_ctg_arr = realloc(ass_ctg_arr, ass_ctg_mall * sizeof(int));
                            }
                            ass_ctg_arr[ass_ctg_num] = struct_path.path_node[j];
                            ass_ctg_num++;
                        }
                        log_info("———————————————————————————————————————\n");

                        ps_num++;
                        if (ps_num > 100) {
                            ps_struct = realloc(ps_struct, (ps_num + 1) * sizeof(pathScore));
                        }
                        ps_struct[ps_num - 1] = struct_path;
                    }


                } else if (temp_main_num == 0 && temp_mainseeds_num == 1) {
                    log_info("Structure %d: \n", struc);
                    log_info("———————————————————————————————————————\n");
                    log_info(" M-path  Length (bp)  Depth (x)  Score\n");
                    log_info("-------  -----------  ---------  ------\n");
                    log_info(" %-7s %-12ld %-10.2f %-5.2f\n", "L", ctgdepth[temp_mainseeds[0] - 1].len, ctgdepth[temp_mainseeds[0] - 1].depth, 100.00);
                    log_info("\n");
                    log_info("** %d (+)\n", temp_mainseeds[0]);
                    log_info("———————————————————————————————————————\n");
                    if (ass_ctg_num > ass_ctg_mall) {
                        ass_ctg_mall += 100;
                        ass_ctg_arr = realloc(ass_ctg_arr, ass_ctg_mall * sizeof(int));
                    }
                    ass_ctg_arr[ass_ctg_num] = temp_mainseeds[0];
                            ass_ctg_num++;

                    pathScore struct_path;
                    struct_path.type = 1;
                    struct_path.path_len = ctgdepth[temp_mainseeds[0] - 1].len;
                    struct_path.node_num = 1;
                    struct_path.path_node = (int*) malloc(1 * sizeof(int));
                    struct_path.path_node[0] = temp_mainseeds[0];
                    struct_path.path_utr = (int*) malloc(1 * sizeof(int));
                    struct_path.path_utr[0] = 5;
                    ps_num++;
                    if (ps_num > 100) {
                        ps_struct = realloc(ps_struct, (ps_num + 1) * sizeof(pathScore));
                    }
                    ps_struct[ps_num - 1] = struct_path;

                } else if (temp_main_num == 0 && temp_mainseeds_num == 0) {
                    struc--;
                }

                // for (int j = 0; j < temp_main_num; j++) 
                // {
                //     free(temp_mainlinks[j].lutr);
                //     free(temp_mainlinks[j].rutr);
                //     free(temp_mainlinks[j].lctg);
                //     free(temp_mainlinks[j].rctg);
                // }
                free(temp_mainlinks);
                free(temp_mainseeds);

            } else {
                log_message(ERROR, "Failed to find structure %d", i);
                exit(EXIT_FAILURE);
            }
        }
        for (i = 0; i < *num_bfslinks; i++)
        {
            if (findint(*mainseeds, *mainseeds_num, (*bfslinks)[i].lctgsmp) == 1 &&
                findint(*mainseeds, *mainseeds_num, (*bfslinks)[i].rctgsmp) == 1) {
                    copy_BFSlinks(&mainlinks[main_num], &(*bfslinks)[i]);
                    main_num++;
                }
        }
    }

    khash_t(Ha_nodeseq)* h_nodeseq = kh_init(Ha_nodeseq);

    if (*mainseeds_num > 0) {
        int numseeds = 0;
        for (i = 0; i < *mainseeds_num; i++) {
            if (findint(interfering_ctg, interfering_ctg_num, (*mainseeds)[i]) == 0) {
                numseeds++;
            }
        }
        log_message(INFO, "Main seeds (%s): %d", organelles_type, numseeds);
    
        /* Assign sequences to main seeds for mt */
        if (strcmp(organelles_type, "mt") == 0) {
            orgAss(exe_path, all_fna, ctgdepth ,output, ass_ctg_arr, ass_ctg_num, "mt", taxo);
        }
        /* main graph */
        FILE* fpmaingfa = fopen(maingfa, "w");
        if (fpmaingfa == NULL) {
            log_message(ERROR, "Failed to open maingfa file %s", maingfa);
            exit(EXIT_FAILURE);
        }
        
        int ret;
        khiter_t iter;
        for (i = 0; i < *mainseeds_num; i++) {
            int seed = (*mainseeds)[i];
            if (rm_flag == 1 && findint(interfering_ctg, interfering_ctg_num, seed) == 1) continue;
            int ctg_RC = ctgdepth[seed - 1].len * ctgdepth[seed - 1].depth;

            for (j = 0; j < num_dynseeds; j++) {
                if (seed == fnainfos[j].ctg) {
                    fprintf(fpmaingfa, "S\t%d\t%s\tLN:i:%d\tRC:i:%d\n", seed, fnainfos[j].seq, ctgdepth[seed - 1].len, ctg_RC);
                    iter = kh_get(Ha_nodeseq, h_nodeseq, seed);
                    if (iter == kh_end(h_nodeseq)) {
                        iter = kh_put(Ha_nodeseq, h_nodeseq, seed, &ret);
                        if (ret == 0) {
                            fprintf(stderr, "Failed to insert into hash table\n");
                            exit(EXIT_FAILURE);
                        }
                        kh_value(h_nodeseq, iter) = strdup(fnainfos[j].seq);
                    }
                }
            }
        }
        for (i = 0; i < main_num; i++) {
            if (rm_flag == 1 && findint(interfering_ctg, interfering_ctg_num, mainlinks[i].lctgsmp) == 1 && findint(interfering_ctg, interfering_ctg_num, mainlinks[i].rctgsmp) == 1) continue;
            fprintf(fpmaingfa, "L\t%d\t", mainlinks[i].lctgsmp);
            if (mainlinks[i].lutrsmp == mainlinks[i].rutrsmp) {
                if (mainlinks[i].lutrsmp == 3) {
                    fprintf(fpmaingfa, "-\t%d\t+\t0M\n", mainlinks[i].rctgsmp);
                } else if (mainlinks[i].lutrsmp == 5) {
                    fprintf(fpmaingfa, "+\t%d\t-\t0M\n", mainlinks[i].rctgsmp);
                }
            } else {
                if (mainlinks[i].lutrsmp == 3) {
                    fprintf(fpmaingfa, "-\t%d\t-\t0M\n", mainlinks[i].rctgsmp);
                } else if (mainlinks[i].lutrsmp == 5) {
                    fprintf(fpmaingfa, "+\t%d\t+\t0M\n", mainlinks[i].rctgsmp);
                }
            }
        }
        fclose(fpmaingfa);
    } else {
        log_message(WARNING, "No main seeds found.");
    }
    /* path to fasta */
    size_t pathfa_out_len = snprintf(NULL, 0, "%s/PMAT_%s.fa", gfa_output, organelles_type) + 1;
    char pathfa[pathfa_out_len];
    snprintf(pathfa, pathfa_out_len, "%s/PMAT_%s.fa", gfa_output, organelles_type);
    path2fa(ps_struct, ps_num, h_nodeseq, pathfa);
    /* free memory */
    // for (int i = 0; i < main_num; i++) {
    //     free(mainlinks[i].lutr);
    //     free(mainlinks[i].rutr);
    //     free(mainlinks[i].lctg);
    //     free(mainlinks[i].rctg);
    // }

    // for (i = 0; i < ps_num; i++) {
    //     free(ps_struct[i].path_node);
    //     free(ps_struct[i].path_utr);
    // }
    // free(ps_struct);
    
    for (k = kh_begin(h_nodeseq); k != kh_end(h_nodeseq); ++k) {
        if (kh_exist(h_nodeseq, k)) {
            free(kh_value(h_nodeseq, k));
        }
    }
    kh_destroy(Ha_nodeseq, h_nodeseq);
    // free(mainlinks);
    // free(mainseeds);
    // free(gfa_output);
    // free(fnainfos);
}