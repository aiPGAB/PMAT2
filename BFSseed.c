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
#include <time.h>
#include <unistd.h> // access
#include <math.h> // sqrt

#include "khash.h"
#include "BFSseed.h"
#include "hitseeds.h"
#include "graphtools.h"
#include "log.h"
#include "misc.h"


typedef struct {
    int* node3;
    int num3;
    int* node5;
    int num5;
} nodeArr;

KHASH_MAP_INIT_INT(Ha_node, nodeArr)

static void BFSmain(const char* type, int num_links, int num_ctg, Ctglinks* ctglinks, CtgDepth* ctgdepth, 
            int* num_dynseeds, int** dynseeds, float nucl_depth , float filter_depth);


void BFSseeds(const char* type, int num_links, int num_ctg, Ctglinks* ctglinks, CtgDepth* ctgdepth, 
            int* num_dynseeds, int** dynseeds, float nucl_depth, float filter_depth, BFSlinks** bfslinks, int* num_bfslinks) {
    
    clock_t start, end;
    int tm = 1;
    
    log_message(INFO, "BFS algorithm starts...");
    log_info("        num\n");
    log_info("  times seeds\n");
    log_info("-------------\n");
    while (1) {
        int temp_num_dynseeds = *num_dynseeds;
        start = clock();
        log_info("  No.%-3d %-4d| ", tm, *num_dynseeds);
        BFSmain(type, num_links, num_ctg, ctglinks, ctgdepth, num_dynseeds, dynseeds, nucl_depth, filter_depth);
        end = clock();
        tm++;
        log_info(" %.2fs\n", (double)(end - start) / CLOCKS_PER_SEC);
        if (tm > 100) {
            log_message(ERROR, "BFS algorithm ends abnormally.");
            exit(EXIT_FAILURE);
        }
        if (temp_num_dynseeds == *num_dynseeds) {
            break;
        }
    }
    log_info("-------------\n");
    
    uint64_t i, j;
    for (i = 0; i < num_links; i++) {
        if (findint(*dynseeds, *num_dynseeds, ctglinks[i].lctg) == 1 && 
            findint(*dynseeds, *num_dynseeds, ctglinks[i].rctg) == 1 &&
            ctglinks[i].linkdepth > 0.3*MIN(ctgdepth[ctglinks[i].lctg - 1].depth, ctgdepth[ctglinks[i].rctg - 1].depth)) 
        {
            (*bfslinks) = realloc(*bfslinks, (*num_bfslinks + 1)*sizeof(BFSlinks));
            if (*bfslinks == NULL) {
                log_message(ERROR, "Failed to allocate memory for BFSlinks");
                exit(EXIT_FAILURE);
            }
            // (*bfslinks)[*num_bfslinks].lctg = strdup(ctgdepth[ctglinks[i].lctg - 1].ctg);
            (*bfslinks)[*num_bfslinks].lutrsmp = remove_quote(ctglinks[i].lutr);
            (*bfslinks)[*num_bfslinks].lctgsmp = ctglinks[i].lctg;
            (*bfslinks)[*num_bfslinks].lctglen = ctgdepth[ctglinks[i].lctg - 1].len;
            
            // (*bfslinks)[*num_bfslinks].rctg = strdup(ctgdepth[ctglinks[i].rctg - 1].ctg);
            (*bfslinks)[*num_bfslinks].rutrsmp = remove_quote(ctglinks[i].rutr);
            (*bfslinks)[*num_bfslinks].rctgsmp = ctglinks[i].rctg;
            (*bfslinks)[*num_bfslinks].rctglen = ctgdepth[ctglinks[i].rctg - 1].len;
            
            (*bfslinks)[*num_bfslinks].linkdepth = ctglinks[i].linkdepth;
            (*num_bfslinks)++;
        }
    }

    while (*num_bfslinks > 0 && *num_dynseeds > 0) {
        /* Remove the bubbles （contg length <= 50bp）*/
        khash_t(Ha_node) *node_hash = kh_init(Ha_node);

        for (i = 0; i < *num_dynseeds; i++) 
        {
            nodeArr node_arr;
            node_arr.num3 = 0;
            node_arr.num5 = 0;
            node_arr.node3 = (int*) malloc((*num_dynseeds + 1)*sizeof(int));
            node_arr.node5 = (int*) malloc((*num_dynseeds + 1)*sizeof(int));

            int ret;
            khint_t k = kh_put(Ha_node, node_hash, (*dynseeds)[i], &ret);

            for (j = 0; j < *num_bfslinks; j++) 
            {
                if ((*bfslinks)[j].lctgsmp == (*dynseeds)[i]) {
                    if ((*bfslinks)[j].lutrsmp == 3) {
                        node_arr.node3[node_arr.num3] = (*bfslinks)[j].rctgsmp;
                        node_arr.num3++;
                    } else if ((*bfslinks)[j].lutrsmp == 5) {
                        node_arr.node5[node_arr.num5] = (*bfslinks)[j].rctgsmp;
                        node_arr.num5++;
                    } else {
                        log_message(ERROR, "Wrong link type: %d", (*bfslinks)[j].lutrsmp);
                    }
                } else if ((*bfslinks)[j].rctgsmp == (*dynseeds)[i]) {
                    if ((*bfslinks)[j].rutrsmp == 3) {
                        node_arr.node3[node_arr.num3] = (*bfslinks)[j].lctgsmp;
                        node_arr.num3++;
                    } else if ((*bfslinks)[j].rutrsmp == 5) {
                        node_arr.node5[node_arr.num5] = (*bfslinks)[j].lctgsmp;
                        node_arr.num5++;
                    } else {
                        log_message(ERROR, "Wrong link type: %d", (*bfslinks)[j].rutrsmp);
                    }
                }
            }
            kh_value(node_hash, k) = node_arr;
        }

        int* del_node = (int*) malloc((*num_dynseeds + 1)*sizeof(int));
        int del_num = 0;

        for (i = 0; i < *num_dynseeds; i++) 
        {
            if (ctgdepth[(*dynseeds)[i] - 1].len <= 50) {
                nodeArr temp_node_arr = kh_value(node_hash, kh_get(Ha_node, node_hash, (*dynseeds)[i]));
                if (temp_node_arr.num3 == 1 && temp_node_arr.num5 == 1) {
                    nodeArr temp_node_arr3 = kh_value(node_hash, kh_get(Ha_node, node_hash, temp_node_arr.node3[0]));
                    nodeArr temp_node_arr5 = kh_value(node_hash, kh_get(Ha_node, node_hash, temp_node_arr.node5[0]));
                    if (
                        ((findint(temp_node_arr3.node3, temp_node_arr3.num3, (*dynseeds)[i]) == 1 && findint(temp_node_arr3.node3, temp_node_arr3.num3, temp_node_arr.node5[0]) == 1) ||
                        (findint(temp_node_arr3.node5, temp_node_arr3.num5, (*dynseeds)[i]) == 1 && findint(temp_node_arr3.node5, temp_node_arr3.num5, temp_node_arr.node5[0]) == 1)) &&

                        ((findint(temp_node_arr5.node3, temp_node_arr5.num3, (*dynseeds)[i]) == 1 && findint(temp_node_arr5.node3, temp_node_arr5.num3, temp_node_arr.node3[0]) == 1) ||
                        (findint(temp_node_arr5.node5, temp_node_arr5.num5, (*dynseeds)[i]) == 1 && findint(temp_node_arr5.node5, temp_node_arr5.num5, temp_node_arr.node3[0]) == 1))
                        
                    ) {
                        del_node[del_num] = (*dynseeds)[i];
                        del_num++;
                    }
                }
            }
        }

        for (i = 0; i < *num_dynseeds; i++) {
            nodeArr temp_node_arr = kh_value(node_hash, kh_get(Ha_node, node_hash, (*dynseeds)[i]));
            if (temp_node_arr.node3 != NULL) {
                free(temp_node_arr.node3);
                temp_node_arr.node3 = NULL;
            }
            if (temp_node_arr.node5 != NULL) {
                free(temp_node_arr.node5);
                temp_node_arr.node5 = NULL;
            }
        }
        kh_destroy(Ha_node, node_hash);

        if (del_num > 0) {
            for (i = 0; i < del_num; i++) 
            {
                remove_element(*dynseeds, num_dynseeds, del_node[i]);
                int temp_num_bfslinks = 0;

                for (j = 0; j < *num_bfslinks; j++) 
                {
                    if ((*bfslinks)[j].lctgsmp != del_node[i] && (*bfslinks)[j].rctgsmp != del_node[i]) {
                        (*bfslinks)[temp_num_bfslinks] = (*bfslinks)[j];
                        temp_num_bfslinks++;
                    }
                }
                *num_bfslinks = temp_num_bfslinks;
            }
        }

        free(del_node);
        if (del_num == 0) break; 
    }
    
    log_message(INFO, "BFS algorithm ends."); 

}


static void BFSmain(const char* type, int num_links, int num_ctg, Ctglinks* ctglinks, CtgDepth* ctgdepth, 
            int* num_dynseeds, int** dynseeds, float nucl_depth, float filter_depth) {

    int progress_step = num_links / 36;
    if (progress_step == 0) {
        progress_step = 1;
    }

    uint64_t i;
    for (i = 0; i < num_links; i++) {
        if (ctgdepth[ctglinks[i].lctg - 1].depth > filter_depth && 
            ctgdepth[ctglinks[i].rctg - 1].depth > filter_depth &&
            ctglinks[i].linkdepth > 0.5*MIN(ctgdepth[ctglinks[i].lctg - 1].depth, ctgdepth[ctglinks[i].rctg - 1].depth)) 
        {
            // log_info("link %d: %s %d %s %g\n", ctglinks[i].lctg, ctglinks[i].lutr, ctglinks[i].rctg, ctglinks[i].rutr, ctglinks[i].linkdepth);
            // log_info("depth of %s: %g\n", ctgdepth[ctglinks[i].lctg - 1].ctg, ctgdepth[ctglinks[i].lctg - 1].depth);
            // log_info("depth of %s: %g\n", ctgdepth[ctglinks[i].rctg - 1].ctg, ctgdepth[ctglinks[i].rctg - 1].depth);

            if (findint(*dynseeds, *num_dynseeds, ctglinks[i].lctg) == 1 &&
                findint(*dynseeds, *num_dynseeds, ctglinks[i].rctg) == 0) 
            {
                *num_dynseeds += 1;
                *dynseeds = realloc(*dynseeds, (*num_dynseeds)*sizeof(int));
                (*dynseeds)[*num_dynseeds - 1] = ctglinks[i].rctg;

            } else if (findint(*dynseeds, *num_dynseeds, ctglinks[i].rctg) == 1 &&
                       findint(*dynseeds, *num_dynseeds, ctglinks[i].lctg) == 0) {
                *num_dynseeds += 1;
                *dynseeds = realloc(*dynseeds, (*num_dynseeds)*sizeof(int));
                (*dynseeds)[*num_dynseeds - 1] = ctglinks[i].lctg;
            }
        }
        if (i % progress_step == 0) {
            // 1000 * 10000
            sleep_ms(10);
            log_info("#");
            fflush(stdout);
        }
    }
}



#ifndef BFSSEED_MAIN
int main(int argc, char* argv[]) {
    if (argc < 7) {
        printf("Usage: %s type all_grap all_fna [output] [cpu] [nucl_depth]\n", argv[0]);
        return 1;
    }
    char *exe_path = realpath(argv[0], NULL);

    const char* all_graph = argv[2];
    if (access(all_graph, F_OK) != 0) {
        log_message(ERROR, "Failed to open file %s", all_graph);
        return -1;    
    }    
    uint64_t i;
    FILE *fp = fopen(all_graph, "r");
    if (fp == NULL) {
        log_message(ERROR, "Failed to open file %s", all_graph);
        return -1;
    }

    int num_links = 0;
    size_t len = 0;
    char *line = NULL;
    int flag_C = 0;
    int num_ctg = 0;
    while (getline(&line, &len, fp) != -1) {
        if (line[0] == 'C') {
            num_links++;
            flag_C = 1;
        } else if (flag_C) {
            break;
        } else {
            num_ctg++;
        }
    }

    rewind(fp);
    Ctglinks* ctglinks = calloc(num_links, sizeof(Ctglinks));
    CtgDepth* ctgdepth = calloc(num_ctg, sizeof(CtgDepth));
    int ctglink_idx = 0;
    int ctg_idx = 0;
    flag_C = 0;
    while (getline(&line, &len, fp) != -1) {
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
            ctgdepth[ctg_id - 1].score = sqrt(sqrt(ctgdepth[ctg_idx].depth) * ctgdepth[ctg_idx].len);
            ctg_idx += 1;
        }
    }
    fclose(fp);
    free(line);

    int num_dynseeds = 0;
    int* dynseeds = calloc(6, sizeof(int));
    hitseeds(exe_path, argv[1], argv[3], argv[4], atoi(argv[5]), num_ctg, ctgdepth, dynseeds, 6);

    for (i = 0; i < 6; i++) {
        if (dynseeds[i] != 0) {
            num_dynseeds++;
        }
    }

    BFSlinks* bfslinks = NULL;
    int num_BFSlinks = 0;
    BFSseeds(argv[1], num_links, num_ctg, ctglinks, ctgdepth, &num_dynseeds, &dynseeds, atof(argv[6]), 2, &bfslinks, &num_BFSlinks);

    // optgfa(num_dynseeds, &dynseeds, bfslinks, num_BFSlinks, ctgdepth, argv[4], argv[3], all_graph);

    /* free memory */
    for (i = 0; i < num_ctg; i++) {
        free(ctgdepth[i].ctg);
    }
    free(ctgdepth);

    for (i = 0; i < num_links; i++) {
        free(ctglinks[i].lutr);
        free(ctglinks[i].rutr);
    }
    
    // for (int i = 0; i < num_BFSlinks; i++) {
    //     free(bfslinks[i].lctg);
    //     free(bfslinks[i].lutr);
    //     free(bfslinks[i].rctg);
    //     free(bfslinks[i].rutr);
    // }

    free(ctglinks);
    free(bfslinks);
    free(dynseeds);
    free(exe_path);

    return 0;
}
#endif