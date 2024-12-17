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

#include "misc.h"
#include "log.h"
#include "graphtools.h"
#include "khash.h"


/*  Function to reverse complement a sequence */
char *convert_seq(const char *seq, int flag) {
    uint64_t i;
    if (flag == 3) {
        return strdup(seq);
    } else if (flag == 5) {
        size_t len = strlen(seq);
        char *rev_comp = (char *)malloc(len + 1);
        if (!rev_comp) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }

        for (i = 0; i < len; i++) {
            char base = seq[len - 1 - i]; // Reverse the sequence
            switch (base) {
                case 'A': rev_comp[i] = 'T'; break;
                case 'T': rev_comp[i] = 'A'; break;
                case 'C': rev_comp[i] = 'G'; break;
                case 'G': rev_comp[i] = 'C'; break;
                case 'a': rev_comp[i] = 'T'; break;
                case 't': rev_comp[i] = 'A'; break;
                case 'c': rev_comp[i] = 'G'; break;
                case 'g': rev_comp[i] = 'C'; break;
                default:  rev_comp[i] = base; break; 
            }
        }
        rev_comp[len] = '\0';
        return rev_comp;
    }
    return NULL;
}


void path2fa(pathScore *path, int ps_num, khash_t(Ha_nodeseq)* node_seq, const char *output) {
    FILE *ow = fopen(output, "w");
    if (!ow) {
        fprintf(stderr, "Failed to open output file\n");
        exit(EXIT_FAILURE);
    }
    uint64_t i, j;
    for (i = 0; i < ps_num; i++) {
        pathScore *ps = &path[i];
        char *final_seq = (char *)malloc(1024);
        if (!final_seq) {
            log_message(ERROR, "Failed to process path %d", i);
            fclose(ow);
            exit(EXIT_FAILURE);
        }
        size_t final_len = 0;

        int ps_node_num = ps->node_num;
        if (ps->type == 0) {
            ps_node_num--; // Remove the start node
        }
        for (j = 0; j < ps_node_num; j++) {
            int node_id = ps->path_node[j];
            int orientation = ps->path_utr[j];

            khiter_t iter = kh_get(Ha_nodeseq, node_seq, node_id);
            if (iter == kh_end(node_seq)) {
                log_message(ERROR, "Node sequence not found for node %d", node_id);
                fclose(ow);
                free(final_seq);
                exit(EXIT_FAILURE);
            }
            char *nodeseq = kh_value(node_seq, iter);
            char *converted_seq = convert_seq(nodeseq, orientation);
            // printf("Converted sequence %s\n", converted_seq);
            if (!converted_seq) {
                log_message(ERROR, "Failed to process node %d", node_id);
                fclose(ow);
                free(final_seq);
                exit(EXIT_FAILURE);
            }

            size_t converted_len = strlen(converted_seq);
            final_seq = (char *)realloc(final_seq, final_len + converted_len + 1);
            if (!final_seq) {
                fprintf(stderr, "Memory allocation failed for final_seq\n");
                fclose(ow);
                free(converted_seq);
                exit(EXIT_FAILURE);
            }
            memcpy(final_seq + final_len, converted_seq, converted_len);
            final_len += converted_len;
            final_seq[final_len] = '\0';

            free(converted_seq);
        }

        size_t buffer_size = 256;
        size_t offset = 0;
        char *seq_id = (char *)malloc(buffer_size);
        if (!seq_id) {
            fprintf(stderr, "Memory allocation failed for seq_id\n");
            fclose(ow);
            free(final_seq);
            exit(EXIT_FAILURE);
        }

        const char *prefix = (ps->type == 1) ? "L" : "C";
        snprintf(seq_id, buffer_size, "%s_", prefix);
        offset = strlen(seq_id);

        for (j = 0; j < ps->node_num; j++) {
            int node_id = ps->path_node[j];
            char utr_flag = (ps->path_utr[j] == 5) ? '+' : '-'; // 5 -> +, 3 -> -
            size_t needed_space = snprintf(NULL, 0, "%d%c", node_id, utr_flag);

            if (offset + needed_space + 1 > buffer_size) {
                buffer_size *= 2;
                char *new_seq_id = (char *)realloc(seq_id, buffer_size);
                if (!new_seq_id) {
                    fprintf(stderr, "Memory reallocation failed for seq_id\n");
                    fclose(ow);
                    free(seq_id);
                    free(final_seq);
                    exit(EXIT_FAILURE);
                }
                seq_id = new_seq_id;
            }

            int written = snprintf(seq_id + offset, buffer_size - offset, "%d%c", node_id, utr_flag);
            if (written < 0) {
                fprintf(stderr, "Error writing to seq_id\n");
                fclose(ow);
                free(seq_id);
                free(final_seq);
                exit(EXIT_FAILURE);
            }

            offset += written;
        }

        seq_id[offset] = '\0';

        // Write to the output file
        fprintf(ow, ">%s\n%s\n", seq_id, final_seq);

        // Clean up memory for this path's final sequence and ID
        free(seq_id);
        free(final_seq);
    }

    // Clean up the output file
    fclose(ow);
}