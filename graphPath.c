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
#include <stdbool.h>
#include <pthread.h>

#include "khash.h"
#include "BFSseed.h"
#include "hitseeds.h"
#include "graphtools.h"
#include "misc.h"
#include "log.h"

typedef struct {
    int* node;
    int* utr;
    int type;
    uint32_t nodenum;
    uint64_t nodelen;
    uint64_t pathlen;
} nodePath;

typedef struct {
    int node;
    int utr;
} NodeUtrPair;

typedef struct {
    int* node3;
    int* node3utr;
    int* node5;
    int* node5utr;
    int node3num;
    int node5num;
} nodeLink;

KHASH_MAP_INIT_INT(node_num, int)
KHASH_MAP_INIT_INT(Ha_nodelink, nodeLink)
KHASH_MAP_INIT_INT(Ha_nodedepth, int)
static khash_t(Ha_nodedepth) *g_sort_hash = NULL;
pthread_mutex_t mutex;
uint32_t max_node = 1;
uint32_t mt_uniq = 0;
uint64_t mt_uniq_len = 0;
uint8_t taxo_index = 0;

typedef struct {
    int* path;
    int* utr;
    khash_t(node_num) *h_mito;
    khash_t(node_num) *h_chloro;
    int type;
    uint32_t nodenum;
    uint64_t pathlen;
} mpath;

typedef struct {
    int node1;
    int utr1;
    int node2;
    int utr2;
    int main_num;
    BFSlinks* mainlinks;
    CtgDepth* ctg_depth;
    nodePath* current_path;
    pathScore* path_score;
    int* mt_contigs;
    int mt_num;
    khash_t(node_num) *h_mito;
    int* pt_contigs;
    int pt_num;
    khash_t(node_num) *h_chloro;
    khash_t(Ha_nodelink) *h_links;
    int stop_flag;
} bfs_m_args;


void copy_BFSlinks(BFSlinks* dest, const BFSlinks* src) {
    dest->lctgsmp = src->lctgsmp;
    dest->lctgdepth = src->lctgdepth;
    dest->lctglen = src->lctglen;
    dest->rctgsmp = src->rctgsmp;
    dest->rctgdepth = src->rctgdepth;
    dest->rctglen = src->rctglen;
    dest->linkdepth = src->linkdepth;
    dest->lutrsmp = src->lutrsmp;
    dest->rutrsmp = src->rutrsmp;

    // dest->lctg = src->lctg ? strdup(src->lctg) : NULL;
    // dest->rctg = src->rctg ? strdup(src->rctg) : NULL;
    // dest->lutr = src->lutr ? strdup(src->lutr) : NULL;
    // dest->rutr = src->rutr ? strdup(src->rutr) : NULL;
}

// void free_BFSlinks(BFSlinks *links, int num_links) {
//     for (int i = 0; i < num_links; i++) {
//         if (links[i].lctg) free(links[i].lctg);
//         if (links[i].rctg) free(links[i].rctg);
//         if (links[i].lutr) free(links[i].lutr);
//         if (links[i].rutr) free(links[i].rutr);
//     }
//     free(links);
// }

static void node_recursive(int node, bool* link_used, int* visited_nodes, int* visited_num, int link_num, BFSlinks* links, BFSlinks* tempBFSlinks, int* temp_link_num, int* temp_node, int* temp_node_num) {
    visited_nodes[*visited_num] = node;
    (*visited_num)++;

    uint64_t i;
    for (i = 0; i < link_num; i++) {
        if (link_used[i]) continue;

        if (links[i].lctgsmp == node) {
            copy_BFSlinks(&tempBFSlinks[*temp_link_num], &links[i]);
            (*temp_link_num)++;
            link_used[i] = true;

            if (findint(temp_node, *temp_node_num, links[i].rctgsmp) == 0) {
                temp_node[*temp_node_num] = links[i].rctgsmp;
                (*temp_node_num)++;
            }
            node_recursive(links[i].rctgsmp, link_used, visited_nodes, visited_num, link_num, links, tempBFSlinks, temp_link_num, temp_node, temp_node_num);
        } else if (links[i].rctgsmp == node) {
            copy_BFSlinks(&tempBFSlinks[*temp_link_num], &links[i]);
            (*temp_link_num)++;
            link_used[i] = true;

            if (findint(temp_node, *temp_node_num, links[i].lctgsmp) == 0) {
                temp_node[*temp_node_num] = links[i].lctgsmp;
                (*temp_node_num)++;
            }
            node_recursive(links[i].lctgsmp, link_used, visited_nodes, visited_num, link_num, links, tempBFSlinks, temp_link_num, temp_node, temp_node_num);
        }
    }
}

uint32_t bfs_structure(int node_num, int link_num, BFSlinks* links, int* node_arry, khash_t(Ha_structures)* h_structures) {
    int* visited_nodes = (int*)malloc(2 * (node_num + 1) * sizeof(int));
    bool* link_used = (bool*)malloc((link_num + 1) * sizeof(bool));
    memset(link_used, 0, (link_num + 1) * sizeof(bool));
    int visited_num = 0;
    uint32_t structure_num = 0;

    int ret;
    khint_t k;
    uint64_t i, j;
    
    for (i = 0; i < node_num; i++) {
        int node = node_arry[i];

        if (findint(visited_nodes, visited_num, node) == 0) {
            BFSlinks* tempBFSlinks = (BFSlinks*)malloc((link_num + 1) * sizeof(BFSlinks));
            int temp_link_num = 0;
            int* temp_node = (int*)malloc(node_num * sizeof(int));
            int temp_node_num = 1;
            temp_node[0] = node;

            k = kh_put(Ha_structures, h_structures, structure_num, &ret);
            structure_num++;
            if (ret) {
                node_recursive(node, link_used, visited_nodes, &visited_num, link_num, links, tempBFSlinks, &temp_link_num, temp_node, &temp_node_num);
                BFSstructure *structure = (BFSstructure*)malloc(sizeof(BFSstructure));
                structure->num_links = temp_link_num;
                structure->links = (BFSlinks*)malloc((temp_link_num + 1) * sizeof(BFSlinks));

                for (j = 0; j < temp_link_num; j++) {
                    copy_BFSlinks(&structure->links[j], &tempBFSlinks[j]);
                }
                structure->num_nodes = temp_node_num;
                structure->node = (int*)malloc(temp_node_num * sizeof(int));
                for (j = 0; j < temp_node_num; j++) {
                    structure->node[j] = temp_node[j];
                }
                kh_value(h_structures, k) = structure;
            }
            free(tempBFSlinks);
        }
    }
    free(visited_nodes);
    free(link_used);
    return structure_num;
}

static void bfs_algorithm(int node_s, int s_utr, int node_t, int t_utr, int main_num, BFSlinks* mainlinks, CtgDepth *ctg_depth, nodePath* current_path, nodePath** all_paths, int* path_count) 
{
    
    /* check if the current node is the target node */
    if (node_s == node_t && s_utr != t_utr) {
        *all_paths = (nodePath*)realloc(*all_paths, (*path_count + 1) * sizeof(nodePath));
        (*all_paths)[*path_count] = *current_path;
        (*all_paths)[*path_count].pathlen = current_path->nodelen;
        (*all_paths)[*path_count].type = 0;
        (*path_count)++;
        return;
    }

    /* find all paths from node_s to node_t */
    uint64_t i;
    for (i = 0; i < main_num; i++) {
        if (mainlinks[i].lctgsmp == node_s && s_utr != mainlinks[i].lutrsmp) {
            current_path->nodenum++;
            current_path->node = (int*)realloc(current_path->node, current_path->nodenum * sizeof(int));
            current_path->node[current_path->nodenum - 1] = mainlinks[i].rctgsmp;
            current_path->utr = (int*)realloc(current_path->utr, current_path->nodenum * sizeof(int));
            current_path->utr[current_path->nodenum - 1] = mainlinks[i].rutrsmp;
            current_path->nodelen += ctg_depth[mainlinks[i].rctgsmp - 1].len;

            bfs_algorithm(mainlinks[i].rctgsmp, mainlinks[i].rutrsmp, node_t, t_utr, main_num, mainlinks, ctg_depth, current_path, all_paths, path_count);
            
            /* backtrack */
            current_path->nodenum--;
            current_path->nodelen -= ctg_depth[mainlinks[i].rctgsmp - 1].len;
        }
        else if (mainlinks[i].rctgsmp == node_s && s_utr != mainlinks[i].rutrsmp) {
            current_path->nodenum++;
            current_path->node = (int*)realloc(current_path->node, current_path->nodenum * sizeof(int));
            current_path->node[current_path->nodenum - 1] = mainlinks[i].lctgsmp;
            current_path->utr = (int*)realloc(current_path->utr, current_path->nodenum * sizeof(int));
            current_path->utr[current_path->nodenum - 1] = mainlinks[i].lutrsmp;
            current_path->nodelen += ctg_depth[mainlinks[i].lctgsmp - 1].len;

            bfs_algorithm(mainlinks[i].lctgsmp, mainlinks[i].lutrsmp, node_t, t_utr, main_num, mainlinks, ctg_depth, current_path, all_paths, path_count);

            /* backtrack */
            current_path->nodenum--;
            current_path->nodelen -= ctg_depth[mainlinks[i].lctgsmp - 1].len;
        }
    }
}



void findSpath(int node1, int node1utr, int node2, int node2utr, int main_num, BFSlinks* mainlinks, CtgDepth *ctg_depth) 
{
    nodePath* all_paths = NULL;
    int path_count = 0;

    nodePath current_path;
    current_path.nodenum = 1;
    current_path.node = (int*)malloc(current_path.nodenum * sizeof(int));
    current_path.node[0] = node1;
    current_path.utr = (int*)malloc(current_path.nodenum * sizeof(int));
    current_path.utr[0] = node1utr;
    current_path.nodelen = ctg_depth[node1 - 1].len;
    current_path.pathlen = 0;

    bfs_algorithm(node1, node1utr, node2, node2utr, main_num, mainlinks, ctg_depth, &current_path, &all_paths, &path_count);

    nodePath* shortest_path = &all_paths[0];
    uint64_t i;
    for (i = 1; i < path_count; i++) {
        if (all_paths[i].pathlen < shortest_path->pathlen) {
            shortest_path = &all_paths[i];
        }
    }

    if (shortest_path->nodenum > 1) {
        log_info("-- %d", shortest_path->pathlen);
        for (i = 0; i < shortest_path->nodenum; i++) {
            log_info("%d %d -> ", shortest_path->node[i], shortest_path->utr[i]);
        }
        log_info("\n");
    }

    for (i = 0; i < path_count; i++) {
        free(all_paths[i].node);
        free(all_paths[i].utr);
    }
    free(all_paths);
    free(current_path.node);
    free(current_path.utr);
}


static int path_up(pathScore* path_score, nodePath current_path, CtgDepth *ctg_depth, int mt_ctg, int* mt_ctgarr, int pt_ctg, int* pt_ctgarr)
{
    int mt_num = 0;
    int uniq_mtnum = 0;
    uint64_t uniq_mtlen = 0;
    uint64_t uniq_ptlen = 0;
    int pt_num = 0;
    int uniq_ptnum = 0;
    uint32_t node_num = current_path.nodenum;
    uint64_t path_len = current_path.nodelen;
    int* pathutr = current_path.utr;
    int* path_node = current_path.node;
    int type = current_path.type;

    int* temp_mtarr = (int*)malloc(current_path.nodenum * sizeof(int));
    int* temp_ptarr = (int*)malloc(current_path.nodenum * sizeof(int));

    uint64_t i;
    for (i = 0; i < current_path.nodenum; i++) 
    {
        if (findint(pt_ctgarr, pt_ctg, path_node[i]) == 1 && findint(temp_ptarr, uniq_ptnum, path_node[i]) == 0) {
            temp_ptarr[uniq_ptnum] = path_node[i];
            uniq_ptlen += ctg_depth[path_node[i] - 1].len;
            uniq_ptnum++;
            pt_num++;
        } else if (findint(pt_ctgarr, pt_ctg, path_node[i]) == 1 && findint(temp_ptarr, uniq_ptnum, path_node[i]) == 1) {
            pt_num++;
        } else if (findint(mt_ctgarr, mt_ctg, path_node[i]) == 1 && findint(temp_mtarr, uniq_mtnum, path_node[i]) == 0) {
            temp_mtarr[uniq_mtnum] = path_node[i];
            uniq_mtlen += ctg_depth[path_node[i] - 1].len;
            uniq_mtnum++;
            mt_num++;
        } else if (findint(mt_ctgarr, mt_ctg, path_node[i]) == 1 && findint(temp_mtarr, uniq_mtnum, path_node[i]) == 1) {
            mt_num++;
        }
    }
    free(temp_mtarr);
    free(temp_ptarr);

    pthread_mutex_lock(&mutex);
    path_score->inval_num++;
    int update_flag = 0;
    if (taxo_index == 1) {
        if (uniq_mtlen > path_score->uniq_mt_pathlen) {
            update_flag = 1;
        } else if (uniq_mtlen == path_score->uniq_mt_pathlen && uniq_mtnum > path_score->uniq_mt_nodenum) {
            update_flag = 1;
        } else if (uniq_mtlen == path_score->uniq_mt_pathlen && uniq_mtnum == path_score->uniq_mt_nodenum && type < path_score->type) {
            update_flag = 1;
        } else if (uniq_mtlen == path_score->uniq_mt_pathlen && uniq_mtnum == path_score->uniq_mt_nodenum && type == path_score->type && mt_num > path_score->mt_nodenum) {
            update_flag = 1;
        } else if (uniq_mtlen == path_score->uniq_mt_pathlen && uniq_mtnum == path_score->uniq_mt_nodenum && type == path_score->type && mt_num == path_score->mt_nodenum && path_len < path_score->path_len) {
            update_flag = 1;
        }
    } else {
        if (uniq_mtlen > path_score->uniq_mt_pathlen) {
            update_flag = 1;
        } else if (uniq_mtlen == path_score->uniq_mt_pathlen && uniq_mtnum > path_score->uniq_mt_nodenum) {
            update_flag = 1;
        } else if (uniq_mtlen == path_score->uniq_mt_pathlen && uniq_mtnum == path_score->uniq_mt_nodenum && type < path_score->type) {
            update_flag = 1;
        } else if (uniq_mtlen == path_score->uniq_mt_pathlen && uniq_mtnum == path_score->uniq_mt_nodenum && type == path_score->type && mt_num > path_score->mt_nodenum) {
            update_flag = 1;
        } else if (uniq_mtlen == path_score->uniq_mt_pathlen && uniq_mtnum == path_score->uniq_mt_nodenum && type == path_score->type && mt_num == path_score->mt_nodenum && path_len < path_score->path_len) {
            update_flag = 1;
        } else if (uniq_mtlen == path_score->uniq_mt_pathlen && uniq_mtnum == path_score->uniq_mt_nodenum && type == path_score->type && path_len == path_score->path_len && mt_num == path_score->mt_nodenum && uniq_ptnum < path_score->uniq_pt_nodenum) {
            update_flag = 1;
        } else if (uniq_mtlen == path_score->uniq_mt_pathlen && uniq_mtnum == path_score->uniq_mt_nodenum && type == path_score->type && path_len == path_score->path_len && mt_num == path_score->mt_nodenum && uniq_ptnum == path_score->uniq_pt_nodenum && pt_num < path_score->pt_nodenum) {
            update_flag = 1;
        }
    }

    if (update_flag) {
        path_score->inval_num = 0;
        path_score->uniq_mt_pathlen = uniq_mtlen;
        path_score->uniq_mt_nodenum = uniq_mtnum;
        path_score->mt_nodenum = mt_num;
        path_score->pt_nodenum = pt_num;
        path_score->uniq_pt_nodenum = uniq_ptnum;
        path_score->node_num = node_num;
        float rato = (float)path_score->uniq_mt_pathlen / mt_uniq_len;
        path_score->path_len = path_len;
        path_score->type = type;
        // printf("The length of pt %ld\n", uniq_ptlen);
        // printf("Graph path %s %.2f%% %ld bp:\n-- ", path_score->type == 0 ? "C" : "L", rato * 100, path_score->path_len);
        // printf("Path upate: ");
        for (i = 0; i < node_num; i++) {
            path_score->path_node[i] = path_node[i];
            path_score->path_utr[i] = pathutr[i];
            // printf("%d %d -> ", path_node[i], pathutr[i]);
        }
        // printf("\n");
    }
    pthread_mutex_unlock(&mutex);

    if (path_score->inval_num > 10000000) {
        return 1;
    }
    return 0;
}

static void bfs_m(int node_s, int s_utr, int node_t, int t_utr, int main_num, BFSlinks* mainlinks, CtgDepth *ctg_depth, nodePath* current_path, pathScore* path_score,
    int* mt_contigs, int mt_num, khash_t(node_num) *h_mito, int* pt_contigs, int pt_num, khash_t(node_num) *h_chloro, khash_t(Ha_nodelink) *h_links, int* stop_flag) 
{
    if (*stop_flag) return;
    /* check if the current node is the target node */
    if (current_path->nodenum > 1 && node_s == node_t && s_utr != t_utr) {
        current_path->type = 0;
        current_path->nodelen -= ctg_depth[node_t - 1].len;
        *stop_flag = path_up(path_score, *current_path, ctg_depth, mt_num, mt_contigs, pt_num, pt_contigs);
        current_path->nodelen += ctg_depth[node_t - 1].len;
        return;
    }

    int* temp_mtarr = (int*)malloc(current_path->nodenum * sizeof(int));
    int uniq_mtnum = 0;
    uint64_t uniq_mtlen = 0;
    uint64_t i;
    for (i = 0; i < current_path->nodenum; i++) 
    {
        if (findint(mt_contigs, mt_num, current_path->node[i]) == 1 
            && findint(temp_mtarr, uniq_mtnum, current_path->node[i]) == 0
            && findint(pt_contigs, pt_num, current_path->node[i]) == 0) {
            temp_mtarr[uniq_mtnum] = current_path->node[i];
            uniq_mtlen += ctg_depth[current_path->node[i] - 1].len;
            uniq_mtnum++;
        }
    }
    free(temp_mtarr);
    float rato = (float)uniq_mtlen / mt_uniq_len;
    if (rato > 0.5) {
        *stop_flag = path_up(path_score, *current_path, ctg_depth, mt_num, mt_contigs, pt_num, pt_contigs);
    }

    /* find all paths from node_s to node_t */
    khint_t is_node_s = kh_get(Ha_nodelink, h_links, node_s);
    if (is_node_s == kh_end(h_links)) return;

    nodeLink* link_s = &kh_value(h_links, is_node_s);

    /* Direction checking logic */
    int* node_array = s_utr == 3 ? link_s->node5 : link_s->node3;
    int* utr_array = s_utr == 3 ? link_s->node5utr : link_s->node3utr;
    int node_num = s_utr == 3 ? link_s->node5num : link_s->node3num;
    bool path_stop = true;
    for (i = 0; i < node_num; i++) 
    {
        int next_node = node_array[i];
        int next_utr = utr_array[i];
        bool is_chloro = kh_get(node_num, h_chloro, next_node) != kh_end(h_chloro);
        bool is_mito = kh_get(node_num, h_mito, next_node) != kh_end(h_mito);

        if (is_chloro && kh_value(h_chloro, kh_get(node_num, h_chloro, next_node)) > 0) continue;
        if (is_mito && kh_value(h_mito, kh_get(node_num, h_mito, next_node)) == 0) continue;

        /* Add node to path */
        current_path->node[current_path->nodenum] = next_node;
        current_path->utr[current_path->nodenum] = next_utr;
        current_path->nodelen += ctg_depth[next_node - 1].len;
        current_path->nodenum++;

        if (is_mito) kh_value(h_mito, kh_get(node_num, h_mito, next_node))--;

        if (is_chloro) kh_value(h_chloro, kh_get(node_num, h_chloro, next_node))++;

        path_stop = false;
        /* Recursively call bfs_m */
        bfs_m(next_node, next_utr, node_t, t_utr, main_num, mainlinks, ctg_depth, current_path, path_score, mt_contigs, mt_num, h_mito, pt_contigs, pt_num, h_chloro, h_links, stop_flag);

        /* Backtrack */
        current_path->type = 1;
        current_path->nodenum--;
        current_path->nodelen -= ctg_depth[next_node - 1].len;
        if (is_mito) kh_value(h_mito, kh_get(node_num, h_mito, next_node))++;
        if (is_chloro) kh_value(h_chloro, kh_get(node_num, h_chloro, next_node))--;
    }
    current_path->type = 1;
    /* End if no path is found */
    if (path_stop) *stop_flag = path_up(path_score, *current_path, ctg_depth, mt_num, mt_contigs, pt_num, pt_contigs);
}

static int compare_by_hash_value(const void* a, const void* b) {
    const NodeUtrPair* pair_a = (const NodeUtrPair*)a;
    const NodeUtrPair* pair_b = (const NodeUtrPair*)b;

    int key_a = pair_a->node;
    int key_b = pair_b->node;

    khint_t ka = kh_get(Ha_nodedepth, g_sort_hash, key_a);
    khint_t kb = kh_get(Ha_nodedepth, g_sort_hash, key_b);

    if (ka == kh_end(g_sort_hash) || kb == kh_end(g_sort_hash)) {
        return 0;
    }

    int value_a = kh_value(g_sort_hash, ka);
    int value_b = kh_value(g_sort_hash, kb);

    return value_b - value_a;
}


int kh_copy(khash_t(node_num) *dst, khash_t(node_num) *src) {
    khiter_t k;
    int ret;

    for (k = kh_begin(src); k != kh_end(src); ++k) {
        if (kh_exist(src, k)) {
            int key = kh_key(src, k);
            int value = kh_val(src, k);
            khiter_t d = kh_put(node_num, dst, key, &ret);
            if (ret != -1) {
                kh_val(dst, d) = value;
            }
        }
    }

    return 0;
}

void bfsMap(int nodes, int utrs, int nodet, int utrt, CtgDepth *ctg_depth, mpath* current_path, khash_t(Ha_nodelink) *h_links, int max_paths, khash_t(node_num) *h_mito, khash_t(node_num) *h_chloro, int* path_t)
{
    
    typedef struct {
        mpath mpath[2*max_paths];
        int flag[2*max_paths];
        int path_count;
    } bfsPath;

    bfsPath bfs_path;
    bfs_path.path_count = 1;
    bfs_path.mpath[0].nodenum = 1;
    bfs_path.mpath[0].path = (int*)malloc(max_node * sizeof(int));
    bfs_path.mpath[0].utr = (int*)malloc(max_node * sizeof(int));
    bfs_path.mpath[0].type = 1;
    bfs_path.mpath[0].pathlen = ctg_depth[nodes - 1].len;
    bfs_path.mpath[0].path[0] = nodes;
    bfs_path.mpath[0].utr[0] = utrs;

    bfs_path.mpath[0].h_mito = kh_init(node_num);
    bfs_path.mpath[0].h_chloro = kh_init(node_num);

    kh_copy(bfs_path.mpath[0].h_mito, h_mito);
    kh_copy(bfs_path.mpath[0].h_chloro, h_chloro);

    uint64_t i, j;
    for (i = 0; i < 2*max_paths; i++) {
        bfs_path.flag[i] = 1;
    }

    int dy_path_count = 1;
    int dy_temp_end = 0;
    while (1)
    {
        int temp_end = 0;
        for (i = 0; i < dy_path_count; i++)
        {
            if (bfs_path.flag[i] == 0) continue;
            temp_end = 1;
            int endutr = bfs_path.mpath[i].utr[bfs_path.mpath[i].nodenum - 1];
            int endnode = bfs_path.mpath[i].path[bfs_path.mpath[i].nodenum - 1];
            khiter_t k = kh_get(Ha_nodelink, h_links, endnode);
            if (k == kh_end(h_links)) {
                continue;
            }
            nodeLink* node_links = &kh_val(h_links, k);

            int* tempmap = endutr == 3 ? node_links->node5 : node_links->node3;
            int* temputr = endutr == 3 ? node_links->node5utr : node_links->node3utr;
            int tempmapnum = endutr == 3 ? node_links->node5num : node_links->node3num;

            if (tempmapnum == 1 && bfs_path.flag[i] == 1) {

                if (tempmap[0] == nodet && temputr[0] == utrt) {
                    bfs_path.flag[i] = 0;
                    continue;
                }

                bool is_chloro = kh_get(node_num, bfs_path.mpath[i].h_chloro, tempmap[0]) != kh_end(bfs_path.mpath[i].h_chloro);
                bool is_mito = kh_get(node_num, bfs_path.mpath[i].h_mito, tempmap[0]) != kh_end(bfs_path.mpath[i].h_mito);

                if (is_chloro && kh_value(bfs_path.mpath[i].h_chloro, kh_get(node_num, bfs_path.mpath[i].h_chloro, tempmap[0])) > 1) {bfs_path.flag[i] = 0; continue;}
                if (is_mito && kh_value(bfs_path.mpath[i].h_mito, kh_get(node_num, bfs_path.mpath[i].h_mito, tempmap[0])) == 0) {bfs_path.flag[i] = 0; continue;}
                if (is_mito) kh_value(bfs_path.mpath[i].h_mito, kh_get(node_num, bfs_path.mpath[i].h_mito, tempmap[0]))--;
                if (is_chloro) kh_value(bfs_path.mpath[i].h_chloro, kh_get(node_num, bfs_path.mpath[i].h_chloro, tempmap[0]))++;
                
                bfs_path.mpath[i].nodenum++;
                bfs_path.mpath[i].path = (int*)realloc(bfs_path.mpath[i].path, bfs_path.mpath[i].nodenum * sizeof(int));
                bfs_path.mpath[i].utr = (int*)realloc(bfs_path.mpath[i].utr, bfs_path.mpath[i].nodenum * sizeof(int));
                bfs_path.mpath[i].path[bfs_path.mpath[i].nodenum - 1] = tempmap[0];
                bfs_path.mpath[i].pathlen += ctg_depth[tempmap[0] - 1].len;
                bfs_path.mpath[i].utr[bfs_path.mpath[i].nodenum - 1] = temputr[0];

                if (tempmap[0] == nodet) {
                    bfs_path.flag[i] = 0;
                    bfs_path.mpath[i].type = 0;
                    continue;
                }

                bfs_path.mpath[i].type = 1;

            } else if (tempmapnum > 1 && bfs_path.flag[i] == 1) {
                mpath temp_mpath = bfs_path.mpath[i];
                int dy_path_index = i;
                for (j = 0; j < tempmapnum; j++) 
                {
                    if (tempmap[0] == nodet && temputr[0] == utrt) {
                        bfs_path.flag[i] = 0;
                        break;
                    }
                
                    bool is_chloro = kh_get(node_num, bfs_path.mpath[i].h_chloro, tempmap[j]) != kh_end(bfs_path.mpath[i].h_chloro);
                    bool is_mito = kh_get(node_num, bfs_path.mpath[i].h_mito, tempmap[j]) != kh_end(bfs_path.mpath[i].h_mito);

                    if (is_chloro && kh_value(bfs_path.mpath[i].h_chloro, kh_get(node_num, bfs_path.mpath[i].h_chloro, tempmap[j])) > 1) {bfs_path.flag[dy_path_index] = 0; continue;}
                    if (is_mito && kh_value(bfs_path.mpath[i].h_mito, kh_get(node_num, bfs_path.mpath[i].h_mito, tempmap[j])) == 0) {bfs_path.flag[dy_path_index] = 0; continue;}

                    bfs_path.mpath[dy_path_index].nodenum = temp_mpath.nodenum;
                    bfs_path.mpath[dy_path_index].path = (int*)malloc(max_node * sizeof(int));
                    memcpy(bfs_path.mpath[dy_path_index].path, temp_mpath.path, temp_mpath.nodenum * sizeof(int));
                    bfs_path.mpath[dy_path_index].utr = (int*)malloc(max_node * sizeof(int));
                    memcpy(bfs_path.mpath[dy_path_index].utr, temp_mpath.utr, temp_mpath.nodenum * sizeof(int));
                    bfs_path.mpath[dy_path_index].pathlen = temp_mpath.pathlen;
                    bfs_path.mpath[dy_path_index].type = 1;

                    bfs_path.mpath[dy_path_index].h_mito = kh_init(node_num);
                    bfs_path.mpath[dy_path_index].h_chloro = kh_init(node_num);
                    kh_copy(bfs_path.mpath[dy_path_index].h_mito, temp_mpath.h_mito);
                    kh_copy(bfs_path.mpath[dy_path_index].h_chloro, temp_mpath.h_chloro);

                    if (is_mito) kh_value(bfs_path.mpath[dy_path_index].h_mito, kh_get(node_num, bfs_path.mpath[dy_path_index].h_mito, tempmap[j]))--;
                    if (is_chloro) kh_value(bfs_path.mpath[dy_path_index].h_chloro, kh_get(node_num, bfs_path.mpath[dy_path_index].h_chloro, tempmap[j]))++;
                    
                    bfs_path.mpath[dy_path_index].nodenum++;
                    bfs_path.mpath[dy_path_index].path[bfs_path.mpath[dy_path_index].nodenum - 1] = tempmap[j];
                    bfs_path.mpath[dy_path_index].utr[bfs_path.mpath[dy_path_index].nodenum - 1] = temputr[j];
                    bfs_path.mpath[dy_path_index].pathlen += ctg_depth[tempmap[j] - 1].len;

                    if (tempmap[j] == nodet) {
                        bfs_path.flag[dy_path_index] = 0;
                        bfs_path.mpath[dy_path_index].type = 0;
                    }
                    
                    dy_path_index = bfs_path.path_count;
                    bfs_path.path_count++;
                }
                bfs_path.path_count--;
            } else if (tempmapnum == 0 && bfs_path.flag[i] == 1) {
                bfs_path.flag[i] = 0;
            }
        }

        if (bfs_path.path_count >= max_paths || temp_end == 0) {
            break;
        }
        dy_path_count = bfs_path.path_count;
    }

    *path_t = 0;
    for (i = 0; i < bfs_path.path_count; i++) {
        // int sk = 0;

        // for (int j = 0; j < bfs_path.mpath[i].nodenum; j++)
        // {
        //     bool is_chloro = kh_get(node_num, bfs_path.mpath[i].h_chloro, bfs_path.mpath[i].path[j]) != kh_end(bfs_path.mpath[i].h_chloro);
        //     bool is_mito = kh_get(node_num, bfs_path.mpath[i].h_mito, bfs_path.mpath[i].path[j]) != kh_end(bfs_path.mpath[i].h_mito);

        //     if (is_chloro && kh_value(bfs_path.mpath[i].h_chloro, kh_get(node_num, bfs_path.mpath[i].h_chloro, bfs_path.mpath[i].path[j])) > 1) {
        //         sk = 1;
        //         break;
        //     }
        //     if (is_mito && (kh_value(bfs_path.mpath[i].h_mito, kh_get(node_num, bfs_path.mpath[i].h_mito, bfs_path.mpath[i].path[j])) < 0)) {
        //         sk = 1;
        //         break;
        //     }
        // }
        // if (sk == 1) continue;
        current_path[*path_t].nodenum = bfs_path.mpath[i].nodenum;
        current_path[*path_t].pathlen = bfs_path.mpath[i].pathlen;
        current_path[*path_t].type = bfs_path.mpath[i].type;

        current_path[*path_t].path = (int*)malloc(max_node * sizeof(int));
        current_path[*path_t].utr = (int*)malloc(max_node * sizeof(int));
        
        memcpy(current_path[*path_t].path, bfs_path.mpath[i].path, bfs_path.mpath[i].nodenum * sizeof(int));
        memcpy(current_path[*path_t].utr, bfs_path.mpath[i].utr, bfs_path.mpath[i].nodenum * sizeof(int));

        current_path[*path_t].h_mito = kh_init(node_num);
        current_path[*path_t].h_chloro = kh_init(node_num);
        kh_copy(current_path[*path_t].h_mito, bfs_path.mpath[i].h_mito);
        kh_copy(current_path[*path_t].h_chloro, bfs_path.mpath[i].h_chloro);
        (*path_t)++;
    }

    for (i = 0; i < bfs_path.path_count; i++) {
            free(bfs_path.mpath[i].path);
            free(bfs_path.mpath[i].utr);
            kh_destroy(node_num, bfs_path.mpath[i].h_mito);
            kh_destroy(node_num, bfs_path.mpath[i].h_chloro);
    }
}

void* thread_bfs_m(void* args) {
    bfs_m_args* bfs_args = (bfs_m_args*)args;

    bfs_m(bfs_args->node1, bfs_args->utr1, bfs_args->node2, bfs_args->utr2,
          bfs_args->main_num, bfs_args->mainlinks, bfs_args->ctg_depth,
          bfs_args->current_path, bfs_args->path_score,
          bfs_args->mt_contigs, bfs_args->mt_num, bfs_args->h_mito,
          bfs_args->pt_contigs, bfs_args->pt_num, bfs_args->h_chloro,
          bfs_args->h_links, &bfs_args->stop_flag);
    return NULL;
}

void findMpath(int node1, int node1utr, int node2, int node2utr, int main_num, BFSlinks* mainlinks, CtgDepth *ctg_depth, 
    int* mt_contigs, int mt_num, int* pt_contigs, int pt_num, int* flag_err, float* mt_ratio, int taxo, pathScore* struc_path)
{
    taxo_index = taxo;
    max_node = 1;
    mt_uniq = 0;
    mt_uniq_len = 0;

    int mt_depth = ctg_depth[node1 - 1].depth;
    float min_depth = mt_depth;
    uint64_t i, j;
    for (i = 0; i < mt_num; i++) 
    {
        if (findint(pt_contigs, pt_num, mt_contigs[i]) == 0) {
            if (ctg_depth[mt_contigs[i] - 1].depth < min_depth && ctg_depth[mt_contigs[i] - 1].len > 1000 && ctg_depth[mt_contigs[i] - 1].depth > 0.5*mt_depth) {
                min_depth = ctg_depth[mt_contigs[i] - 1].depth;
            }
        }
    }

    int ret;
    khint_t k;
    
    /* Initialize hash tables for mitochondria and chloroplast contigs */
    khash_t(node_num) *h_mito = kh_init(node_num);
    khash_t(node_num) *h_chloro = kh_init(node_num);
    khash_t(Ha_nodelink) *h_links = kh_init(Ha_nodelink);
    khash_t(Ha_nodedepth) *h_depth = kh_init(Ha_nodedepth);

    /* Calculate max pass count for mitochondria contigs */
    int maxnum_mito = 0;
    // for (int i = 0; i < mt_num; i++) 
    // {
    //     if (findint(pt_contigs, pt_num, mt_contigs[i]) == 0) {
    //         k = kh_put(node_num, h_mito, mt_contigs[i], &ret);
    //         kh_value(h_mito, k) = (int)(ctg_depth[mt_contigs[i] - 1].depth / ((mt_depth + min_depth) / 2) + 0.5);
    //         if (kh_value(h_mito, k) < 1) {
    //             kh_value(h_mito, k) = 1;  
    //         } else if (kh_value(h_mito, k) > 6) {
                
    //             printf("Error: too many mitochondria contigs for node %d\n", node1);
    //             *flag_err = 1;
    //             return;
    //         }
    //         maxnum_mito += kh_value(h_mito, k);
    //         mt_uniq_len += ctg_depth[mt_contigs[i] - 1].len;
    //         // printf("ctg: %d, depth: %f, mt_uniq_len: %ld\n", mt_contigs[i], ctg_depth[mt_contigs[i] - 1].depth, mt_uniq_len);
    //         mt_uniq++;
    //     }
    // }
    for (i = 0; i < mt_num; i++) {
        if (findint(pt_contigs, pt_num, mt_contigs[i]) == 0) {
            int k = kh_put(node_num, h_mito, mt_contigs[i], &ret);
            
            if (ret == 0) {
                continue;
            }
            float denom = (mt_depth + min_depth) / 2;
            if (denom <= 0) {
                log_message(ERROR, "Invalid denominator for contig depth calculation.\n");
                *flag_err = 1;
                return;
            }

            kh_value(h_mito, k) = (int)(ctg_depth[mt_contigs[i] - 1].depth / denom + 0.5);

            if (kh_value(h_mito, k) < 1) {
                kh_value(h_mito, k) = 1;
            // } else if (kh_value(h_mito, k) > 6) {
            //     printf("Error: Too many mitochondria contigs for node %d\n", node1);
            //     *flag_err = 1;
            //     return;
            }

            maxnum_mito += kh_value(h_mito, k);
            mt_uniq_len += ctg_depth[mt_contigs[i] - 1].len;
            mt_uniq++;
        }
    }
    max_node = pt_num*5 + maxnum_mito + 1;


    // for (k = kh_begin(h_mito); k != kh_end(h_mito); k++) 
    // {
    //     if (kh_exist(h_mito, k)) {
    //         log_info("%d: %d; ", kh_key(h_mito, k), kh_value(h_mito, k));
    //     }
    // }
    // log_info("\n");

    /* Initialize chloroplast contig pass count to 0 */
    for (i = 0; i < pt_num; i++) 
    {
        k = kh_put(node_num, h_chloro, pt_contigs[i], &ret);
        kh_value(h_chloro, k) = 0;  // Initially, no chloroplast contig is passed
    }

    for (i = 0; i < mt_num; i++) 
    {
        k = kh_put(Ha_nodedepth, h_depth, mt_contigs[i], &ret);
        kh_value(h_depth, k) = ctg_depth[mt_contigs[i] - 1].depth;
    }
    g_sort_hash = h_depth;

    /* Initialize h_links hash table */
    if (mt_num == 1 && main_num == 1) {
        nodeLink node_links;
        node_links.node3num = 1;
        node_links.node5num = 1;
        node_links.node3 = (int*) malloc(1 * sizeof(int));
        node_links.node3[0] = mt_contigs[0];
        node_links.node3utr = (int*) malloc(1 * sizeof(int));
        node_links.node3utr[0] = 5;
        node_links.node5 = (int*) malloc(1 * sizeof(int));
        node_links.node5[0] = mt_contigs[0];
        node_links.node5utr = (int*) malloc(1 * sizeof(int));
        node_links.node5utr[0] = 3;
        k = kh_put(Ha_nodelink, h_links, mt_contigs[0], &ret);
        kh_value(h_links, k) = node_links;
    } else {
        for (i = 0; i < mt_num; i++) 
        {
            nodeLink node_links;
            node_links.node3num = 0;
            node_links.node5num = 0;
            node_links.node3 = (int*) malloc(10 * sizeof(int));
            node_links.node3utr = (int*) malloc(10 * sizeof(int));
            node_links.node5 = (int*) malloc(10 * sizeof(int));
            node_links.node5utr = (int*) malloc(10 * sizeof(int));

            k = kh_put(Ha_nodelink, h_links, mt_contigs[i], &ret);
            for (j = 0; j < main_num; j++) 
            {
                if (mainlinks[j].lctgsmp == mt_contigs[i] || mainlinks[j].rctgsmp == mt_contigs[i]) {
                    if (mainlinks[j].lctgsmp == mt_contigs[i] && mainlinks[j].lutrsmp == 3) {
                        node_links.node3[node_links.node3num] = mainlinks[j].rctgsmp;
                        node_links.node3utr[node_links.node3num] = mainlinks[j].rutrsmp;
                        node_links.node3num++;
                        if (node_links.node3num >= 10) {
                            *flag_err = 1;
                            return;
                        }
                    } else if (mainlinks[j].rctgsmp == mt_contigs[i] && mainlinks[j].rutrsmp == 3) {
                        node_links.node3[node_links.node3num] = mainlinks[j].lctgsmp;
                        node_links.node3utr[node_links.node3num] = mainlinks[j].lutrsmp;
                        node_links.node3num++;
                        if (node_links.node3num >= 10) {
                            *flag_err = 1;
                            return;
                        }
                    } else if (mainlinks[j].lctgsmp == mt_contigs[i] && mainlinks[j].lutrsmp == 5) {
                        node_links.node5[node_links.node5num] = mainlinks[j].rctgsmp;
                        node_links.node5utr[node_links.node5num] = mainlinks[j].rutrsmp; 
                        node_links.node5num++;
                        if (node_links.node5num >= 10) {
                            *flag_err = 1;
                            return;
                        }
                    } else {
                        node_links.node5[node_links.node5num] = mainlinks[j].lctgsmp;
                        node_links.node5utr[node_links.node5num] = mainlinks[j].lutrsmp;
                        node_links.node5num++;
                        if (node_links.node5num >= 10) {
                            *flag_err = 1;
                            return;
                        }
                    }
                }
            }
            
            int* mttempnode5 = (int*)malloc(node_links.node5num * sizeof(int));
            int* mttempnode5utr = (int*)malloc(node_links.node5num * sizeof(int));
            int mtnode5num = 0;
            int* pttempnode5 = (int*)malloc(node_links.node5num * sizeof(int));
            int* pttempnode5utr = (int*)malloc(node_links.node5num * sizeof(int));
            int ptnode5num = 0;

            int* mttempnode3 = (int*)malloc(node_links.node3num * sizeof(int));
            int* mttempnode3utr = (int*)malloc(node_links.node3num * sizeof(int));
            int mtnode3num = 0;
            int* pttempnode3 = (int*)malloc(node_links.node3num * sizeof(int));
            int* pttempnode3utr = (int*)malloc(node_links.node3num * sizeof(int));
            int ptnode3num = 0;

            for (j = 0; j < node_links.node3num; j++) {
                if (findint(pt_contigs, pt_num, node_links.node3[j]) == 0) {
                    mttempnode3[mtnode3num] = node_links.node3[j];
                    mttempnode3utr[mtnode3num] = node_links.node3utr[j];
                    mtnode3num++;
                } else {
                    pttempnode3[ptnode3num] = node_links.node3[j];
                    pttempnode3utr[ptnode3num] = node_links.node3utr[j];
                    ptnode3num++;
                }
            }


            for (j = 0; j < node_links.node5num; j++) {
                if (findint(pt_contigs, pt_num, node_links.node5[j]) == 0) {
                    mttempnode5[mtnode5num] = node_links.node5[j];
                    mttempnode5utr[mtnode5num] = node_links.node5utr[j];
                    mtnode5num++;
                } else {
                    pttempnode5[ptnode5num] = node_links.node5[j];
                    pttempnode5utr[ptnode5num] = node_links.node5utr[j];
                    ptnode5num++;
                }
            }
            if (mtnode3num > 0) {
                NodeUtrPair* mttempnode3pairs = (NodeUtrPair*)malloc(mtnode3num * sizeof(NodeUtrPair));
                for (j = 0; j < mtnode3num; j++) {
                    mttempnode3pairs[j].node = mttempnode3[j];
                    mttempnode3pairs[j].utr = mttempnode3utr[j];
                }

                qsort(mttempnode3pairs, mtnode3num, sizeof(NodeUtrPair), compare_by_hash_value);

                for (j = 0; j < mtnode3num; j++) {
                    mttempnode3[j] = mttempnode3pairs[j].node;
                    mttempnode3utr[j] = mttempnode3pairs[j].utr;
                }
                free(mttempnode3pairs);
            }
            if (mtnode5num > 0) {
                NodeUtrPair* mttempnode5pairs = (NodeUtrPair*)malloc(mtnode5num * sizeof(NodeUtrPair));
                for (j = 0; j < mtnode5num; j++) {
                    mttempnode5pairs[j].node = mttempnode5[j];
                    mttempnode5pairs[j].utr = mttempnode5utr[j];
                }

                qsort(mttempnode5pairs, mtnode5num, sizeof(NodeUtrPair), compare_by_hash_value);

                for (j = 0; j < mtnode5num; j++) {
                    mttempnode5[j] = mttempnode5pairs[j].node;
                    mttempnode5utr[j] = mttempnode5pairs[j].utr;
                }
                free(mttempnode5pairs);
            }

            for (j = 0; j < node_links.node3num; j++) {
                if (j < mtnode3num) {
                    node_links.node3[j] = mttempnode3[j];
                    node_links.node3utr[j] = mttempnode3utr[j];
                } else {
                    node_links.node3[j] = pttempnode3[j - mtnode3num];
                    node_links.node3utr[j] = pttempnode3utr[j - mtnode3num];
                }
            }
            for (j = 0; j < node_links.node5num; j++) {
                if (j < mtnode5num) {
                    node_links.node5[j] = mttempnode5[j];
                    node_links.node5utr[j] = mttempnode5utr[j];
                } else {
                    node_links.node5[j] = pttempnode5[j - mtnode5num];
                    node_links.node5utr[j] = pttempnode5utr[j - mtnode5num];
                }
            }
            kh_value(h_links, k) = node_links;
        }
    }

    pathScore path_score;
    path_score.mt_nodenum = 0;
    path_score.node_num = 0;
    path_score.uniq_mt_nodenum = 0;
    path_score.path_len = 0;
    path_score.uniq_mt_pathlen = 0;
    path_score.pt_nodenum = 10 * pt_num;
    path_score.uniq_pt_nodenum = pt_num;
    path_score.path_node = (int*)malloc((max_node) * sizeof(int));
    path_score.path_utr = (int*)malloc((max_node) * sizeof(int));
    path_score.type = 1;
    path_score.inval_num = 0;

    // nodePath current_path;
    // current_path.nodenum = 1;
    // current_path.node = (int*)malloc((max_node) * sizeof(int));
    // current_path.node[0] = node1;
    // current_path.utr = (int*)malloc((max_node) * sizeof(int));
    // current_path.utr[0] = node1utr;
    // current_path.nodelen = ctg_depth[node1 - 1].len;
    // current_path.pathlen = 0;
    // current_path.type = 1;

    int cpu = 8;
    mpath* bfs_path = (mpath*)malloc(2*cpu * sizeof(mpath));
    int path_t = 0;
    bfsMap(node1, node1utr, node2, node2utr, ctg_depth, bfs_path, h_links, cpu, h_mito, h_chloro, &path_t);
    
    // printf("Path threads: %d\n", path_t);
    // for (int i = 0; i < path_t; i++) 
    // {
    //     printf("Path threads %d: ", i);
    //     for (int j = 0; j < bfs_path[i].nodenum; j++)
    //     {
    //         printf("%d %d ", bfs_path[i].path[j], bfs_path[i].utr[j]);
    //     }
    //     printf("\n");

    //     // for (k = kh_begin(bfs_path[i].h_mito); k != kh_end(bfs_path[i].h_mito); k++) 
    //     // {
    //     //     if (kh_exist(bfs_path[i].h_mito, k)) {
    //     //         log_info("%d: %d; ", kh_key(bfs_path[i].h_mito, k), kh_value(bfs_path[i].h_mito, k));
    //     //     }
    //     // }
    //     // log_info("\n");

    // }

    pthread_t threads[path_t];
    bfs_m_args thread_args[path_t];

    nodePath current_path[path_t];
    
    if (pthread_mutex_init(&mutex, NULL) != 0) {
        log_message(ERROR, "Error: mutex init failed\n");
        *flag_err = 1;
        return;
    }

    for (i = 0; i < path_t; i++) {
        int node_s = bfs_path[i].path[bfs_path[i].nodenum - 1];
        int utr_s = bfs_path[i].utr[bfs_path[i].nodenum - 1];
        current_path[i].nodenum = bfs_path[i].nodenum;
        current_path[i].node = (int*)malloc((max_node) * sizeof(int));
        current_path[i].utr = (int*)malloc((max_node) * sizeof(int));
        for (j = 0; j < bfs_path[i].nodenum; j++) {
            current_path[i].node[j] = bfs_path[i].path[j];
            current_path[i].utr[j] = bfs_path[i].utr[j];
        }
        current_path[i].nodelen = bfs_path[i].pathlen;
        current_path[i].type = bfs_path[i].type;

        thread_args[i].node1 = node_s;
        thread_args[i].utr1 = utr_s;
        thread_args[i].node2 = node2;
        thread_args[i].utr2 = node2utr;
        thread_args[i].main_num = main_num;
        thread_args[i].mainlinks = mainlinks;
        thread_args[i].ctg_depth = ctg_depth;
        thread_args[i].current_path = &(current_path[i]);
        thread_args[i].path_score = &path_score;
        thread_args[i].mt_contigs = mt_contigs;
        thread_args[i].mt_num = mt_num;
        thread_args[i].h_mito = bfs_path[i].h_mito;
        thread_args[i].pt_contigs = pt_contigs;
        thread_args[i].pt_num = pt_num;
        thread_args[i].h_chloro = bfs_path[i].h_chloro;
        thread_args[i].h_links = h_links;
        thread_args[i].stop_flag = 0;

        pthread_create(&threads[i], NULL, thread_bfs_m, (void*)&thread_args[i]);

    }
    for (i = 0; i < path_t; i++) {
        pthread_join(threads[i], NULL);
    }
    free(bfs_path);

    pthread_mutex_destroy(&mutex);
    // bfs_m(node1, node1utr, node2, node2utr, main_num, mainlinks, ctg_depth, &current_path, &path_score, mt_contigs, mt_num, h_mito, pt_contigs, pt_num, h_chloro, h_links);
    
    if (path_score.node_num == 0) {
        *flag_err = 1;
        // free(current_path.node);
        // free(current_path.utr);
        free(path_score.path_node);
        free(path_score.path_utr);
        kh_destroy(node_num, h_mito);
        kh_destroy(node_num, h_chloro);
        for (i = 0; i < path_t; i++)
        {
            free(current_path[i].node);
            free(current_path[i].utr);
        }
        return;
    }

    for (i = 0; i < path_t; i++)
    {
        free(current_path[i].node);
        free(current_path[i].utr);
    }

    /* print the paths */
    float ratio = (float)path_score.uniq_mt_pathlen / mt_uniq_len;
    *mt_ratio = ratio;
    if (ratio >= 0.1) {
        // log_info("————————————————————————————\n");
        // log_info(" M-path  Length (bp)  Score\n");
        // log_info("----------------------------\n");
        // log_info(" %-7s %-12ld %-5.2f\n", path_score.type == 0 ? "C" : "L", path_score.path_len, ratio * 100);
        // // printf("Graph path %s %.2f%% %ld bp:\n", path_score.type == 0 ? "C" : "L", ratio * 100, path_score.path_len);
        // log_info("\n");
        // log_info("** %d %d -> ", path_score.path_node[0], node1utr == 3 ? 5 : 3);
        // for (int i = 1; i < (path_score.node_num - 1); i++) {
        //     log_info("%d %d -> ", path_score.path_node[i], path_score.path_utr[i]);
        // }
        // log_info("%d %d\n", path_score.path_node[path_score.node_num - 1], path_score.path_utr[path_score.node_num - 1]);
        // log_info("————————————————————————————\n");

        struc_path->type = path_score.type;
        struc_path->path_len = path_score.path_len;
        struc_path->node_num = path_score.node_num;
        struc_path->uniq_mt_pathlen = path_score.uniq_mt_pathlen;
        struc_path->uniq_mt_nodenum = path_score.uniq_mt_nodenum;
        struc_path->uniq_pt_nodenum = path_score.uniq_pt_nodenum;
        struc_path->pt_nodenum = path_score.pt_nodenum;
        struc_path->mt_nodenum = path_score.mt_nodenum;
        
        struc_path->path_node = (int*)malloc((path_score.node_num) * sizeof(int));
        struc_path->path_utr = (int*)malloc((path_score.node_num) * sizeof(int));
        for (i = 0; i < path_score.node_num; i++) {
            if (i == 0) {
                struc_path->path_node[i] = node1;
                struc_path->path_utr[i] = node1utr;
            } else {
                struc_path->path_node[i] = path_score.path_node[i];
                struc_path->path_utr[i] = path_score.path_utr[i];
            }
        }

        // for (int i = 0; i < mt_num; i++) {
        //     if (findint(path_score.path_node, path_score.node_num, mt_contigs[i]) == 0 && 
        //         findint(pt_contigs, pt_num, mt_contigs[i]) == 0) {
        //         printf("%d ", mt_contigs[i]);
        //     }
        // }
        // printf("\n");
    }

    /* Free memory */
    // free(current_path.node);
    // free(current_path.utr);
    free(path_score.path_node);
    free(path_score.path_utr);
    kh_destroy(node_num, h_mito);
    kh_destroy(node_num, h_chloro);
    kh_destroy(Ha_nodelink, h_links);
    kh_destroy(Ha_nodedepth, h_depth);
}

