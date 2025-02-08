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
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>
#include <libgen.h>

#include "log.h"
#include "misc.h"
#include "seqtools.h"

// void remove_file(const char *filepath);
int remove_dir(const char *path);
void remove_prefix_files(const char *dir, const char *prefix);
void clean_directory(const char *dir);


void run_Assembly(const char *sif_path, int cpu, const char *assembly_seq, const char *output_path, int mi, int ml, int mem, float genomesize_bp) {

    log_message(INFO, "Reads assembly start...");

    mkdirfiles(output_path);

    // Allocate memory for runAssembly_output
    char* runAssembly_output = (char*)malloc(sizeof(*runAssembly_output) * (snprintf(NULL, 0, "%s/assembly_result", output_path) + 1));
    sprintf(runAssembly_output, "%s/assembly_result", output_path);

    // Resolve absolute path
    char *absolute_assembly_seq = realpath(assembly_seq, NULL);
    if (!absolute_assembly_seq) {
        log_message(ERROR, "Error resolving absolute path of assembly_seq");
        free(runAssembly_output);
        exit(EXIT_FAILURE);
    }

    // Allocate memory for mount_output
    char* mount_output = (char*)malloc(sizeof(*mount_output) * (snprintf(NULL, 0, "/data/%s", output_path) + 1));
    sprintf(mount_output, "/data/%s", output_path);

    // Allocate memory for bindpath
    char* bindpath = (char*)malloc(sizeof(*bindpath) * (snprintf(NULL, 0, "%s,%s:%s", 
                     absolute_assembly_seq, output_path, mount_output) + 1));
    sprintf(bindpath, "%s,%s:%s", absolute_assembly_seq, output_path, mount_output);

    char* command = NULL;
    if (which_executable("apptainer") == 1) {
        if (setenv("APPTAINER_BINDPATH", bindpath, 1) != 0) {
            log_message(ERROR, "Error setting APPTAINER_BINDPATH");
            exit(EXIT_FAILURE);
        }

        size_t cmd_len = snprintf(NULL, 0, "setsid apptainer exec %s runAssembly -cpu %d -het -force -sio -urt -large -s 100 -m -nobig -mi %d -ml %d -o %s %s", 
            sif_path, cpu, mi, ml, runAssembly_output, absolute_assembly_seq) + 1;
        command = malloc(cmd_len);

        if (mem) {
            log_info("Running command:\n"
                    " apptainer exec %s \\\n"
                    "     runAssembly \\\n"
                    "     -cpu %d -het -force -sio -m \\\n"
                    "     -urt -large -s 100 -nobig -mi %d \\\n"
                    "     -ml %d -o %s %s\n\n",
                    sif_path, cpu, mi, ml, runAssembly_output, absolute_assembly_seq);
            snprintf(command, cmd_len,
                    "setsid apptainer exec %s runAssembly -cpu %d -het -force -sio -m -urt -large -s 100 -nobig -mi %d -ml %d -o %s %s",
                    sif_path, cpu, mi, ml, runAssembly_output, absolute_assembly_seq);
        } else {
            log_info("Running command:\n"
                    " apptainer exec %s \\\n"
                    "     runAssembly \\\n"
                    "     -cpu %d -het -force -sio -urt -large -s 100 -nobig -mi %d \\\n"
                    "     -ml %d -o %s %s\n\n",
                    sif_path, cpu, mi, ml, runAssembly_output, absolute_assembly_seq);
            snprintf(command, cmd_len,
                    "setsid apptainer exec %s runAssembly -cpu %d -het -force -sio -urt -large -s 100 -nobig -mi %d -ml %d -o %s %s",
                    sif_path, cpu, mi, ml, runAssembly_output, absolute_assembly_seq);
        }

    } else if (which_executable("singularity") == 1) {
        if (setenv("SINGULARITY_BINDPATH", bindpath, 1) != 0) {
            log_message(ERROR, "Error setting SINGULARITY_BINDPATH");
            exit(EXIT_FAILURE);
        }

        size_t cmd_len = snprintf(NULL, 0, "setsid singularity exec %s runAssembly -cpu %d -het -force -sio -urt -large -s 100 -m -nobig -mi %d -ml %d -o %s %s", 
            sif_path, cpu, mi, ml, runAssembly_output, absolute_assembly_seq) + 1;
        command = malloc(cmd_len);
        if (mem) {
            log_info("Running command:\n"
                    " singularity exec %s \\\n"
                    "     runAssembly \\\n"
                    "     -cpu %d -het -force -sio -m \\\n"
                    "     -urt -large -s 100 -nobig -mi %d \\\n"
                    "     -ml %d -o %s %s\n\n",
                    sif_path, cpu, mi, ml, runAssembly_output, absolute_assembly_seq);
            snprintf(command, cmd_len,
                    "setsid singularity exec %s runAssembly -cpu %d -het -force -sio -m -urt -large -s 100 -nobig -mi %d -ml %d -o %s %s",
                    sif_path, cpu, mi, ml, runAssembly_output, absolute_assembly_seq);
        } else {
            log_info("Running command:\n"
                    " singularity exec %s \\\n"
                    "     runAssembly \\\n"
                    "     -cpu %d -het -force -sio -urt -large -s 100 -nobig -mi %d \\\n"
                    "     -ml %d -o %s %s\n\n",
                    sif_path, cpu, mi, ml, runAssembly_output, absolute_assembly_seq);
            snprintf(command, cmd_len,
                    "setsid singularity exec %s runAssembly -cpu %d -het -force -sio -urt -large -s 100 -nobig -mi %d -ml %d -o %s %s",
                    sif_path, cpu, mi, ml, runAssembly_output, absolute_assembly_seq);
        }
    } else {
        log_message(ERROR, "Neither apptainer nor singularity is installed.");
        exit(EXIT_FAILURE);
    }

    // int go_flag = ass_command(command, 0, 1);
    // int tm = 0;
    // while(1)
    // {
    //     if (go_flag == 0 || tm > 3) break;
        
    //     long longassembly_bp = getFileSize(assembly_seq);
    //     float seq_depth = (float)longassembly_bp / (float)genomesize_bp;
    //     float subsize = 0.0;
    //     if (seq_depth >= 5) {
    //         subsize = (seq_depth + 5) / (2 * seq_depth);
    //     } else if (seq_depth >= 2) {
    //         subsize = (seq_depth + 2) / (2 * seq_depth);
    //     } else {
    //         log_message(ERROR, "The sequence depth is too low.");
    //         delete_directory(runAssembly_output);
    //         free(dir);
    //         exit(EXIT_FAILURE);
    //     }

    //     char* new_cut_seq = malloc(strlen(assembly_seq) + strlen(".bak") + 1);
    //     if (new_cut_seq == NULL) {
    //         fprintf(stderr, "Memory allocation failed!\n");
    //         exit(EXIT_FAILURE);
    //     }
    //     snprintf(new_cut_seq, strlen(assembly_seq) + strlen(".bak") + 1, 
    //         "%s.bak", assembly_seq);
    //     subsample(new_cut_seq, assembly_seq, subsize, 6);
    //     tm++;
    //     remove_file(assembly_seq);
    //     rename_file(new_cut_seq, assembly_seq);
    //     free(new_cut_seq);

    //     go_flag = ass_command(command, 0, 0);
    // }
    // free(command);

    int go_flag = ass_command(command, 0, 1);
    int tm = 0;
    while(1)
    {
        if (go_flag != -2 || tm > 3) break;
        float subsize = 0.7;

        char* new_cut_seq = (char*)malloc(sizeof(*new_cut_seq) * (snprintf(NULL, 0, "%s.bak", assembly_seq) + 1));
        sprintf(new_cut_seq, "%s.bak", assembly_seq);
        subsample(new_cut_seq, assembly_seq, subsize, 6);
        tm++;
        remove_file(assembly_seq);
        rename_file(new_cut_seq, assembly_seq);
        free(new_cut_seq);

        go_flag = ass_command(command, 0, 0);
    }
    free(command);

    // 
    size_t contig_fna_len = snprintf(NULL, 0, "%s/454AllContigs.fna", runAssembly_output) + 1;
    char contig_fna[contig_fna_len];
    size_t contig_graph_len = snprintf(NULL, 0, "%s/454ContigGraph.txt", runAssembly_output) + 1;
    char contig_graph[contig_graph_len];
    snprintf(contig_fna, sizeof(contig_fna), "%s/454AllContigs.fna", runAssembly_output);
    snprintf(contig_graph, sizeof(contig_graph), "%s/454ContigGraph.txt", runAssembly_output);

    if (access(contig_fna, F_OK) == 0 && access(contig_graph, F_OK) == 0) {
        size_t new_fna_len = snprintf(NULL, 0, "%s/PMATAllContigs.fna", runAssembly_output) + 1;
        char new_fna[new_fna_len];
        size_t new_graph_len = snprintf(NULL, 0, "%s/PMATContigGraph.txt", runAssembly_output) + 1;
        char new_graph[new_graph_len];
        snprintf(new_fna, sizeof(new_fna), "%s/PMATAllContigs.fna", runAssembly_output);
        snprintf(new_graph, sizeof(new_graph), "%s/PMATContigGraph.txt", runAssembly_output);

        rename(contig_fna, new_fna);
        rename(contig_graph, new_graph);

        remove_prefix_files(runAssembly_output, "454");
        
        char* sff_dir = (char*)malloc(sizeof(*sff_dir) * (snprintf(NULL, 0, "%s/sff", runAssembly_output) + 1));
        sprintf(sff_dir, "%s/sff", runAssembly_output);
        rmdir(sff_dir);
        
    } else {
        remove_prefix_files(runAssembly_output, "454");
        
        char* sff_dir = (char*)malloc(sizeof(*sff_dir) * (snprintf(NULL, 0, "%s/sff", runAssembly_output) + 1));
        sprintf(sff_dir, "%s/sff", runAssembly_output);
        rmdir(sff_dir);
        
        clean_directory(runAssembly_output);
        log_message(ERROR, "The assembly failed.");
        exit(EXIT_FAILURE);
    }

    free(runAssembly_output);
    free(absolute_assembly_seq);
    free(mount_output);
    free(bindpath);

    log_message(INFO, "Reads assembly end.");
}

// void remove_file(const char *filepath) {
//     if (remove(filepath) != 0) {
//         perror("Error deleting file");
//     }
// }

int remove_dir(const char *path) {
    struct dirent *entry;
    DIR *dir = opendir(path);

    if (dir == NULL) {
        log_message( ERROR, "Unable to open directory");
        free(dir);
        exit(EXIT_FAILURE);
    }

    while ((entry = readdir(dir)) != NULL) {
        char full_path[4096];

        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }

        snprintf(full_path, sizeof(full_path), "%s/%s", path, entry->d_name);

        struct stat statbuf;
        if (stat(full_path, &statbuf) == 0) {
            if (S_ISDIR(statbuf.st_mode)) {
                remove_dir(full_path);
            } else {
                remove_file(full_path);
            }
        }
    }

    closedir(dir);
    // rmdir(path);  // Remove the directory itself

    return 0;
}

void remove_prefix_files(const char *dir, const char *prefix) {
    DIR *dirp = opendir(dir);
    if (!dirp) {
        return;
    }

    struct dirent *entry;
    while ((entry = readdir(dirp)) != NULL) {
        if (strncmp(entry->d_name, prefix, strlen(prefix)) == 0) {
            char filepath[4096];
            snprintf(filepath, sizeof(filepath), "%s/%s", dir, entry->d_name);
            remove_file(filepath);
        }
    }

    closedir(dirp);
}

void clean_directory(const char *dir) {
    DIR *dirp = opendir(dir);
    if (!dirp) {
        return;
    }

    struct dirent *entry;
    while ((entry = readdir(dirp)) != NULL) {
        char filepath[4096];

        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }

        snprintf(filepath, sizeof(filepath), "%s/%s", dir, entry->d_name);
        struct stat statbuf;
        if (stat(filepath, &statbuf) == 0) {
            if (S_ISDIR(statbuf.st_mode)) {
                remove_dir(filepath);
            } else {
                remove_file(filepath);
            }
        }
    }

    closedir(dirp);
}

#ifndef RUNASSEMBLY_MAIN
int main(int argc, char *argv[]) {
    if (argc < 4) {
        log_message(ERROR, "Usage: %s input_seq output_path cpu", argv[0]);
        exit(EXIT_FAILURE);
    }

    char *exe_path = realpath(argv[0], NULL);
    char *assembly_seq = argv[1];
    char *output_path = argv[2];
    int cpu = atoi(argv[3]);
    int mi = 90;
    int ml = 40;
    int mem = 0;

    // run_Assembly(exe_path, cpu, assembly_seq, output_path, mi, ml, mem);

    free(exe_path);

    return 0;
}
#endif