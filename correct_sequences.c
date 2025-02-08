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
#include <sys/stat.h>
#include <dirent.h>
#include <libgen.h>

#include "log.h"
#include "misc.h"
#include "seqtools.h"


static void canu_trim(const char* canu_path, const char* input_seq, int genomeSize, const char* output_dir, const char* readstype, int cpu);
static void config_info(const char* cfgr, const char* cfgw, char* workdir, int cfg_flag, char* readstype, char* input_fofn, int genomeSize, int cpu);
static void cns_files(const char* cns_dir, const char* output_file);

void canu_correct(const char* canu_path, const char* input_seq, int genomeSize, const char* output_dir, const char* readstype, int cpu) {
    log_message(INFO, "Start canu correction...");

    char* output_correct = (char*)malloc(sizeof(*output_correct) * (snprintf(NULL, 0, "%s/correct_out", output_dir) + 1));
    sprintf(output_correct, "%s/correct_out", output_dir);

    mkdirfiles(output_correct);

    char* canu_cmd = (char*)malloc(sizeof(*canu_cmd) * (snprintf(NULL, 0, "%s -correct -p PMAT -d %s genomeSize=%d useGrid=false maxThreads=%d -%s %s", 
                                  canu_path, output_correct, genomeSize, cpu, readstype, input_seq) + 1));
    sprintf(canu_cmd, "%s -correct -p PMAT -d %s genomeSize=%d useGrid=false maxThreads=%d -%s %s", 
            canu_path, output_correct, genomeSize, cpu, readstype, input_seq);

    execute_command(canu_cmd, 0, 1);

    char* corrected_reads_path = (char*)malloc(sizeof(*corrected_reads_path) * (snprintf(NULL, 0, "%s/PMAT.correctedReads.fasta.gz", output_correct) + 1));
    sprintf(corrected_reads_path, "%s/PMAT.correctedReads.fasta.gz", output_correct);

    if (access(corrected_reads_path, F_OK) != 0) {
        log_message(ERROR, "An error occurred during the correction process?");
        free(output_correct);
        free(canu_cmd);
        free(corrected_reads_path);
        exit(EXIT_FAILURE);
    }

    log_message(INFO, "canu correction done.");

    free(output_correct);
    free(canu_cmd);

    // canu_trim(canu_path, corrected_reads_path, genomeSize, output_dir, readstype, cpu);
    // 
    free(corrected_reads_path);
    
}


void nextdenovo_correct(const char* nextdenovo_path, const char* canu_path, const char* input_seq, 
                        const char* cfg, int cfg_flag, const char* output_dir, char* readstype, char* seqtype, int cpu, int genomeSize) {

    char* output_abs = realpath(output_dir, NULL);
    char* output_correct = (char*)malloc(sizeof(*output_correct) * (snprintf(NULL, 0, "%s/correct_out", output_abs) + 1));
    sprintf(output_correct, "%s/correct_out", output_abs);

    mkdirfiles(output_correct);

    char* cfgw = (char*)malloc(sizeof(*cfgw) * (snprintf(NULL, 0, "%s/correct_out/run.cfg", output_abs) + 1));
    sprintf(cfgw, "%s/correct_out/run.cfg", output_abs);

    char* fofn = (char*)malloc(sizeof(*fofn) * (snprintf(NULL, 0, "%s/correct_out/input.fofn", output_abs) + 1));
    sprintf(fofn, "%s/correct_out/input.fofn", output_abs);
    FILE* fofn_fp = fopen(fofn, "w");
    if (fofn_fp == NULL) {
        log_message(ERROR, "Failed to open fofn file: %s", fofn);
        free(output_correct);
        free(cfgw);
        free(fofn);
        exit(EXIT_FAILURE);
    } else {
        char* inputseq_abs = realpath(input_seq, NULL);
        fprintf(fofn_fp, "%s\n", inputseq_abs);
        fclose(fofn_fp);
        free(inputseq_abs);
    }
    config_info(cfg, cfgw, output_correct, cfg_flag, readstype, fofn, genomeSize, cpu);
    free(fofn);
    
    log_message(INFO, "Start nextdenovo correction...");

    char* nextdenovo_cmd = (char*)malloc(sizeof(*nextdenovo_cmd) * (snprintf(NULL, 0, "%s %s", nextdenovo_path, cfgw) + 1));
    sprintf(nextdenovo_cmd, "%s %s", nextdenovo_path, cfgw);

    execute_command(nextdenovo_cmd, 0, 1);

    // 
    free(nextdenovo_cmd);
    free(cfgw);

    char* corrected_reads_path = (char*)malloc(sizeof(*corrected_reads_path) * (snprintf(NULL, 0, "%s/PMAT.correctedReads.fasta", output_correct) + 1));
    sprintf(corrected_reads_path, "%s/PMAT.correctedReads.fasta", output_correct);

    char* cns_path = (char*)malloc(sizeof(*cns_path) * (snprintf(NULL, 0, "%s/02.cns_align/01.seed_cns.sh.work", output_correct) + 1));
    sprintf(cns_path, "%s/02.cns_align/01.seed_cns.sh.work", output_correct);
    checkfile(cns_path);
    
    cns_files(cns_path, corrected_reads_path);
    free(cns_path);

    log_message(INFO, "nextdenovo correction done.");

    // canu_trim(canu_path, corrected_reads_path, genomeSize, output_abs, seqtype, cpu);
    // 
    free(corrected_reads_path);
    free(output_correct);

}

static void cns_files(const char* cns_dir, const char* output_file) {
    struct dirent* entry;
    DIR* dp = opendir(cns_dir);

    if (dp == NULL) {
        log_message(ERROR, "Failed to open cns directory: %s", cns_dir);
        exit(EXIT_FAILURE);
    }

    FILE* fp = fopen(output_file, "w");
    if (fp == NULL) {
        log_message(ERROR, "Failed to open output file: %s", output_file);
        closedir(dp);
        exit(EXIT_FAILURE);
    }

    while ((entry = readdir(dp)) != NULL) {
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue; // Skip current and parent directories
        }

        // Construct full path
        char full_path[4046];
        snprintf(full_path, sizeof(full_path), "%s/%s", cns_dir, entry->d_name);

        struct stat statbuf;
        if (stat(full_path, &statbuf) == 0 && S_ISDIR(statbuf.st_mode)) {

            char cns_path[4096];
            snprintf(cns_path, sizeof(cns_path), "%s/cns.fasta", full_path);

            FILE* cns_file = fopen(cns_path, "r");
            if (cns_file == NULL) {
                log_message(WARNING, "Failed to open cns file: %s", cns_path);
                continue;
            }

            char* line = NULL;
            size_t line_len = 0;
            while (getline(&line, &line_len, cns_file) != -1) {
                fputs(line, fp);
            }

            free(line);
            fclose(cns_file);
        }
    }

    fclose(fp);
    closedir(dp);
}


static void canu_trim(const char* canu_path, const char* input_seq, int genomeSize, const char* output_dir, const char* readstype, int cpu) {
    log_message(INFO, "Start canu trimming...");

    char* output_trim = (char*)malloc(sizeof(*output_trim) * (snprintf(NULL, 0, "%s/correct_out/trim_out", output_dir) + 1));
    sprintf(output_trim, "%s/correct_out/trim_out", output_dir);

    mkdirfiles(output_trim);

    char* canu_cmd = (char*)malloc(sizeof(*canu_cmd) * (snprintf(NULL, 0, "%s -trim -p PMAT -d %s genomeSize=%d useGrid=false maxThreads=%d -corrected -%s %s", 
                                  canu_path, output_trim, genomeSize, cpu, readstype, input_seq) + 1));
    sprintf(canu_cmd, "%s -trim -p PMAT -d %s genomeSize=%d useGrid=false maxThreads=%d -corrected -%s %s", 
            canu_path, output_trim, genomeSize, cpu, readstype, input_seq);

    execute_command(canu_cmd, 0, 1);

    char* trimmed_reads_path = (char*)malloc(sizeof(*trimmed_reads_path) * (snprintf(NULL, 0, "%s/PMAT.trimmedReads.fasta.gz", output_trim) + 1));
    sprintf(trimmed_reads_path, "%s/PMAT.trimmedReads.fasta.gz", output_trim);

    if (access(trimmed_reads_path, F_OK) != 0) {
        log_message(ERROR, "An error occurred during the trim process?");
        free(output_trim);
        free(canu_cmd);
        free(trimmed_reads_path);
        exit(EXIT_FAILURE);
    }

    log_message(INFO, "canu trimming done.");

    // 
    free(output_trim);
    free(canu_cmd);
    free(trimmed_reads_path);
}



static void config_info(const char* cfgr, const char* cfgw, char* workdir, int cfg_flag, char* readstype, char* input_fofn, int genomeSize, int cpu) {
    FILE* fin = fopen(cfgr, "r");
    if (fin == NULL) {
        log_message(ERROR, "Failed to open config file: %s", cfgr);
        exit(EXIT_FAILURE);
    }

    FILE* fout = fopen(cfgw, "w");
    if (fout == NULL) {
        log_message(ERROR, "Failed to open config file: %s", cfgw);
        fclose(fin);
        exit(EXIT_FAILURE);
    }

    if (cfg_flag == 1) {
        char* line = NULL;
        size_t len = 0;
        while (getline(&line, &len, fin) != -1) {
            if (strncmp(line, "workdir", 7) == 0) {
                fprintf(fout, "workdir = %s\n", workdir);
            } else if (strncmp(line, "task", 4) == 0) {
                fprintf(fout, "task = correct\n");
            } else if (strncmp(line, "input_fofn", 10) == 0) {
                fprintf(fout, "input_fofn = %s\n", input_fofn);
            } else {
                fprintf(fout, "%s", line);
            }
        }
        free(line);
    } else {
        char* line = NULL;
        size_t len = 0;
        while (getline(&line, &len, fin) != -1) {
            if (strncmp(line, "read_type", 9) == 0) {
                fprintf(fout, "read_type = %s\n", readstype);
            } else if (strncmp(line, "parallel_jobs", 13) == 0) {
                fprintf(fout, "parallel_jobs = %d\n", cpu/3);
            } else if (strncmp(line, "correction_options", 18) == 0) {
                fprintf(fout, "correction_options = -p %d\n", cpu/3);
            } else if (strncmp(line, "input_fofn", 10) == 0) {
                fprintf(fout, "input_fofn = %s\n", input_fofn);
            } else if (strncmp(line, "workdir", 7) == 0) {
                fprintf(fout, "workdir = %s\n", workdir);
            } else if (strncmp(line, "genome_size", 11) == 0) {
                fprintf(fout, "genome_size = %d\n", genomeSize);
            } else {
                fprintf(fout, "%s", line);
            }
        }
        free(line);
    }

    fclose(fin);
    fclose(fout);
}
