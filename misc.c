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
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>   // gzFile, gzgets, gzopen, gzclose
#include <ctype.h>  // tolower
#include <errno.h> // errno
#include <sys/wait.h> // waitpid
#include <fcntl.h> // fcntl, F_GETFL, F_SETFL, O_NONBLOCK
#include <unistd.h> // access
#include <sys/stat.h> // mkdir
#include <dirent.h> // opendir, readdir, closedir, rewinddir, struct dirent
#include <time.h>
#include <signal.h>


#include "log.h"
#include "misc.h"


/* find a string in an array of strings */
int findstr(const char* array[], int array_size, char* target) {
    uint64_t i;
    for (i = 0; i < array_size; i++) {
        if (strcmp(array[i], target) == 0) {
            return 1;
        }
    }
    return 0;
}


int findint(const int array[], int array_size, int target) {
    uint64_t i;
    for (i = 0; i < array_size; i++) {
        if (array[i] == target) {
            return 1;
        }
    }
    return 0;
}

/* check if a string is a number */
int is_numeric(const char *str) {
    uint64_t i;
    if (str == NULL || *str == '\0') return 0;
    for (i = 0; i < strlen(str); i++) {
        if (!isdigit((unsigned char)str[i])) {
            return 0;
        }
    }
    return 1;
}


/* get the absolute path of a file */
char* abspath(const char *path) {
    char *absolute_path = realpath(path, NULL);
    return absolute_path;
}

/* get the software path from an environment variable or a path */
char* pmat_path(const char *prog_name) {
    char *resolved_path = NULL;

    if (strchr(prog_name, '/')) {
        resolved_path = realpath(prog_name, NULL);
        if (!resolved_path) {
            log_message(ERROR, "Error resolving path '%s': %s", prog_name, strerror(errno));
            exit(EXIT_FAILURE);
        }
        return resolved_path;
    }

    const char *path_env = getenv("PATH");
    if (!path_env) {
        log_message(ERROR, "PATH environment variable not set");
        return NULL;
    }

    char *paths = strdup(path_env);
    char *dir = strtok(paths, ":");

    while (dir) {
        char pmat_path[4098];
        snprintf(pmat_path, sizeof(pmat_path), "%s/%s", dir, prog_name);
        char pmat_c[4098];
        snprintf(pmat_c, sizeof(pmat_c), "%s/PMAT.c", dir);

        if (access(pmat_path, X_OK) == 0 && access(pmat_c, F_OK) == 0) {
            resolved_path = realpath(pmat_path, NULL);
            break;
        }

        dir = strtok(NULL, ":");
    }

    free(paths);
    return resolved_path;
}

/* convert a string to lower case */
void to_lower(char *str) {
    uint64_t i;
    for (i = 0; str[i]; i++) {
        str[i] = tolower(str[i]);
    }
}

void to_upper(char *str) {
    uint64_t i;
    for (i = 0; str[i]; i++) {
        str[i] = toupper(str[i]);
    }
}

/* remove the "contig0" */
int rm_contig(char *str) {
    char *p = str;

    if (strncmp(str, "contig", 6) == 0) {
        p += 6;
    }

    while (*p == '0') {
        p++;
    }

    memmove(str, p, strlen(p) + 1);
    return atoi(str);
}

void sleep_ms(long milliseconds) {
    struct timespec req = {0};
    req.tv_sec = milliseconds / 1000;
    req.tv_nsec = (milliseconds % 1000) * 1000000;
    nanosleep(&req, NULL);
}

void removeUnique(int arr[], int *size) {
    if (*size <= 1) return;
    uint64_t i, j;
    bool *isDuplicate = (bool *)calloc(*size, sizeof(bool));
    int newSize = 0;

    for (i = 0; i < *size - 1; i++) {
        if (isDuplicate[i]) continue;
        for (j = i + 1; j < *size; j++) {
            if (arr[i] == arr[j]) {
                isDuplicate[i] = true;
                isDuplicate[j] = true;
                break;
            }
        }
    }

    for (i = 0; i < *size; i++) {
        if (isDuplicate[i]) {
            arr[newSize++] = arr[i];
        }
    }

    *size = newSize;
    free(isDuplicate);
}


void removeDup(int arr[], int *size) {
    if (*size <= 1) return;

    int newSize = 1;
    uint64_t i, j;
    for (i = 1; i < *size; i++) {
        bool isDuplicate = false;
        
        for (j = 0; j < newSize; j++) {
            if (arr[i] == arr[j]) {
                isDuplicate = true;
                break;
            }
        }
        
        if (!isDuplicate) {
            arr[newSize] = arr[i];
            newSize++;
        }
    }

    *size = newSize;
}


void remove_commas(char *str) {
    char *src, *dst;
    for (src = dst = str; *src != '\0'; src++) {
        if (*src != ',') {
            *dst++ = *src;
        }
    }
    *dst = '\0';
}


int remove_quote(const char *str) {
    int len = strlen(str);
    
    char *temp_str = (char *)malloc(len + 1);
    if (!temp_str) {
        fprintf(stderr, "Memory allocation failed\n");
        return -1;
    }
    strcpy(temp_str, str);
    
    if (temp_str[len - 1] == '\'') {
        temp_str[len - 1] = '\0';
    }
    
    int result = atoi(temp_str);
    
    free(temp_str);
    
    return result;
}

void remove_element(int arr[], int *size, int value) {
    int new_size = 0;
    uint64_t i;
    for (i = 0; i < *size; i++) {
        if (arr[i] != value) {
            arr[new_size++] = arr[i];
        }
    }

    *size = new_size;
}


long getFileSize(const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        log_message(ERROR, "Failed to open file: %s", filename);
        exit(EXIT_FAILURE);
    }

    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    fclose(file);

    return size;
}

void check_executable(const char *exe) {
    if (access(exe, X_OK) == 0) {
        return;
    } else {
        if (strchr(exe, '/') == NULL) {
            char command[1024];
            snprintf(command, sizeof(command), "which %s > /dev/null 2>&1", exe);
            
            if (system(command) == 0) {
                return;
            } else {
                log_message(ERROR, "Executable not found: %s", exe);
                exit(EXIT_FAILURE);
            }
        }
        
        log_message(ERROR, "Executable not found: %s", exe);
        exit(EXIT_FAILURE);
    }
}

int is_digits(const char *str) {
    if (str == NULL || *str == '\0') {
        return 0;
    }
    uint64_t i;
    for (i = 0; i < strlen(str); i++) {
        if (!isdigit((unsigned char)str[i])) {
            return 0;
        }
    }
    return 1;
}

int which_executable(const char *exe) {
    char command[4096];
    snprintf(command, sizeof(command), "command -v %s > /dev/null 2>&1", exe);
    
    if (system(command) == 0) {
        return 1;
    }
    
    return 0;
}

/* create and check files and directories */
void mkdirfiles(const char *dir_path) {
    if (access(dir_path, F_OK) != 0) {
        if (mkdir(dir_path, 0777) != 0) {
            log_message(ERROR, "Error creating output directory");
            exit(EXIT_FAILURE);
        }
    }
}

void rename_file(const char *old_name, const char *new_name) {
    if (rename(old_name, new_name) != 0) {
        log_message(ERROR, "Error renaming file: %s", old_name);
    }
}

void remove_file(const char *filename) {
    if (remove(filename) != 0) {
        if (errno != ENOENT) {
            log_message(ERROR, "Error removing file: %s", filename);
            exit(EXIT_FAILURE); 
        }
    }
}

void checkfile(const char *dir_path) {
    if (access(dir_path, F_OK) != 0) {
        log_message(ERROR, "File not found: %s", dir_path);
        exit(EXIT_FAILURE);
    }
}

int is_file(const char *path) {
    struct stat statbuf;
    if (stat(path, &statbuf) == -1) {
        return 0;
    }

    if (S_ISREG(statbuf.st_mode)) {
        return 1;
    } else {
        return 0;
    }
}


int delete_directory(const char *path) {
    struct dirent *entry;
    DIR *dir = opendir(path);

    if (dir == NULL) {
        return -1;
    }

    while ((entry = readdir(dir)) != NULL) {
        char full_path[4096];
        struct stat statbuf;

        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }

        snprintf(full_path, sizeof(full_path), "%s/%s", path, entry->d_name);

        if (stat(full_path, &statbuf) == -1) {
            closedir(dir);
            return -1;
        }

        if (S_ISDIR(statbuf.st_mode)) {
            if (delete_directory(full_path) == -1) {
                closedir(dir);
                return -1;
            }
        } else {
            if (remove(full_path) == -1) {
                closedir(dir);
                return -1;
            }
        }
    }

    closedir(dir);

    if (rmdir(path) == -1) {
        return -1;
    }

    return 0;
}


int is_gzipped_file(const char* filename) {
    const char* ext = strrchr(filename, '.');
    return ext && strcmp(ext, ".gz") == 0;
}


int compare(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

double findMedian(int arr[], int size) {
    qsort(arr, size, sizeof(int), compare);

    if (size % 2 == 0) {
        return (arr[size / 2 - 1] + arr[size / 2]) / 2.0;
    } else {
        return arr[size / 2];
    }
}

/* execute a command and log its output */
void execute_command(const char* command, int verbose, int log_output) {
    if (verbose) {
        log_info("Running command: %s\n", command);
    }

    int pipefd[2];
    if (pipe(pipefd) == -1) {
        log_message(ERROR, "Failed to create pipe: %s", strerror(errno));
        exit(EXIT_FAILURE);
    }

    pid_t pid = fork();
    if (pid == -1) {
        log_message(ERROR, "Failed to fork: %s", strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (pid == 0) {  // Child process
        close(pipefd[0]);  // Close unused read end
        dup2(pipefd[1], STDOUT_FILENO);
        dup2(pipefd[1], STDERR_FILENO);
        close(pipefd[1]);

        signal(SIGHUP, SIG_IGN);
        execl("/bin/sh", "sh", "-c", command, (char *)NULL);

        log_message(ERROR, "Failed to execute command: %s", strerror(errno));
        exit(EXIT_FAILURE);
    } else {  // Parent process
        close(pipefd[1]);  // Close unused write end

        // Set the pipe to non-blocking mode
        int flags = fcntl(pipefd[0], F_GETFL, 0);
        fcntl(pipefd[0], F_SETFL, flags | O_NONBLOCK);

        char buffer[4096];
        ssize_t bytes_read;
        char line_buffer[4096] = {0};
        size_t line_buffer_pos = 0;
        size_t i;
        
        char filtered_output[65536] = {0};  // To hold only filtered, meaningful lines
        size_t filtered_pos = 0;

        while (1) {
            bytes_read = read(pipefd[0], buffer, sizeof(buffer) - 1);
            if (bytes_read > 0) {
                buffer[bytes_read] = '\0';
                for (i = 0; i < bytes_read; i++) {
                    if (buffer[i] == '\n' || line_buffer_pos == sizeof(line_buffer) - 1) {
                        line_buffer[line_buffer_pos] = '\0';

                        // Filter unwanted log lines
                        if (strstr(line_buffer, "Warning: [blastn]") == NULL &&
                            strstr(line_buffer, "Examining 5 or more matches is recommended") == NULL &&
                            strstr(line_buffer, "GenomeScope analyzing") == NULL &&
                            strstr(line_buffer, "Model converged") == NULL) {

                            // Remove carriage returns
                            char* cr = strchr(line_buffer, '\r');
                            while (cr != NULL) {
                                memmove(cr, cr + 1, strlen(cr + 1) + 1);
                                cr = strchr(cr, '\r');
                            }

                            if (log_output) {
                                log_info("%s\n", line_buffer);
                            }

                            // Save filtered line for potential error report
                            if (filtered_pos + line_buffer_pos + 2 < sizeof(filtered_output)) {
                                memcpy(filtered_output + filtered_pos, line_buffer, line_buffer_pos);
                                filtered_pos += line_buffer_pos;
                                filtered_output[filtered_pos++] = '\n';
                                filtered_output[filtered_pos] = '\0';
                            }
                        }

                        line_buffer_pos = 0;
                    } else {
                        line_buffer[line_buffer_pos++] = buffer[i];
                    }
                }
            } else if (bytes_read == 0) {
                // EOF
                break;
            } else {
                sleep_ms(10);  // Brief pause
            }
        }

        close(pipefd[0]);

        int status;
        waitpid(pid, &status, 0);
        if (WIFEXITED(status)) {
            int exit_status = WEXITSTATUS(status);
            if (exit_status != 0) {
                log_message(ERROR, "Command failed with status: %d", exit_status);
                if (filtered_pos > 0) {
                    log_message(ERROR, "%s", filtered_output);
                }
                exit(EXIT_FAILURE);
            }
        } else {
            log_message(ERROR, "Command did not exit normally");
            exit(EXIT_FAILURE);
        }
    }
}



int ass_command(const char *command, int verbose, int log_output) {
    if (verbose) {
        log_info("Running command: %s\n", command);
    }

    FILE *fp = popen(command, "r");
    if (fp == NULL) {
        log_message(ERROR, "Failed to open pipe: %s", strerror(errno));
        exit(EXIT_FAILURE);
    }

    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), fp) != NULL) {
        buffer[strcspn(buffer, "\n")] = '\0';

        if (strstr(buffer, "v3.0 (20140410_1040)") != NULL ||
            strstr(buffer, "Warning:  No quality scores file found.") != NULL) {
            continue;
        }

        if (strstr(buffer, "chord->getLength()") != NULL) {
            pclose(fp);
            return -2;
        }

        if (strstr(buffer, "doAsmAlignment traceback") != NULL) {
            log_message(ERROR, "Command did not exit normally");
            pclose(fp);
            return -1;
        }

        if (log_output == 0 && strstr(buffer, "Indexing PMAT_cut_seq.fa...") != NULL) {
            log_output = 1;
        }

        if (log_output) {
            log_info("%s\n", buffer);
        }
        fflush(stdout);
    }

    int status = pclose(fp);
    if (WIFEXITED(status)) {
        int exit_status = WEXITSTATUS(status);
        if (exit_status != 0) {
            log_message(ERROR, "Command failed with status: %d", exit_status);
            exit(EXIT_FAILURE);
        }
    } else {
        log_message(ERROR, "Command did not exit normally");
        exit(EXIT_FAILURE);
    }

    return 0;
}



/* validate a FASTA file */
int validate_fasta_file(const char* filename) {
    if (filename == NULL) {
        return 0;
    }

    // gzipped fasta file
    FILE* file = fopen(filename, "rb");
    if (file == NULL) {
        return 0;
    }

    unsigned char buf[2];
    size_t read_bytes = fread(buf, 1, 2, file);
    fclose(file);

    if (read_bytes == 2 && buf[0] == 0x1F && buf[1] == 0x8B) {
        // This is a gzipped file, let's check if it contains a FASTA format
        gzFile gzfile = gzopen(filename, "rb");
        if (gzfile == NULL) {
            return 0;
        }

        char first_line[1024];
        if (gzgets(gzfile, first_line, sizeof(first_line)) == NULL) {
            gzclose(gzfile);
            return 0;
        }

        gzclose(gzfile);

        // Check if the first line starts with '@'
        if (first_line[0] == '>') {
            return 1;  // gzipped FASTA
        }

        return 0;  // gzipped, but not FASTA
    }

    // regular fasta file
    file = fopen(filename, "r");
    if (file == NULL) {
        return 0;
    }

    char first_line[1024];
    if (fgets(first_line, sizeof(first_line), file) == NULL) {
        fclose(file);
        return 0;
    }

    fclose(file);

    if (first_line[0] == '>') {
        return 2;  // plain FASTA
    }

    return 0;  // not a FASTA file
}


/* validate a FASTQ file */
int validate_fastq_file(const char* filename) {
    if (filename == NULL) {
        return 0;
    }

    // gzipped fastq file
    FILE* file = fopen(filename, "rb");
    if (file == NULL) {
        return 0;
    }

    unsigned char buf[2];
    size_t read_bytes = fread(buf, 1, 2, file);
    fclose(file);

    if (read_bytes == 2 && buf[0] == 0x1F && buf[1] == 0x8B) {
        // This is a gzipped file, let's check if it contains a FASTQ format
        gzFile gzfile = gzopen(filename, "rb");
        if (gzfile == NULL) {
            return 0;
        }

        char first_line[1024];
        if (gzgets(gzfile, first_line, sizeof(first_line)) == NULL) {
            gzclose(gzfile);
            return 0;
        }

        gzclose(gzfile);

        // Check if the first line starts with '@'
        if (first_line[0] == '@') {
            return 1;  // gzipped FASTQ
        }

        return 0;  // gzipped, but not FASTQ
    }

    // regular fastq file
    file = fopen(filename, "r");
    if (file == NULL) {
        return 0;
    }

    char first_line[1024];
    if (fgets(first_line, sizeof(first_line), file) == NULL) {
        fclose(file);
        return 0;
    }

    fclose(file);

    if (first_line[0] == '@') {
        return 2;  // plain FASTQ
    }

    return 0;  // not a FASTQ file
}


void maparr_100(int32_t *data, int size, uint8_t *mapped_data) {
    int32_t min = data[0], max = data[0];
    uint64_t i;
    for (i = 1; i < size; i++) {
        if (data[i] < min) min = data[i];
        if (data[i] > max) max = data[i];
    }
    
    if (max == min) {
        for (i = 0; i < size; i++) {
            mapped_data[i] = 50;
        }
    } else {
        for (i = 0; i < size; i++) {
            double mapped = (data[i] - min) * 99.0 / (max - min) + 1;
            mapped_data[i] = (uint8_t)(mapped < 1 ? 1 : (mapped > 100 ? 100 : mapped));
        }
    }
}