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
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include "log.h"
#include "version.h"

#define MAX_LOG_MESSAGE_LENGTH 4096

FILE* log_output_stream = NULL;
FILE* log_file_stream = NULL;

static int atexit_registered = 0;

static const char *log_level_strings[] = {
    "INFO",
    "WARNING",
    "ERROR"
};

void set_log_output(FILE* stream) {
    log_output_stream = stream;
}

int init_log_file(const char* filepath) {
    if (log_file_stream != NULL) {
        close_log_file();
    }
    if (filepath == NULL || filepath[0] == '\0') {
        return -1;
    }
    log_file_stream = fopen(filepath, "w");
    if (log_file_stream == NULL) {
        fprintf(stderr, "[WARNING] Failed to open log file for writing: %s\n", filepath);
        return -1;
    }
    if (!atexit_registered) {
        atexit(close_log_file);
        atexit_registered = 1;
    }
    return 0;
}

void close_log_file(void) {
    if (log_file_stream != NULL) {
        time_t rawtime;
        struct tm timeinfo;
        char time_buffer[32];
        time(&rawtime);
        localtime_r(&rawtime, &timeinfo);
        strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", &timeinfo);
        fprintf(log_file_stream, "\n================================================================================\n");
        fprintf(log_file_stream, "End Time   : %s\n", time_buffer);
        fprintf(log_file_stream, "================================================================================\n");
        fflush(log_file_stream);
        fclose(log_file_stream);
        log_file_stream = NULL;
    }
}

void log_header(int argc, char* argv[]) {
    if (log_file_stream == NULL) {
        return;
    }
    time_t rawtime;
    struct tm timeinfo;
    char time_buffer[32];
    time(&rawtime);
    localtime_r(&rawtime, &timeinfo);
    strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", &timeinfo);

    fprintf(log_file_stream, "================================================================================\n");
    fprintf(log_file_stream, "PMAT v%s Run Log\n", VERSION_PMAT);
    fprintf(log_file_stream, "Start Time : %s\n", time_buffer);
    if (argc > 0 && argv != NULL) {
        fprintf(log_file_stream, "Command    : ");
        for (int i = 0; i < argc; i++) {
            fprintf(log_file_stream, "%s%s", argv[i], (i + 1 < argc) ? " " : "");
        }
        fprintf(log_file_stream, "\n");
    }
    fprintf(log_file_stream, "================================================================================\n\n");
    fflush(log_file_stream);
}

void log_section_header(const char* message) {
    printf("** %s \n", message);
    fflush(stdout);
    if (log_file_stream != NULL) {
        fprintf(log_file_stream, "** %s \n", message);
        fflush(log_file_stream);
    }
}

void log_section_tail(const char* message) {
    printf("** %s \n", message);
    fflush(stdout);
    if (log_file_stream != NULL) {
        fprintf(log_file_stream, "** %s \n", message);
        fflush(log_file_stream);
    }
}

void log_info(const char* format, ...) {
    FILE* console = (log_output_stream != NULL) ? log_output_stream : stdout;

    va_list args;
    va_start(args, format);
    vfprintf(console, format, args);
    va_end(args);
    fflush(console);

    if (log_file_stream != NULL) {
        va_list args_file;
        va_start(args_file, format);
        vfprintf(log_file_stream, format, args_file);
        va_end(args_file);
        fflush(log_file_stream);
    }
}

void log_message(int level, const char *fmt, ...) {
    time_t rawtime;
    struct tm timeinfo;
    char time_buffer[32];
    char message[MAX_LOG_MESSAGE_LENGTH];

    if (level < INFO || level > ERROR) {
        fprintf(stderr, "Invalid log level\n");
        return;
    }

    time(&rawtime);
    localtime_r(&rawtime, &timeinfo);
    strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", &timeinfo);

    va_list args;
    va_start(args, fmt);
    vsnprintf(message, sizeof(message), fmt, args);
    va_end(args);

    FILE* console = (log_output_stream != NULL) ? log_output_stream : ((level == INFO) ? stdout : stderr);
    fprintf(console, "[%s] %s: %s\n", time_buffer, log_level_strings[level], message);
    fflush(console);

    if (log_file_stream != NULL) {
        fprintf(log_file_stream, "[%s] %s: %s\n", time_buffer, log_level_strings[level], message);
        fflush(log_file_stream);
    }
}
