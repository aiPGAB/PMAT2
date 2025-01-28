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


#include <stdarg.h>
#include <time.h>
#include "log.h"

#define MAX_LOG_MESSAGE_LENGTH 4096

FILE* log_output_stream = NULL;

static const char *log_level_strings[] = {
    "INFO",
    "WARNING",
    "ERROR"
};

void set_log_output(FILE* stream) {
    log_output_stream = stream;
}

void log_section_header(const char* message) {
    printf("** %s \n", message);
}

void log_section_tail(const char* message) {
    printf("** %s \n", message);
}

void log_info(const char* format, ...) {
    va_list args;
    va_start(args, format);
    if (log_output_stream == NULL) {
        log_output_stream = stdout;
    }
    vfprintf(log_output_stream, format, args);
    va_end(args);
    fflush(log_output_stream);
}

void log_message(int level, const char *fmt, ...) {
    va_list args;
    time_t rawtime;
    struct tm timeinfo;
    char time_buffer[32];
    char message[MAX_LOG_MESSAGE_LENGTH];

    if (log_output_stream == NULL) {
        log_output_stream = (level == INFO) ? stdout : stderr;
    }

    if (level < INFO || level > ERROR) {
        fprintf(stderr, "Invalid log level\n");
        return;
    }

    time(&rawtime);
    localtime_r(&rawtime, &timeinfo);
    strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", &timeinfo);

    va_start(args, fmt);
    vsnprintf(message, sizeof(message), fmt, args);
    va_end(args);

    fprintf(log_output_stream, "[%s] %s: %s\n", time_buffer, log_level_strings[level], message);
    fflush(log_output_stream);
}
