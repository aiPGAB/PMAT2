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


#ifndef CHECK_FILE_H
#define CHECK_FILE_H

#include <stdint.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

void to_lower(char *str); /* convert all characters to lower case */
void to_upper(char *str); /* convert all characters to upper case */
double findMedian(int arr[], int size); /* find median of array */

int rm_contig(char *str); /* remove contig0 */
void removeUnique(int arr[], int *size); /* remove unique elements in array */
void removeDup(int arr[], int *size); /* remove duplicate elements in array */
void remove_commas(char *str); /* remove commas in string */
void remove_element(int arr[], int *size, int value); /* remove element from array */
int remove_quote(const char *str); /* remove quotes in string */

int findstr(const char* array[], int size, char* target); /* return index of target in array, -1 if not found */
int findint(const int array[], int array_size, int target); /* return index of target in array, -1 if not found */
char* abspath(const char *path);                /* return absolute path */
char* pmat_path(const char *input);             /* return path of PMAT2 */
int is_numeric(const char *str);                /* check if string is numeric */ 

void sleep_ms(long milliseconds);

void mkdirfiles(const char *dir_path);          /* create directory and all intermediate directories if not exist */
int delete_directory(const char *path);         /* delete directory recursively */
void checkfile(const char *dir_path);           /* check if directory exists, create directory if not exist */
int is_file(const char *path);
void remove_file(const char *filename);        /* remove file */
void rename_file(const char *old_name, const char *new_name);
int is_gzipped_file(const char* filename);      /* check if file is gzipped */

long getFileSize(const char *filename);         /* get file size */
void execute_command(const char* command, int verbose, int log_output); // execute command and print output if verbose is true
int ass_command(const char* command, int verbose, int log_output);
int which_executable(const char *exe);          /* return 1 if executable exists, 0 otherwise */
int validate_fastq_file(const char *filename);  /* check if file is a valid fastq file */
int validate_fasta_file(const char *filename);  /* check if file is a valid fasta file */
int is_digits(const char *str);                 /* check if string contains only digits */

void maparr_100(int32_t *data, int size, uint8_t *mapped_data); /* normalize data to 0-100 scale */

#endif
