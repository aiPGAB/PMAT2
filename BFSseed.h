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


#ifndef BFSSEED_H
#define BFSSEED_H

#include "hitseeds.h"

typedef struct {
    int lctg;
    int rctg;
    char *lutr;
    char *rutr;
    float linkdepth;
} Ctglinks;


// typedef struct {
//     char *lctg;
//     int lctgsmp;
//     float lctgdepth;
//     int lctglen;
//     char *rctg;
//     int rctgsmp;
//     float rctgdepth;
//     int rctglen;
//     char *lutr;
//     char *rutr;
//     float linkdepth;
// } BFSlinks;
typedef struct {
    int lctgsmp;
    float lctgdepth;
    int lctglen;
    int rctgsmp;
    float rctgdepth;
    int rctglen;
    int lutrsmp;
    int rutrsmp;
    float linkdepth;
} BFSlinks;

void BFSseeds(const char* type, int num_links, int num_ctg, 
            Ctglinks* ctglinks, CtgDepth* ctgdepth, 
            int* num_dynseeds, int** dynseeds, float nucl_depth, 
            float filter_depth, BFSlinks** bfslinks, int* num_bfslinks);


#endif