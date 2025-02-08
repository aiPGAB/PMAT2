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
#include <getopt.h>
#include <libgen.h> // dirname
#include <signal.h>


#include "log.h"
#include "version.h"
#include "misc.h"
#include "pmat.h"


void usage() {
    fprintf(stdout, 
        "usage: PMAT <command> <arguments>\n"
        "\n"
        "  ______     ___           __        ____       _____________     \n"
        " |   __  \\  |   \\        /   |      / __ \\     |_____   _____| \n"
        " |  |__)  | | |\\ \\      / /| |     / /  \\ \\          | |      \n"
        " |   ____/  | | \\ \\    / / | |    / /____\\ \\         | |      \n"
        " |  |       | |  \\ \\  / /  | |   / /______\\ \\        | |      \n"
        " |  |       | |   \\ \\/ /   | |  / /        \\ \\       | |      \n"
        " |__|       |_|    \\__/    |_| /_/          \\_\\      |_|       \n"
        "\n"
        "PMAT2           an efficient assembly toolkit for organellar genome\n"
        "Contributors    Bi,C. and Han,F.\n"
        "Email           bichwei@njfu.edu.cn, hanfc@caf.ac.cn\n"
        "Version         PMAT v%s\n"
        "\n"
        "For more information about PMAT, see https://github.com/aiPGAB/PMAT2\n"
        "\n"
        "Commands:\n"
        "\n"
        "    autoMito    One-step de novo assembly of organellar genomes.    \n"
        "                This command processes raw ONT/CLR data or uses     \n"
        "                corrected data or HiFi reads for direct assembly.   \n"
        "                Based on the assembly result, it automatically      \n"
        "                selects seeds for extension and filters false       \n"
        "                positives to obtain the complete organellar         \n"
        "                genome sequence.                                    \n"
        "\n"
        "    graphBuild  If the autoMito command fails to generate the   \n"
        "                complete organellar genome sequence in one-step \n"
        "                assembly, you can use this command to manually  \n"
        "                select seeds for assembly.                      \n"
        "\n"
        "Optional options:\n"
        "   -v, --version   show program's version number and exit\n"
        "   -h, --help      show this help message and exit\n"
        , VERSION_PMAT
    );
}


void autoMito_usage() {
    fprintf(stdout,
        "Usage: PMAT autoMito [-i INPUT] [-o OUTPUT] [-t SEQTYPE] [options]\n"
        "Example:\n"
        "       PMAT autoMito -i hifi.fastq.gz -o hifi_assembly -t hifi -m -T 8\n"
        "       PMAT autoMito -i ont.fastq.gz -o ont_assembly -t ont -S nextdenovo -C canu -N nextdenovo\n"
        "       PMAT autoMito -i clr.fastq.gz -o clr_assembly -t clr -S canu -C canu\n\n"
        

        "Required options:\n"
        "   -i, --input          Input sequence file (fasta/fastq)\n"
        "   -o, --output         Output directory\n"
        "   -t, --seqtype        Sequence type (hifi/ont/clr)\n"
        "\n"

        "Optional options:\n"
        "   -k, --kmer           kmer size for estimating genome size (default: 31)\n"
        "   -g, --genomesize     Genome size (g/m/k), skip genome size estimation if set\n"
        "   -p, --task           Task type (0/1), skip error correction for ONT/CLR by selecting 0, otherwise 1 (default: 1)\n"
        "   -G, --organelles     Genome organelles (mt/pt/all, default: mt)\n"
        "   -x, --taxo           Specify the organism type (0/1/2), 0: plants, 1: animals, 2: Fungi (default: 0)\n"
        "   -S, --correctsoft    Error correction software (canu/nextdenovo, default: nextdenovo)\n"
        "   -C, --canu           Canu path\n"
        "   -N, --nextdenovo     NextDenovo path\n"
        "   -n, --cfg            Config file for nextdenovo (default: temprun.cfg)\n"
        "   -F, --factor         Subsample factor (default: 1)\n"
        "   -D, --subseed        Random number seeding when extracting subsets (default: 6)\n"
        "   -K, --breaknum       Break long reads (>30k) with this (default: 20000)\n"
        "   -I, --minidentity    Set minimum overlap identity (default: 90)\n"
        "   -L, --minoverlaplen  Set minimum overlap length (default: 40)\n"
        "   -T, --cpu            Number of threads (default: 8)\n"
        "   -m, --mem            Keep sequence data in memory to speed up computation\n"
        "   -h, --help           Show this help message and exit\n"
        
    );
}


void graphBuild_usage() {
    fprintf(stdout,
        "Usage: PMAT graphBuild [-i SUBSAMPLE] [-a ASSEMBLY] [-o OUTPUT] [options]\n"
        "Example:\n"
        "       PMAT graphBuild -i assembly_test1/subsample -a assembly_test1/assembly_result -o graphBuild_result -s 1 312 356 -T 8\n"
        "       PMAT graphBuild -i assembly_test1/subsample -a assembly_test1/assembly_result -o graphBuild_result -d 5 -s 1 312 356 -T 8\n\n"
        

        "Required options:\n"
        "   -i, --subsample     Input subsample directory (assembly_test1/subsample)\n"
        "   -a, --graphinfo     Input assembly result directory (assembly_test1/assembly_result)\n"
        "   -o, --output        Output directory\n"
        "\n"

        "Optional options:\n"
        "   -G, --organelles    Genome organelles (mt/pt, default: mt)\n"
        "   -x, --taxo          Specify the organism type (0/1), 0: plants, 1: animals, 2: Fungi (default: 0)\n"
        "   -d, --depth         Contig depth threshold\n"
        "   -s, --seeds         ContigID for extending. Multiple contigIDs should be separated by space. For example: 1 312 356\n"
        "   -T, --cpu           Number of threads (default: 8)\n"
        "   -h, --help          Show this help message and exit\n"
    );
}



void autoMito_arguments(int argc, char *argv[], char* exe_path, autoMitoArgs *opts) {
    
    int option_index = 0;
    static struct option long_options[] = {
        {"input", 1, 0, 'i'},
        {"output", 1, 0, 'o'},
        {"seqtype", 1, 0, 't'},
        {"kmer", 1, 0, 'k'},
        {"genomesize", 1, 0, 'g'},
        {"task", 1, 0, 'p'},
        {"taxo", 1, 0, 'x'},
        {"organelles", 1, 0, 'G'},
        {"correctsoft", 1, 0, 'S'},
        {"canu", 1, 0, 'C'},
        {"nextdenovo", 1, 0, 'N'},
        {"cfg", 1, 0, 'n'},
        {"factor", 1, 0, 'F'},
        {"subseed", 1, 0, 'D'},
        {"breaknum", 1, 0, 'K'},
        {"minidentity", 1, 0, 'I'},
        {"minoverlaplen", 1, 0, 'L'},
        {"cpu", 1, 0, 'T'},
        {"mem", 0, 0, 'm'},
        {"help", 0, 0, 'h'},
        {"version", 0, 0, 'v'},
        {0, 0, 0, 0}
    };

    while (1) {
        int c = getopt_long(argc, argv, ":i:o:t:k:g:p:G:x:S:C:N:n:F:D:K:I:L:T:mhv", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 'i': opts->input_file = optarg; break;
            case 'o': opts->output_file = optarg; break;
            case 't': opts->seqtype = optarg; break;
            case 'k': opts->kmersize = atoi(optarg); break;
            case 'g': 
                if (optind < argc && argv[optind][0] != '-') {
                    log_message(ERROR, "Missing value for genomesize");
                    autoMito_usage();
                    exit(EXIT_FAILURE);
                }
                opts->genomesize = optarg;
                break;
            case 'p': opts->task = atoi(optarg); break;
            case 'G': opts->organelles = optarg; break;
            case 'x': opts->taxo = atoi(optarg); break;
            case 'S': opts->correct_software = optarg; break;
            case 'C': opts->canu_path = optarg; break;
            case 'N': opts->nextdenovo_path = optarg; break;
            case 'n': opts->cfg_file = optarg; opts->cfg_flag = 1; break;
            case 'F': opts->factor = atof(optarg); break;
            case 'D': opts->seed = atoi(optarg); break;
            case 'K': opts->breaknum = atoi(optarg); break;
            case 'I': opts->mi = atoi(optarg); break;
            case 'L': opts->ml = atoi(optarg); break;
            case 'T': opts->cpu = atoi(optarg); break;
            case'm': opts->mem = 1; break;
            case 'h': autoMito_usage(); exit(EXIT_SUCCESS);
            case 'v': log_info("PMAT v%s\n", VERSION_PMAT); exit(EXIT_SUCCESS);
            default: log_message(ERROR, "Invalid option: %c", c); autoMito_usage(); exit(EXIT_FAILURE);
        }
    }
    if (opts->input_file == NULL || opts->output_file == NULL || opts->seqtype == NULL) {
        log_message(ERROR, "Missing required options");
        autoMito_usage();
        exit(EXIT_FAILURE);
    } 
    
    if (is_file(opts->input_file) == 0) {
        log_message(ERROR, "Input file does not exist: %s", opts->input_file);
        autoMito_usage();
        exit(EXIT_FAILURE);
    }

    if (opts->taxo != 0 && opts->taxo != 1 && opts->taxo != 2) {
        log_message(ERROR, "Invalid taxo type: %d", opts->taxo);
        exit(EXIT_FAILURE);
    }

    if (opts->task != 0 && opts->task != 1) {
        log_message(ERROR, "Invalid task type: %d", opts->task);
        exit(EXIT_FAILURE);
    }
    if (opts->factor < 0 || opts->factor > 1) {
        log_message(ERROR, "Invalid factor: %f", opts->factor);
        exit(EXIT_FAILURE);
    }
    if (opts->seed < 0) {
        log_message(ERROR, "Invalid subseed: %d", opts->seed);
        exit(EXIT_FAILURE);
    }
    if (opts->breaknum < 100) {
        log_message(ERROR, "Invalid breaknum: %d", opts->breaknum);
        exit(EXIT_FAILURE);
    }
    if (opts->cpu < 1) {
        log_message(ERROR, "Invalid cpu: %d", opts->cpu);
        exit(EXIT_FAILURE);
    }

    if (opts->kmersize < 1 || opts->kmersize > 31) {
        log_message(ERROR, "Invalid kmer size (k<=31): %d", opts->kmersize);
        exit(EXIT_FAILURE);
    }

    /* get sequence type */
    to_lower(opts->seqtype);
    if (strcmp(opts->seqtype, "hifi") != 0 && strcmp(opts->seqtype, "clr") != 0 && strcmp(opts->seqtype, "ont") != 0) {
        log_message(ERROR, "Invalid sequence type (hifi/ont/clr): %s", opts->seqtype);
        exit(EXIT_FAILURE);
    }

    if ((strcmp(opts->seqtype, "clr") == 0 || strcmp(opts->seqtype, "ont") == 0) && opts->task == 1) {
        if (opts->correct_software != NULL) {
            to_lower(opts->correct_software);
        } else {
            opts->correct_software = strdup("nextdenovo");
        }

        if (strcmp(opts->correct_software, "canu") == 0) {
            if (opts->canu_path == NULL) {
                if (which_executable("canu") == 0) {
                    log_message(ERROR, "Can't find Canu, please specify the path with -C");
                    exit(EXIT_FAILURE);
                } else {
                    opts->canu_path = strdup("canu");
                }
            }
            if (which_executable(opts->canu_path) == 0) {
                log_message(ERROR, "Can't find Canu, please specify the correct path with -C");
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(opts->correct_software, "nextdenovo") == 0) {
            if (opts->nextdenovo_path == NULL) {
                if (which_executable("nextDenovo") == 0) {
                    log_message(ERROR, "Can't find NextDenovo, please specify the path with -N");
                    exit(EXIT_FAILURE);
                } else {
                    opts->nextdenovo_path = strdup("nextDenovo");
                }
            }
            if (which_executable(opts->nextdenovo_path) == 0) {
                log_message(ERROR, "Can't find NextDenovo, please specify the correct path with -N");
                exit(EXIT_FAILURE);
            }

        } else {
            log_message(ERROR, "Invalid error correction software (canu/nextdenovo) : %s", opts->correct_software);
            exit(EXIT_FAILURE);
        }

        if (opts->cfg_flag == 1) {
            checkfile(opts->cfg_file);
        } else {
            char *dir = dirname(strdup(exe_path));
            char *cfg_path = malloc(strlen(dir) + strlen("/temprun.cfg") + 1);
            strcpy(cfg_path, dir);
            strcat(cfg_path, "/temprun.cfg");
            opts->cfg_file = strdup(cfg_path);
            checkfile(opts->cfg_file);
            free(cfg_path);
            free(dir);
        }
    }

    if (opts->organelles != NULL) {
        if (strcmp(opts->organelles, "mt") != 0 && strcmp(opts->organelles, "pt") != 0) {
            log_message(ERROR, "Invalid organelles type (mt/pt)");
            exit(EXIT_FAILURE);
        }
        // if taxo is not plant, organelles should not be pt
        if (opts->taxo != 0 && strcmp(opts->organelles, "pt") == 0) {
            log_message(ERROR, "Invalid organelles type: %s for taxo type: %d", opts->organelles, opts->taxo);
            exit(EXIT_FAILURE);
        }
    } else {
        opts->organelles = strdup("mt");
    }

    if (which_executable("blastn") == 0) {
        log_message(ERROR, "Can't find blastn, please install it");
        exit(EXIT_FAILURE);
    }

    if (which_executable("apptainer") == 0 && which_executable("singularity") == 0) {
        log_message(ERROR, "Can't find apptainer or singularity, please install one of them");
        exit(EXIT_FAILURE);
    }
    // if (*organelles != NULL) {
    //     if (strcmp(*organelles, "mt") != 0 && strcmp(*organelles, "pt") != 0 && strcmp(*organelles, "all") != 0) {
    //         log_message(ERROR, "Invalid organelles type (mt/pt/all): %s", organelles);
    //         exit(EXIT_FAILURE);
    //     }
    // }

}


void graphBuild_arguments(int argc, char *argv[], graphBuildArgs *args) {

    int option_index = 0;
    static struct option long_options[] = {
        {"subsample", 1, 0, 'i'},
        {"graphinfo", 1, 0, 'a'},
        {"output", 1, 0, 'o'},
        {"organelles", 1, 0, 'G'},
        {"taxo", 1, 0, 'x'},
        {"depth", 1, 0, 'd'},
        {"seeds", 1, 0, 's'},
        {"cpu", 1, 0, 'T'},
        {"help", 0, 0, 'h'},
        {"version", 0, 0, 'v'},
        {0, 0, 0, 0}
    };

    while (1) {
        int c = getopt_long(argc, argv, ":i:a:o:G:x:d:s:T:mhv", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 'i': args->subsample = optarg; checkfile(args->subsample); break;
            case 'a': args->graphinfo = optarg; checkfile(args->graphinfo); break;
            case 'o': args->output_file = optarg; break;
            case 'G': args->organelles = optarg; break;
            case 'x': args->taxo = atoi(optarg); break;
            case 'd': args->depth = atof(optarg); break;
            case 's':
                while ((optind - 1) < argc && argv[(optind - 1)][0] != '-') {
                    if (is_numeric(argv[optind - 1])) {
                        args->seeds = realloc(args->seeds, (args->seedCount + 1) * sizeof(int));
                        args->seeds[args->seedCount] = atoi(argv[optind - 1]);
                        args->seedCount++;
                        optind++;
                    } else {
                        log_message(ERROR, "Invalid seed: %s", argv[optind - 1]);
                        graphBuild_usage();
                        exit(EXIT_FAILURE);
                    }
                }
                optind--;
                if (args->seedCount < 1) {
                    log_message(ERROR, "No seeds provided for -s option.");
                    graphBuild_usage();
                    exit(EXIT_FAILURE);
                }
                break;
            case 'T': args->cpu = atoi(optarg); break;
            case 'h': graphBuild_usage(); exit(EXIT_SUCCESS);
            case 'v': log_info("PMAT v%s\n", VERSION_PMAT); exit(EXIT_SUCCESS);
            case '?':
                log_message(ERROR, "Invalid option: %s", argv[optind-1]);
                graphBuild_usage();
                exit(EXIT_FAILURE);
            default:
                log_message(ERROR, "Unexpected error while parsing options");
                graphBuild_usage();
                exit(EXIT_FAILURE);
        }
    }

    // Validate required parameters
    if (args->subsample == NULL || args->graphinfo == NULL || args->output_file == NULL) {
        log_message(ERROR, "Missing required options");
        graphBuild_usage();
        exit(EXIT_FAILURE);
    } 

    if (args->subsample != NULL) {
        args->cutseq = malloc(strlen(args->subsample) + strlen("/PMAT_cut_seq.fa") + 1);
        strcpy(args->cutseq, args->subsample);
        strcat(args->cutseq, "/PMAT_cut_seq.fa");
        if (is_file(args->cutseq) == 0) {
            log_message(ERROR, "Input file does not exist: %s", args->cutseq);
            graphBuild_usage();
            exit(EXIT_FAILURE);
        }
    }

    if (args->graphinfo != NULL) {
        args->assembly_graph = malloc(strlen(args->graphinfo) + strlen("/PMATContigGraph.txt") + 1);
        strcpy(args->assembly_graph, args->graphinfo);
        strcat(args->assembly_graph, "/PMATContigGraph.txt");
        if (is_file(args->assembly_graph) == 0) {
            log_message(ERROR, "Input file does not exist: %s", args->assembly_graph);
            graphBuild_usage();
            exit(EXIT_FAILURE);
        }
        args->assembly_fna = malloc(strlen(args->graphinfo) + strlen("/PMATAllContigs.fna") + 1);
        strcpy(args->assembly_fna, args->graphinfo);
        strcat(args->assembly_fna, "/PMATAllContigs.fna");
        if (is_file(args->assembly_fna) == 0) {
            log_message(ERROR, "Input file does not exist: %s", args->assembly_fna);
            graphBuild_usage();
            exit(EXIT_FAILURE);
        }
    }

    if (args->cpu < 1) {
        log_message(ERROR, "Invalid cpu: %d", args->cpu);
        exit(EXIT_FAILURE);
    }

    if (args->organelles != NULL) {
        if (strcmp(args->organelles, "mt") != 0 && strcmp(args->organelles, "pt") != 0) {
            log_message(ERROR, "Invalid organelles type (mt/pt)");
            exit(EXIT_FAILURE);
        }
    } else {
        args->organelles = strdup("mt");
    }

    if (args->taxo != 0 && args->taxo != 1 && args->taxo != 2) {
        log_message(ERROR, "Invalid taxo type: %d", args->taxo);
        exit(EXIT_FAILURE);
    } else if (args->taxo == 1) {
        if (strcmp(args->organelles, "pt") == 0) {
            log_message(ERROR, "Invalid organelles type (pt)");
            exit(EXIT_FAILURE);
        }
    }

    if (which_executable("blastn") == 0) {
        log_message(ERROR, "Can't find blastn, please install it");
        exit(EXIT_FAILURE);
    }

}

int main(int argc, char *argv[]) {
    
    // char *exe_path = realpath(argv[0], NULL);
    // printf("exe_path: %s\n", exe_path);
    char *exe_path = pmat_path(argv[0]);
    if (!exe_path) {
        log_message(ERROR, "Failed to resolve executable path");
        exit(EXIT_FAILURE);
    }

    if (argc >= 2) 
    {
        /* Run autoMito */
        if (strcmp(argv[1], "autoMito") == 0) {

            autoMitoArgs optauto = {0};
            optauto.genomesize = NULL;
            optauto.runassembly = NULL;
            optauto.organelles = NULL;
            optauto.task = 1;
            optauto.taxo = 0;
            optauto.cfg_flag = 0;
            optauto.factor = 1;
            optauto.seed = 6;
            optauto.breaknum = 20000;
            optauto.mi = 90;
            optauto.ml = 40;
            optauto.cpu = 8;
            optauto.mem = 0;
            optauto.kmersize = 31;
            
            autoMito_arguments(argc - 1, argv + 1, exe_path, &optauto);
            
            log_message(INFO, "PMAT v%s", VERSION_PMAT);
            autoMito(exe_path, &optauto);

        } else if (strcmp(argv[1], "graphBuild") == 0) {

            graphBuildArgs optgraph = {0};
            optgraph.subsample = NULL;
            optgraph.graphinfo = NULL;
            optgraph.assembly_graph = NULL;
            optgraph.assembly_fna = NULL;
            optgraph.cutseq = NULL;
            optgraph.output_file = NULL;
            optgraph.organelles = NULL;
            optgraph.depth = -1;
            optgraph.seeds = NULL;
            optgraph.seedCount = 0;
            optgraph.taxo = 0;
            optgraph.cpu = 8;

            graphBuild_arguments(argc - 1, argv + 1, &optgraph);

            log_message(INFO, "PMAT v%s", VERSION_PMAT);
            graphBuild(exe_path, &optgraph);

        } else if (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0) {
            log_info("PMAT v%s\n", VERSION_PMAT);
            exit(EXIT_SUCCESS);
        } else if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
            usage();
            exit(EXIT_SUCCESS);
        } else {
            log_message(ERROR, "Invalid command '%s'", argv[1]);
            log_info("For more information, please visit https://github.com/aiPGAB/PMAT2\n");
            usage();
            exit(EXIT_FAILURE);
        }
    } else {
        usage();
        exit(EXIT_FAILURE);
    }
    
    /* End of program */
    log_message(INFO, "Task over. bye!");

    free(exe_path);
    


    return 0;
}