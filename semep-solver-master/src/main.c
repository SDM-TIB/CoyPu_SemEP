/**
 * Copyright (C) 2012-2017 Universidad Simón Bolívar
 * Copyright (C) 2018 Research Center L3S
 * 
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gvpalma@usb.ve>
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <assert.h>

#include "types.h"
#include "util.h"
#include "input.h"
#include "semEP.h"
#include "link_pred.h"

#define MIN_ENTITY_INPUT   3
#define MIN_TYPE           1
#define LAST_ARGS          1

struct global_args {
     int n;
     char *graph_filename;
     bool constraint;
     bool prediction;
     struct entity_input_array e_in; 
};

static struct global_args g_args;
static const char *optString = "pcn:";

/*********************************
 **  Parse Arguments
 *********************************/

static void display_usage(void)
{
     fatal("Incorrect arguments:\n\tsemEP [-c] [-p] -n <number of vertex types> <vertices type 1> <similarity matrix type 1> <threshold type 1> ... <vertices type n> <similarity matrix type n> <threshold type n> <graph>\n");
}

static char *get_name(const char *filename, int n, char *instance)
{
     char buf[n];
     char *aux = NULL;
     time_t raw_time;
     struct tm * time_info;

     time(&raw_time);
     time_info = localtime(&raw_time);

     strcpy(buf,filename);
     aux = strtok(buf,"/");
     do {
	  strcpy(instance, aux);
     } while((aux = strtok(NULL,"/")) != NULL);

     strcpy(buf,instance);
     aux = strtok(buf,".");
     strcpy(instance, aux);

     if (asprintf(&aux, "%s-%dh-%dm-%ds", instance, time_info->tm_hour, time_info->tm_min, time_info->tm_sec) == -1)
	  fatal("Error in the construction of the name for the output");
     return aux;
}

static void initialize_arguments(void)
{
     g_args.n = MIN_TYPE;
     init_struct_array(g_args.e_in);
     g_args.graph_filename = NULL;
     g_args.constraint = false;
     g_args.prediction = false;
}

static void print_args(void)
{
     unsigned int i, n;
     
     printf("\n**********************************************\n");
     printf("Parameters:\n");
     printf("Number of vertex type: %d\n", g_args.n);
     n = g_args.e_in.nr;
     assert((unsigned int)g_args.n == n);
     for (i = 0; i < n; i++) {
	  printf("\nVertex type: %d\n", (i+1));
	  printf("Vertices file name: %s\n", g_args.e_in.data[i].vertices_filename);
	  printf("Matrix file name: %s\n", g_args.e_in.data[i].matrix_filename);
	  printf("Threshold: %.3f\n", g_args.e_in.data[i].threshold);
     }
     printf("\nGraph file name: %s\n", g_args.graph_filename);
     printf("Relation constraints: %s\n", g_args.constraint ? "true" : "false");
     printf("Get predicted links: %s\n", g_args.prediction ? "true" : "false");
     printf("************************************************\n");
}

static void parse_args(int argc, char **argv)
{
     int i, opt;
     bool ntypes;
     unsigned int j, n;
     struct entity_input etmp;
     
     initialize_arguments();
     opt = getopt(argc, argv, optString);
     ntypes = false;
     i = 0;
     while(opt != -1) {
	  switch(opt) {
	  case 'c':
	       g_args.constraint = true;
	       break;
	  case 'n':
	       g_args.n = (int)strtol(optarg, (char **)NULL, 10);
	       if (g_args.n < MIN_TYPE)
		    display_usage();
	       ntypes = true;
	       i = optind;
	       if ((  argc - optind - (g_args.n * MIN_ENTITY_INPUT)) != LAST_ARGS)
		    display_usage();
	       init_entity_input_array(&g_args.e_in, (unsigned int)g_args.n);
	       n = (unsigned int)g_args.n;
	       for (j = 0; j < n; j++) {
		    etmp.vertices_filename = argv[i++];
		    etmp.matrix_filename = argv[i++];
		    etmp.threshold = strtod(argv[i++], (char **)NULL);
		    ARRAY_PUSH(g_args.e_in, etmp);
	       }
	       break;
	  case 'p':
	       g_args.prediction = true;
	       break;
	  case '?':
	       display_usage();
	       break;
	  default:
	       /* You won't actually get here. */
	       fatal("?? getopt returned character code 0%o ??\n", opt);
	  }
	  opt = getopt(argc, argv, optString);
     }
     if (!ntypes || ((  argc - optind - (g_args.n * MIN_ENTITY_INPUT)) != LAST_ARGS))
	  display_usage();
     g_args.graph_filename = argv[i];
}

/*********************************
 *********************************
 **
 **       Main section
 **
 *********************************
 **********************************/

int main(int argc, char **argv)
{
     int len;
     clock_t ti, tf;
     static char *instance;
     char *name;
     struct problem_data in;
     clusters_t *c;
     
     ti = clock();
     parse_args(argc, argv);
     print_args();
     len = strlen(g_args.graph_filename) + 1;
     instance = xcalloc(len, 1);
     name =  get_name(g_args.graph_filename, len, instance);
     printf("\n**** GO semEP! **** \n");
     in = get_input_data(g_args.e_in, g_args.graph_filename);
     c = semEP_solver(&in.entities, name, g_args.constraint, &in.graph);
     if (g_args.prediction) {
	  printf("\nStarting Link Prediction\n");
	  perform_link_prediction(&in.entities, &in.graph, &c->partitions, name);
     }
     tf = clock();
     printf("Total time %.3f secs\n", (double)(tf-ti)/CLOCKS_PER_SEC);
     free(instance);
     free(name);
     free_input_data(&in);
     free_entity_input_array(&g_args.e_in);
     free_clusters(c);
     free(c);
     
     return 0;
}
