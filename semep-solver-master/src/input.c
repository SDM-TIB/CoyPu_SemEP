/**
 * Copyright (C) 2013-2017 Universidad Simón Bolívar
 * Copyright (C) 2017 Research Center L3S
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
#include <assert.h>

#include "input.h"
#include "string.h"
#include "hash_map.h"
#include "util.h"

/*********************************
 ** Contants
 *********************************/

#define BUFSZ          256
#define MIN_HASH_SZ    389

/*********************************
 ** Structures
 *********************************/

struct concept {
     int etype;
     unsigned int pos;
     struct hash_entry entry;
};

struct edge {
     char *vertex1;
     char *vertex2;
     char *relation;
     double weight;
};

struct graph_data {
     unsigned int nr;
     struct edge *edge_array;
};

/*************************************
 ** Utilities
 ************************************/

#ifdef PRGDEBUG
static void print_hash_term(struct hash_map *term_pos)
{
     struct concept *item;
     struct hash_entry *hentry;

     printf("\nMap term-position\n");
     printf("Number of terms: %u\n", term_pos->fill);
     hmap_for_each(hentry, term_pos) {
	  item = hash_entry(hentry, struct concept, entry);
	  printf("** %s %u %d\n", hentry->key, item->pos, item->etype);
     }
}
#endif

static void init_problem_data(struct problem_data *pd)
{
     init_struct_array(pd->entities);
     init_struct_array(pd->graph);
}

void init_entity_input_array(struct entity_input_array *e, unsigned int n)
{
     e->data = (struct entity_input *)xcalloc(n, sizeof(struct entity_input));
     e->alloc = n;
     e->nr = 0;
}

/*********************************
 ** Vertices processing
 *********************************/

static void terms_load_and_type(const char *filename, struct hash_map *vertex_table,
				struct entity *e, int etype, double threshold)
{
     FILE *f;
     char buf[BUFSZ];
     size_t last, len;
     unsigned int i, n;
     int l;
     struct concept *item;
     
     e->etype = etype;
     e->threshold = threshold; 

     /* Processing of the file with the vertices */
     f = fopen(filename, "r");
     if (!f)
	  fatal("no terms file specified, abort\n");
     if (fgets(buf, sizeof(buf), f) == NULL)
	  fatal("error reading file");
     errno = 0;
     n = strtol(buf, NULL, 10);
     if (errno)
	  fatal("error in the conversion of string to integer\n");
     assert(e->vertices.nr == 0);
     ALLOC_GROW(e->vertices.data, n+e->vertices.nr, e->vertices.alloc);
      for (i = 0; i < n; i++) {
	   /* Reading of the string */
	  if (fgets(buf, sizeof(buf), f) == NULL)
	       fatal("error reading file");

          /* Change of the last character */
	  len = strlen(buf);
	  last = len - 1;
	  if(buf[last] == '\n')
	       buf[last] = 0;

          /* Copying the string to the list entities */
	  l = asprintf(&e->vertices.data[e->vertices.nr],"%s", buf);
	  if (l == -1)
	       fatal("error in term copy");
	  e->vertices.nr++;

          /* Copying the string to the hash table */
	  item = xmalloc(sizeof(struct concept));
	  item->pos = i;
	  item->etype = etype;
	  if (hmap_add_if_not_member(vertex_table, &item->entry, buf, len) != NULL)
	       fatal("the term %s is repeated in the file %s\n", buf, filename);
     
     }
     fclose(f);
}

/*********************************
 ** Matrices processing
 *********************************/

static void similarity_matrix_load(const char *filename, struct matrix *m)
{
     FILE *f;
     struct char_array buf;
     int n, i, j;
     int ch;

     f = fopen(filename, "r");
     if (!f) {
	  fatal("No instance file specified, abort\n");
     }
     n = 0;
     init_char_array(&buf);
     ch = getc(f);
     errno = 0;
     /* read number of nodes and arcs */
     while((ch != '\n') && (ch != EOF)) {
	  add_char(&buf, ch);
	  ch = getc(f);
     }
     if (ch != EOF) {
	  add_char(&buf, '\0');
	  n = strtol(buf.data, NULL, 10);
	  if (errno)
	       fatal("error in the conversion of string to integer\n");
     } else {
	  fatal("error reading the matrix data file\n");
     }
     string_clean(&buf);
     m->n = n;
     m->data = double_matrix(0, n, 0, n);
     i = 0;
     j = 0;
     ch = getc(f);
     if (ch == EOF) {
	  fatal("error reading the matrix data file\n");
     }
     errno = 0;
     while (ch != EOF) {
	  if ((ch != ' ') && (ch != '\n')) {
	       add_char(&buf, ch);
	  } else {
	       add_char(&buf, '\0');
	       m->data[i][j] = strtod(buf.data, NULL);
	       if (errno)
		    fatal("error in the conversion of string to double\n");
	       if (ch == ' ') {
		    j++;
	       } else if (ch == '\n') {
		    i++;
		    j = 0;
	       } else {
		    fatal("unknown character");
	       }
	       string_clean(&buf);
	  }
	  ch = getc(f);
     }
     fclose(f);
     free_array(buf);
}

/*********************************
 ** Loading the entities data
 *********************************/

static void load_entities_data(struct entity_input_array in, struct entity_array *ea,
			       struct hash_map *vt) 
{
     unsigned int i, n;

     n = in.nr; /* Number of entities types */
#ifdef PRGDEBUG
     printf("Number of of entities types %u\n", n);
#endif
     assert(ea->nr == 0);
     ALLOC_GROW(ea->data, n+ea->nr, ea->alloc);
     for (i = 0; i < n; i++) {
	  printf("Proccessing the vertex file %s ...\n", in.data[i].vertices_filename);
	  init_entity(&ea->data[i]);
	  terms_load_and_type(in.data[i].vertices_filename, vt,
			      &ea->data[i], i, in.data[i].threshold);
	  ea->nr++;
#ifdef PRGDEBUG
	  print_entity_data(ea->data[i]);
#endif
     }
}

static void load_similarities_matrices(struct entity_input_array in, struct entity_array *ea)
{
     unsigned int i, n;

     n = in.nr; /* Number of entities types */
     for (i = 0; i < n; i++) {
	  printf("Proccessing the matrix file %s ...\n", in.data[i].matrix_filename);
	  similarity_matrix_load(in.data[i].matrix_filename, &ea->data[i].ematrix);
#ifdef PRGDEBUG
	  print_matrix(ea->data[i].ematrix);
#endif
     }
}

/*********************************
 ** Graph processing
 *********************************/

static void graph_load(const char *filename, const struct hash_map *vt,
		       struct node_ptr_array *g)
{
     FILE *f;
     int n, i, ch, tok;
     size_t len;
     struct char_array buf;
     struct concept *item;
     struct hash_entry *hentry;
          
     f = fopen(filename, "r");
     if (!f) {
	  fatal("no graph file specified, abort\n");
     }
     n = 0;
     init_char_array(&buf);
     ch = getc(f);
     errno = 0;

     /* Reading the number of edges */
     while((ch != '\n') && (ch != EOF)) {
	  add_char(&buf, ch);
	  ch = getc(f);
     }
     if (ch != EOF) {
	  add_char(&buf, '\0');
	  n = strtol(buf.data, NULL, 10);
	  if (errno)
	       fatal("error in the conversion of string to integer\n");
     } else {
	  fatal("error reading the description data file\n");
     }
     string_clean(&buf);
     printf("\nNumber of edges of the input graph: %d\n", n);

     /* Allocating memory to be used */
     assert(g->nr == 0);
     ALLOC_GROW(g->data, n+g->nr, g->alloc);
     
     /* Reading the edges */
     ch = getc(f);
     if (ch == EOF) {
	  fatal("error reading the graph data file\n");
     }
     tok = 1;
     i = 0;
     while ((ch != EOF) && (i < n)) {
	  if ((ch != '\t') && (ch != '\n')) {
	       add_char(&buf, ch);
	  } else {
	       add_char(&buf, '\0');
	       if (ch == '\t') {
		    len = strlen(buf.data);
		    assert(len >= 1);

		    if ( (tok == 1) || (tok == 2) ) {
		    /* Getting the entity data */
			 hentry = hmap_find_member(vt, buf.data, len+1); 
			 if (hentry == NULL)
			      fatal("the vertex %s is not in the list of entities\n", buf.data);
			 item = hash_entry(hentry, struct concept, entry);
			 if (item == NULL)
			      fatal("error getting a valid vertex");
			 
			 /* Adding the entity data to the graph */
			 if (tok == 1) {
			      g->data[i].e1.type = item->etype;
			      g->data[i].e1.pos = item->pos;
			 } else if (tok == 2) {
			      g->data[i].e2.type = item->etype;
			      g->data[i].e2.pos = item->pos;
			 }
		    } else if (tok == 3) {
			 g->data[i].relation = xcalloc(len+1, 1);
			 strcpy(g->data[i].relation, buf.data);
		    } else {
			 fatal("error in the format of the graph %s\n", filename);
		    }
		    tok++;
	       } else {
		    /* Getting the edge weight */
		    assert(ch == '\n');
		    if (tok != 4)
			 fatal("error in the format of the graph %s\n", filename);
		    errno = 0;
		    g->data[i].sim = strtod(buf.data, NULL);
		    if (errno)
			 fatal("error in the conversion of string to integer\n");

		    /* Updating the array counter */
		    assert(g->nr == (unsigned int)i);
		    g->data[i].id = i; 
		    g->data[i].cp = NULL; 
		    g->nr++;

		    /* Updating the token and line counters */
		    i++;
		    tok = 1;
	       }
	       string_clean(&buf);
	  }
	  ch = getc(f);
     }
     fclose(f);
     free_array(buf);     
}

void free_input_data(struct problem_data *in)
{
     free_entity_array(&in->entities);
     free_node_ptr_array(&in->graph);    
}

void free_entity_input_array(struct entity_input_array *e)
{
     if (e->data)
	  free(e->data);
     e->alloc = 0;
     e->nr = 0;
}

static void free_entities_table(struct hash_map *term_pos)
{
     struct concept *item;
     struct hash_entry *hentry;
     struct hlist_node *n;

     hmap_for_each_safe(hentry, n, term_pos) {
	  item = hash_entry(hentry, struct concept, entry);
	  hmap_delete(term_pos, hentry);
	  free(item);
     }
     hmap_destroy(term_pos);
}
 
/*********************************
*********************************
**
** Getting the Problem Data
**
*********************************
*********************************/

struct problem_data get_input_data(struct entity_input_array e_in,
				   const char *graph_filename)
{
     struct problem_data in;
     struct hash_map vertex_table;

     hmap_create(&vertex_table, MIN_HASH_SZ);
     init_problem_data(&in);
     load_entities_data(e_in, &in.entities, &vertex_table);
     load_similarities_matrices(e_in, &in.entities);
     init_struct_array(in.graph);
     graph_load(graph_filename, &vertex_table, &in.graph);
#ifdef PRGDEBUG
     print_hash_term(&vertex_table);
     print_node_prt_array(in.graph);
#endif
     free_entities_table(&vertex_table);
     return in;
}
