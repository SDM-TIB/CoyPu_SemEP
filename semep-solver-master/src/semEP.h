/**
 * Copyright (C) 2013-2017 Universidad Simón Bolívar
 * Copyright (C) 2017 Research Center L3S 
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gvpalma@usb.ve>
 */

#ifndef __SEMEP_H
#define __SEMEP_H

#include "types.h"
#include "string.h"
#include "util.h"
#include "hash_iset.h"
#include "hash_map.h"

#define NO_ETYPE  -1

struct color_entry {
     double sim_entities;
     struct hash_iset entities;
     struct hash_item entry;
};

struct color {
     int id;
     double sim_between;
     double cDensity;
     
     struct int_array id_nodes;
     struct hash_map ce;         /* Hash table of struct color_entry */
};

struct entity_key {
     int type;
     unsigned int pos;
};

typedef struct color_ptr_array {
     unsigned nr;
     unsigned alloc;
     struct color **data;
} color_ptr_array_t;

struct node {
     int id;
     double sim;
     char *relation;
     struct entity_key e1;
     struct entity_key e2;
     struct color *cp;
};

struct node_ptr_array {
     unsigned int nr;
     unsigned int alloc;
     struct node *data;
};

struct entity {
     int etype;
     double threshold;
     struct matrix ematrix;
     struct string_array vertices;
};

struct entity_array {
     unsigned int nr;
     unsigned int alloc;
     struct entity *data;
};

typedef struct clusters {
     double nc;                    /* Value to minimizing */
     color_ptr_array_t partitions; /* Array with the clustes */
} clusters_t;


#ifdef PRGDEBUG
void print_entity_data(struct entity e);

void print_node_prt_array(struct node_ptr_array a);
#endif

void init_entity(struct entity *e);

void free_entity_array(struct entity_array *e);

void free_node_ptr_array(struct node_ptr_array *nptr);

clusters_t *semEP_solver(const struct entity_array *entities, const char *graph_name,
			 bool rel_constraint, struct node_ptr_array *color_nodes);

void free_clusters(clusters_t *c);

#endif /* __SEMEP_H */
