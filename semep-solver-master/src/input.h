/**
 * Copyright (C) 2013-2017 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gvpalma@usb.ve>
 */

#ifndef ___INPUT_H
#define ___INPUT_H

#include "types.h"
#include "memory.h"
#include "semEP.h"

struct entity_input {
     char *vertices_filename;
     char *matrix_filename;
     double threshold;
};

struct entity_input_array {
     unsigned int nr;
     unsigned int alloc;
     struct entity_input *data;
};
     
struct problem_data {
     struct entity_array entities;
     struct node_ptr_array graph;
};

void init_entity_input_array(struct entity_input_array *e, unsigned int n);

struct problem_data get_input_data(struct entity_input_array e_in, const char *graph_filename);

void free_input_data(struct problem_data *in);

void free_entity_input_array(struct entity_input_array *e);

#endif /* ___INPUT_H */
