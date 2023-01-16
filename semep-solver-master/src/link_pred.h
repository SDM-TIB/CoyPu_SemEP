/**
 * Copyright (C) 2013-2017 Universidad Simón Bolívar
 * Copyright (C) 2018 Research Center L3S
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <palma@l3s.de>
 */

#ifndef __LINKPRED_H
#define __LINKPRED_H

#include "hash_map.h"
#include "semEP.h"

void perform_link_prediction(const struct entity_array *entities,
			     const struct node_ptr_array *color_nodes,
			     const color_ptr_array_t *partitions, const char *graph_name);

#endif /* __LINKPRED_H */
