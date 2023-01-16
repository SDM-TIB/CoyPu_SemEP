/**
 * Copyright (C) 2013-2017 Universidad Simón Bolívar 
 * Copyright (C) 2018 Research Center L3S
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <palma@l3s.de>
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "memory.h"
#include "hash_map.h"
#include "semEP.h"
#include "link_pred.h"

/*********************************
 ** Constants
 *********************************/

#define NOCLUSTER    -1

/*********************************
 ** Structures
 *********************************/

typedef struct item_entry {
     struct hash_entry entry;
} item_entry_t;

typedef struct relation_entry {
     char *relation;
     struct hash_entry entry;
} relation_entry_t;

struct vertice {
     int type;
     unsigned int pos;
};

struct vertice_array {
     unsigned int nr;
     unsigned int alloc;
     struct vertice *data;
};

typedef struct prediction {
     int cluster;
     struct vertice entity1;
     struct vertice entity2;
     char *relation;
     double prob;
} prediction_t;

typedef struct prediction_array {
     unsigned nr;
     unsigned alloc;
     prediction_t *data;
} prediction_array_t;

/*************************************
 *************************************
 **
 ** Functions and Procedures
 **
 ************************************
 ************************************/

static void get_list_of_vertices(struct color *cluster,
				 struct vertice_array *lnodes)
{
     unsigned int i, alloc;
     int type;
     struct vertice vtmp;
     struct hash_item *hentry;
     struct color_entry *color_item;
     int *table_nodes;
     
     hmap_for_each_int(hentry, &cluster->ce) {
	  color_item = hash_entry(hentry, struct color_entry, entry);
	  type = hentry->key;
	  table_nodes = color_item->entities.htable; 
	  alloc = color_item->entities.alloc;
	  for (i = 0; i < alloc; i++) {
	       if (table_nodes[i] != NOELEM) {
		    vtmp.type = type;
		    vtmp.pos = table_nodes[i];
		    ARRAY_PUSH(*lnodes, vtmp);
	       }
	  }
     }
}

#ifdef PRGDEBUG
static void check_repeated_entities(const struct vertice_array *lnodes)
{
     unsigned int i, j, n;

     n = lnodes->nr;
     for (i = 0; i < n-1; i++) {
	  for (j = i+1; j < n; j++) {
	       if ( (lnodes->data[i].type == lnodes->data[j].type) &&
		    (lnodes->data[i].pos == lnodes->data[j].pos) ) {
		    fatal("Repetead nodes in the vector");
	       }
	  }
     }
}
#endif

static double get_cluster_probability(unsigned int n_nodes, unsigned int n_edges)
{
     if (n_nodes == 0)
	  fatal("Zero Division");
     return (double) (2.0 * n_edges) / (n_nodes * (n_nodes-1));
}

static void get_predicted_links(const color_ptr_array_t *partitions,
				struct hash_map *index_edges,
				struct hash_map *index_relations,
				prediction_array_t *cluster_pred)
{
     char *key_e, *key_r;
     int  t1, t2, l;
     prediction_t ptemp;
     unsigned int k, i, j, n_nodes, n_clusters, p1, p2, n_edges;
     struct vertice_array lnodes;
     struct color *cluster;
     relation_entry_t *ritem;
     struct hash_entry *hentry;
	  
     n_clusters = partitions->nr;
     for (k = 0; k < n_clusters; k++) {
	  cluster = partitions->data[k];
	  n_edges = cluster->id_nodes.nr; 
	  init_struct_array(lnodes);
	  get_list_of_vertices(cluster, &lnodes);
#ifdef PRGDEBUG
	  check_repeated_entities(&lnodes);
#endif
	  n_nodes = lnodes.nr;
	  for (i = 0; i < n_nodes-1; i++) {
	       t1 = lnodes.data[i].type;
	       p1 = lnodes.data[i].pos;
	       for (j = i+1; j < n_nodes; j++) {
		    t2 = lnodes.data[j].type;
		    p2 = lnodes.data[j].pos;
		    /* Check if the edge exists */
		    l = asprintf(&key_e, "%d-%u-%d-%u", t1, p1, t2, p2);
		    if ( l == -1 )
			 fatal("Error in key creation");
		    if (hmap_is_member(index_edges, key_e, l) == 0) {
			 /* Check if a valid relation */
			 l = asprintf(&key_r, "%d-%d", t1, t2);
			 if ( l == -1 )
			      fatal("Error creating a relation key");
			 hentry = hmap_find_member(index_relations, key_r, l);
			 if (hentry != NULL) {
			      /* we found a new prediction */
			      ritem = hash_entry(hentry, relation_entry_t, entry);
			      assert((unsigned int)cluster->id == k);
			      ptemp.cluster = cluster->id;
			      ptemp.entity1 = lnodes.data[i];
			      ptemp.entity2 = lnodes.data[j];
			      ptemp.relation = ritem->relation;
			      ptemp.prob = get_cluster_probability(n_nodes, n_edges);
			      ARRAY_PUSH(*cluster_pred, ptemp);
			 }
			 free(key_r);
		    }
		    free(key_e);
	       }
	  }
	  free_array(lnodes);
     }
}

static inline void add_relation(struct hash_map *r_index, unsigned int t1,
				unsigned int t2, char *relation)
{
     int l;
     char *buf;
     relation_entry_t *ritem;

     l = asprintf(&buf, "%u-%u", t1, t2);
     if ( l == -1 )
	  fatal("Error in relation key creation");
     if (hmap_find_member(r_index, buf, l) == NULL) {
	  ritem = (relation_entry_t *)xmalloc(sizeof(relation_entry_t));
	  ritem->relation = relation;
	  if (hmap_add(r_index, &ritem->entry, buf, l) != 0)
	       fatal("Error in adding in the relation index\n");
     }
     free(buf);
}

static inline void add_arc(struct hash_map *index_edges, unsigned int t1, int p1,
			   unsigned int t2, int p2)
{
     int l;
     char *buf;
     item_entry_t *item;
     
     l = asprintf(&buf, "%u-%d-%u-%d", t1, p1, t2, p2);
     if ( l == -1 )
	  fatal("Error in edge key creation");
     item = (item_entry_t *)xmalloc(sizeof(item_entry_t));
     if (hmap_add_if_not_member(index_edges, &item->entry, buf, l) != NULL)
	  fatal("Error, repeated edge in the graph\n");
     free(buf);
}

static void build_indexes(const struct node_ptr_array *color_nodes,
			  struct hash_map *index_edges,
			  struct hash_map *index_relations)
{
     unsigned int i, n_edges, t1, t2;
     int p1, p2;
        
     n_edges = color_nodes->nr;
     for (i = 0; i < n_edges; i++) {
	  
	  /* We obtain the edge information */
	  t1 = color_nodes->data[i].e1.type;
	  p1 = color_nodes->data[i].e1.pos;
	  t2 = color_nodes->data[i].e2.type;
	  p2 = color_nodes->data[i].e2.pos;

	   /* We add two arcs to edge index */
	  add_arc(index_edges, t1, p1, t2, p2);
	  add_arc(index_edges, t2, p2, t1, p1);
	  	 
	  /* We add a relationship to relation index */
	  add_relation(index_relations, t1, t2, color_nodes->data[i].relation);
	  add_relation(index_relations, t2, t1, color_nodes->data[i].relation);
     }
}

static void print_predicted_links(const prediction_array_t *cluster_pred,
				  const struct entity_array *entities,
				  const char *name)
{
     unsigned i, n;
     prediction_t ptemp;
     FILE *f;
     char *output, *s1, *s2;
     int current;
   
     if (asprintf(&output, "%s-Predictions.txt", name) == -1)
	  fatal("Error in prediction file");
     f = fopen(output, "w");
     printf("File with the predictions: %s\n", output);
     free(output);
     if (!f)
	  fatal("No descriptor file specified, abort\n");
     n = cluster_pred->nr;
     printf("Number of links predicted: %d\n", n);
     current = NOCLUSTER;
     for (i = 0; i < n; i++) {
	  ptemp = cluster_pred->data[i];
	  if (current != ptemp.cluster) {
	       current = ptemp.cluster;
	       fprintf(f, "Cluster\t%u\n", current);
	  }
	  s1 = entities->data[ptemp.entity1.type].vertices.data[ptemp.entity1.pos];
	  s2 = entities->data[ptemp.entity2.type].vertices.data[ptemp.entity2.pos];
	  fprintf(f, "%s\t%s\t%s\t%.4f\n", s1, ptemp.relation, s2, ptemp.prob);
     }
     fclose(f);
}

/*************************************
 *************************************
 **
 **  Freeing memory
 **
 ************************************
 ************************************/

static void free_edges_table(struct hash_map *edges_obs)
{
     item_entry_t *item;
     struct hash_entry *hentry;
     struct hlist_node *n;

     hmap_for_each_safe(hentry, n, edges_obs) {
	  item = hash_entry(hentry, item_entry_t, entry);
	  hmap_delete(edges_obs, hentry);
	  free(item);
     }
     hmap_destroy(edges_obs);
}

static void free_relation_table(struct hash_map *r_obs)
{
     relation_entry_t *item;
     struct hash_entry *hentry;
     struct hlist_node *n;

     hmap_for_each_safe(hentry, n, r_obs) {
	  item = hash_entry(hentry, relation_entry_t, entry);
	  hmap_delete(r_obs, hentry);
	  free(item);
     }
     hmap_destroy(r_obs);
}

/*************************************
 *************************************
 **
 ** Link Prediction main procedure
 **
 ************************************
 ************************************/

void perform_link_prediction(const struct entity_array *entities,
			     const struct node_ptr_array *color_nodes,
			     const color_ptr_array_t  *partitions, const char *name)
{
     prediction_array_t cluster_pred;
     struct hash_map edges_obs, relations_obs;
     unsigned int n_edges; 

     init_struct_array(cluster_pred);
     n_edges = color_nodes->nr;
     hmap_create(&edges_obs, 2*n_edges);
     hmap_create(&relations_obs, n_edges/2);
     build_indexes(color_nodes, &edges_obs, &relations_obs);
     get_predicted_links(partitions, &edges_obs, &relations_obs, &cluster_pred);
     print_predicted_links(&cluster_pred, entities, name);
     free_array(cluster_pred);
     free_edges_table(&edges_obs);
     free_relation_table(&relations_obs);
}
