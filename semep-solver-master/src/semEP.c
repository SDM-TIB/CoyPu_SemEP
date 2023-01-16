/**
 * Copyright (C) 2013-2017 Universidad Simón Bolívar
 * Copyright (C) 2017-2019 Research Center L3S
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <palma@l3s.de>
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include "dlist.h"
#include "graph_adj.h"
#include "memory.h"
#include "hash_iset.h"
#include "hash_map.h"
#include "util.h"
#include "semEP.h"

/*********************************
 ** Constants
 *********************************/

#define NOCOLOR       -1
#define NONODE        -1
#define NOCLUSTER     -1
#define INFTY         INT_MAX
#define EPSILON       10e-7
#define NO_ETYPE      -1
#define HT_SIM_SZ     7
#define ONE           1.0001
#define NOSIM         -999.999999
#define COLOR_MARK    5000
#define EDGE_MARK     20000000
#define MILLION       1000000

/************************************************************
 ** General macros
 ************************************************************/

#define CHECK_ELEM(e)						\
     do {							\
	  if ((e) == NOELEM)					\
	       fatal("Error in searching a valid element\n");	\
     } while(0)

/************************************************************
 ** Macros used for the binary heap implementation
 ************************************************************/

#define PARENT(i)     ((int)(((i) - 1) / 2))
#define LEFT(i)       (((i) * 2) + 1)
#define RIGHT(i)      (((i) * 2) + 2)

/*********************************
 ** Structures
 *********************************/

typedef struct saturation {
     int node;
     int n_adj_colors;
     struct hash_iset color_used;
} saturation_t;

typedef struct pqueue {
     int size;
     int *node_pos;
     saturation_t **heap;
} pqueue_t; 

typedef struct sim_cluster_entities {
     int type;
     double sim;
} sim_cluster_entities_t;

typedef struct info_partition {
     bool is_new_color;
     int color;
     double nc;
     double sim_between;
     double pairwise_sim;
     double cDensity;
     bool eq_type; /* Pair of entities of a edge has the same type */ 

     /* Maximum number of types on one edge is two */
     sim_cluster_entities_t sim_entity_1; 
     sim_cluster_entities_t sim_entity_2;
} info_partition_t;

typedef struct sim_tuple {
     double total;
     double sim;
} sim_tuple_t;

/*************************************
 *************************************
 **
 ** Utilities 
 **
 ************************************
 ************************************/

void init_entity(struct entity *e)
{
     e->etype = NO_ETYPE;
     e->threshold = 0.0;
     init_struct_array(e->vertices);
     e->ematrix.n  = 0;
     e->ematrix.data = NULL;
}

static void init_clusters(clusters_t *c)
{
     c->partitions.alloc = 0;
     c->partitions.nr = 0.0;
     c->partitions.data = NULL;
}

#ifdef PRGDEBUG
void print_entity_data(struct entity e)
{
     printf("\n-----------\n");
     printf("Entity type %d\n", e.etype);
     printf("Threshold %.3f\n", e.threshold);
     print_string_array(&e.vertices);
     if (e.ematrix.data != NULL)
	  print_matrix(e.ematrix);
     printf("-----------\n");
}

void print_node_prt_array(struct node_ptr_array a)
{
     printf("\nNumber of elements: %u\n", a.nr);
     for (unsigned i = 0; i < a.nr; i++) {
	  printf("id %d type1 %d pos1 %u type2 %d pos2 %u sim %.3f relation %s\n",
		 a.data[i].id, a.data[i].e1.type, a.data[i].e1.pos,
		 a.data[i].e2.type, a.data[i].e2.pos,
		 a.data[i].sim, a.data[i].relation);
     }
}

static void print_coloring(clusters_t *clusters)
{
     struct color *c;
     unsigned int i, n;

     n = clusters->partitions.nr;
     printf("\nResults of the partitioning:\n");
     printf("Number of partitions: %u\n", n);
     printf("nc %.4f\n", clusters->nc);
     for (i = 0; i < n; i++) {
	  c = clusters->partitions.data[i];
	  assert((unsigned)c->id == i);
	  printf("id %d - cDensity %.4f - sim_btw %.4f - id_nodes %u - num. types %u\n",
		 c->id, c->cDensity, c->sim_between , c->id_nodes.nr, c->ce.fill);
     }
     printf("\n");
}
#endif

static inline bool eq_double(double x, double y)
{
     return (fabs(x - y) < EPSILON);
}

static inline void init_info_partition(info_partition_t *ip)
{
     ip->is_new_color = false;
     ip->color = NOCOLOR;
     ip->nc = NOSIM;
     ip->sim_between = NOSIM;
     ip->pairwise_sim = NOSIM;
     ip->cDensity = NOSIM;
     ip->eq_type = false;
     ip->sim_entity_1.type = NONODE;
     ip->sim_entity_1.sim = NOSIM; 
     ip->sim_entity_2.type = NONODE;
     ip->sim_entity_2.sim = NOSIM;
}

/*************************************
 *************************************
 **
 ** Similarity functions
 **
 ************************************
 ************************************/

/**
 *  Function to get the similirity between entities
 */
static inline double similarity(const struct matrix *ematrix,
				unsigned int pos1, unsigned int pos2)
{
     if ( (ematrix->n <= pos1) || (ematrix->n <= pos2) ) 
	  fatal("matrix similarity index out of bounds");
     return ematrix->data[pos1][pos2];
}

/***************************************
****************************************
**
** Contraints functions and 
** Mapping to Graph Coloring Problem
**
**************************************
**************************************/

/**
 *  Determine whether two edges can be in the same cluster, based on the
 *  similarity of the entities. 
 */
static inline bool edge_constraint(const struct entity_array *ea, const struct node *n1,
				   const struct node *n2) 
{
     if (n1->e1.type == n2->e1.type) 
	  if (MAX(ea->data[n1->e1.type].ematrix.data[n1->e1.pos][n2->e1.pos],
		  ea->data[n1->e1.type].ematrix.data[n2->e1.pos][n1->e1.pos]) <=
	      ea->data[n1->e1.type].threshold)
	       return true;
     if (n1->e1.type == n2->e2.type) 
	  if (MAX(ea->data[n1->e1.type].ematrix.data[n1->e1.pos][n2->e2.pos],
		  ea->data[n1->e1.type].ematrix.data[n2->e2.pos][n1->e1.pos]) <=
	      ea->data[n1->e1.type].threshold)
	       return true;
     if (n1->e2.type == n2->e1.type) 
	  if (MAX(ea->data[n1->e2.type].ematrix.data[n1->e2.pos][n2->e1.pos],
		  ea->data[n1->e2.type].ematrix.data[n2->e1.pos][n1->e2.pos]) <=
	      ea->data[n1->e2.type].threshold)
	       return true;
     if (n1->e2.type == n2->e2.type)
	  if (MAX(ea->data[n1->e2.type].ematrix.data[n1->e2.pos][n2->e2.pos],
		  ea->data[n1->e2.type].ematrix.data[n2->e2.pos][n1->e2.pos]) <=
	      ea->data[n1->e2.type].threshold)
	       return true;
     
     return false;
}

/**
 *  Determine whether two edges can be in the same cluster,
 *  only if they have the same relationship type.
 */
static inline bool relation_constraint(const struct node *n1, const struct node *n2)
{
     return (strcmp(n1->relation, n2->relation) != 0);
}

/**
 *  Build the graph to coloring
 */
static void build_graph_to_coloring_matrix(const struct entity_array *entities,
					   bool rel_constr,
					   struct graph_adj *gc,
					   struct node_ptr_array *vn)
{
     long cont;
     unsigned int i, j, n;
     struct node *x, *y;
     bool e_constraint, r_constraint;

     cont = 0;
     n = vn->nr;

     if (rel_constr) {
	  for (i = 0; i < n-1; i++) {
	       x = &vn->data[i];
	       for (j = i+1; j < n; j++) {
		    y = &vn->data[j];
		    e_constraint = edge_constraint(entities, x, y);
		    r_constraint = relation_constraint(x, y);
		    if (e_constraint ||  r_constraint) {
			 add_arc(gc, cont, i, j);
			 add_arc(gc, cont, j, i);
			 cont++;
			 if (((cont % EDGE_MARK) == 0) && (cont >= EDGE_MARK))
			      printf("Number of edges in the graph to coloring so far: %ld MM\n",
				     cont/MILLION);
		    }
	       }
	  }
     } else {
	  for (i = 0; i < n-1; i++) {
	       x = &vn->data[i];
	       for (j = i+1; j < n; j++) {
		    y = &vn->data[j];
		    e_constraint = edge_constraint(entities, x, y);
		    if (e_constraint) {
			 add_arc(gc, cont, i, j);
			 add_arc(gc, cont, j, i);
			 cont++;
			 if (((cont % EDGE_MARK) == 0) && (cont >= EDGE_MARK))
			      printf("Number of edges in the graph to coloring so far: %ld MM\n",
				     cont/MILLION);
		    }
	       }
	  }
     }
}

/*************************************
 *************************************
 **
 ** Checking and control
 **
 ************************************
 ************************************/

/**
 * Pairwise symilarity for elements of a same type
 */
static double pairwise_sim(const struct hash_iset *entities,
			   const struct matrix *ematrix)
{
     double total_sim;
     unsigned int i, j, n, k, l, alloc;
     int a, b;
     int *data;
   
     data = entities->htable;
     n = entities->nr;
     alloc = entities->alloc;
     if (n == 0) {
	  total_sim = 0.0;
     } else if (n == 1) {
	  k = 0;
	  a = NOELEM;
	  for (k = 0; k < alloc; k++) {
	       if (data[k] != NOELEM) {
		    a = data[k];
		    break;
	       }
	  }
	  CHECK_ELEM(a);
	  total_sim = similarity(ematrix, a, a);
     } else {
	  total_sim = 0.0;
	  k = 0;
	  i = 0;
	  while ((i < n-1) && (k < alloc)) {
	       if (data[k] != NOELEM) {
		    a = data[k];
		    j = i+1;
		    l = k+1; 
		    while ((j < n) && (l < alloc)) {
			 if (data[l] != NOELEM) {
			      b = data[l];
			      assert(a != b);
			      total_sim += similarity(ematrix, a, b);
			      j++;
			 }
			 l++;
		    }
		    assert(j == n);
		    i++;
	       }
	       k++;
	  }
	  assert(i == (n-1)); 
     }
     return total_sim;
}


/**
 * Cluster density
 */
static double cDensity(struct color *c, const struct entity_array *ea)
{
     double cdensity, sim_acc, bpe;
     unsigned int n, nr;
     struct color_entry *ctmp;
     struct hash_item *hentry;
     
     n = c->id_nodes.nr;
     if (n == 0) {
	  error("Warning, processing a empty partition\n");
	  cdensity = 0.0;
     } else if (n == 1){
	  cdensity = EPSILON;
     } else {
	  cdensity = 0.0;
	  nr = c->ce.fill;
	  sim_acc = 0.0;
	  hmap_for_each_int(hentry, &c->ce) { 
	       ctmp = hash_entry(hentry, struct color_entry, entry);
	       n = ctmp->entities.nr;
	       assert(n > 0);
	       if (n > 1)
		    sim_acc += pairwise_sim(&ctmp->entities,
					    &ea->data[hentry->key].ematrix)/((n*n-n)/2.0);
	       else 
		    sim_acc += pairwise_sim(&ctmp->entities,
					    &ea->data[hentry->key].ematrix);
	  }
	  bpe = (double)c->sim_between/c->id_nodes.nr;
	  cdensity = (bpe + sim_acc)/(1.0 + nr);
     }
     return cdensity;
}

#ifdef PRGDEBUG
/**
 * Function to minimizing for the clustering
 */
static double color_density(const clusters_t *c, const struct entity_array *ea)
{
     unsigned int i, n_colors;
     double nc;
     
     nc = 0.0;
     n_colors = c->partitions.nr;
     for (i = 0; i < n_colors; i++) {
	  nc += (i+1) - cDensity(c->partitions.data[i], ea);
     }
     return nc;
}
#endif

/*************************************
 *************************************
 **
 **  Priority Queue
 **
 ************************************
 ************************************/

static inline int compare_saturation(const struct graph_adj *g, const saturation_t *a,
				     const saturation_t *b)
{
     int r;

     r = 0;
     if (a->n_adj_colors > b->n_adj_colors) { 
	  r = 1;
     } else if (a->n_adj_colors < b->n_adj_colors) {
	  r = -1;
     } else {
	  if (g->degree[a->node].din > g->degree[b->node].din) {
	       r = 1;
	  } else if (g->degree[a->node].din < g->degree[b->node].din) { 
	       r = -1;
	  } else {
	       /* Numerical order */
	       if (a->node < b->node) {
		    r = 1;
	       } else if (a->node > b->node) {
		    r = -1;
	       } else {
		    fatal("Same pair the nodes in the comparison");
	       }
	  }
     }
     return r;
}

static inline void free_saturation_node(saturation_t *node)
{
     if (node) {
	  free_hash_iset(&node->color_used);
	  free(node);
     }
}

static inline void pq_init(pqueue_t *pq)
{
     pq->size = 0;
     pq->heap = NULL;
     pq->node_pos = NULL;
}

static inline void pq_delete(pqueue_t *pq)
{
     int i;
     
     for(i = 0; i < pq->size; i++)
	  free_saturation_node(pq->heap[i]);
     free(pq->node_pos);
     free(pq->heap);
}

static inline void pq_insert(const struct graph_adj *g, pqueue_t *pq,
			     saturation_t *node)
{
     int i, p;
     saturation_t **tmp;

     tmp = xrealloc(pq->heap, (pq->size+1)*sizeof(saturation_t *));
     pq->heap = tmp;
     pq->heap[pq->size] = node;
     i = pq->size;
     p = PARENT(i);
     while((i > 0) &&  (compare_saturation(g, pq->heap[p], pq->heap[i]) < 0)){
	  SWAP(pq->heap[p], pq->heap[i]);
	  i = p;
	  p = PARENT(i);
     }
     pq->size++;
}

static int extract_max(const struct graph_adj *g, pqueue_t *pq, saturation_t **node)
{
     int i, j, l, r;
     saturation_t *aux;
     saturation_t **tmp;
          
     if(pq->size == 0)
	  return -1;

     *node = pq->heap[0];
     aux =  pq->heap[pq->size-1];
     SWAP(pq->node_pos[pq->heap[0]->node],
	  pq->node_pos[pq->heap[pq->size-1]->node]); /* SWAP the positions*/
     if((pq->size - 1) > 0){
	  tmp = (saturation_t **)xrealloc(pq->heap, (pq->size-1)*sizeof(saturation_t *));
	  pq->heap = tmp;
	  pq->size--;
     } else {
	  free(pq->heap);
	  pq->heap = NULL;
	  free(pq->node_pos);
	  pq->node_pos = NULL;
	  pq->size = 0;
	  return 0;
     }
     pq->heap[0] = aux;
     i = 0;
     while (true) {
	  l = LEFT(i);
	  r = RIGHT(i);
	  if((l < pq->size) && (compare_saturation(g, pq->heap[l], pq->heap[i]) > 0))
	       j = l;
	  else
	       j = i;

	  if((r < pq->size) && (compare_saturation(g, pq->heap[r], pq->heap[j]) > 0))
	       j = r;

	  if( j == i ) {
	       break;
	  } else {
	       SWAP(pq->node_pos[pq->heap[j]->node],
		    pq->node_pos[pq->heap[i]->node]); /* SWAP the positions*/
	       SWAP(pq->heap[j], pq->heap[i]);
	       i = j;
	  }
     }
     return 0;
}

static int increase_key(const struct graph_adj *g, pqueue_t *pq, int node, int color)
{
     int i, p, pos;
	  
     if (pq->size == 0)
	  return -1;

     pos = pq->node_pos[node];
     if (pos >= pq->size)
	  pos = -1;
   
     if (pos == -1)
	  return -2;

     if (insert_hash_iset(&(pq->heap[pos]->color_used), color))
	  pq->heap[pos]->n_adj_colors++;
     else
	  return 0;
     
     i = pos;
     p = PARENT(i);
     while((i > 0) && (compare_saturation(g, pq->heap[p], pq->heap[i]) < 0)){
	  SWAP(pq->node_pos[pq->heap[p]->node],
	       pq->node_pos[pq->heap[i]->node]); /* SWAP the positions*/
	  SWAP(pq->heap[p], pq->heap[i]);
	  i = p;
	  p = PARENT(i);
     }
     return 0;
}

static void init_saturation_pq(const struct graph_adj *g, pqueue_t *pq)
{
     int i, n;
     saturation_t *ns;
     
     n = g->n_nodes;
     for (i = 0; i < n; i++) {
	  ns = (saturation_t *)xmalloc(sizeof(saturation_t));
	  ns->node = i;
	  ns->n_adj_colors = 0;
	  init_hash_iset(&ns->color_used);
	  pq_insert(g, pq, ns);	  
     }
     assert(pq->size == n);
     pq->node_pos = (int *)xmalloc(n*sizeof(int));
     for (i = 0; i < n; i++) {
	  pq->node_pos[pq->heap[i]->node] = i;
     }
}

/*************************************
 *************************************
 **
 **  Coloration Solver
 **
 ************************************
 ************************************/

static inline struct color *new_color(int id, double sim_bt, int node,
				      struct entity_key e1, struct entity_key e2,
				      const struct entity_array *ea)
{
     struct color *new;
     struct color_entry *new_entry;
	  
     new = xcalloc(1, sizeof(struct color));
     new->id = id;
     new->sim_between = sim_bt;
     init_struct_array(new->id_nodes);
     ARRAY_PUSH(new->id_nodes, node);
     hmap_create(&new->ce, HT_SIM_SZ);
     new_entry = xcalloc(1, sizeof(struct color_entry));
     init_hash_iset(&new_entry->entities);
     insert_hash_iset(&new_entry->entities, e1.pos);
     if (e1.type == e2.type) {
          /* We have a only type of entity in the cluster */ 
	  insert_hash_iset(&new_entry->entities, e2.pos);
	  new_entry->sim_entities = pairwise_sim(&new_entry->entities,
						 &ea->data[e1.type].ematrix); 
	  hmap_add_int(&new->ce, &new_entry->entry, e1.type);
     } else {
	  /* Similarity with 1 nodes */
	  new_entry->sim_entities = pairwise_sim(&new_entry->entities,
						 &ea->data[e1.type].ematrix); 
	  hmap_add_int(&new->ce, &new_entry->entry, e1.type);

	  /* We add a second type of entity in the cluster */ 
	  new_entry = xcalloc(1, sizeof(struct color_entry));
	  init_hash_iset(&new_entry->entities);	  
	  insert_hash_iset(&new_entry->entities, e2.pos);
	  /* Similarity with 2 nodes */
	  new_entry->sim_entities = pairwise_sim(&new_entry->entities,
						 &ea->data[e2.type].ematrix);
	  hmap_add_int(&new->ce, &new_entry->entry, e2.type);
     }
     return new;
}

static inline int get_color(const struct node_ptr_array *node_color, int pos)
{
     int color;
   
     if (pos < 0)
	  fatal("Invalid position");
     if (node_color == NULL)
	  fatal("Invalid pointer to color");
     if((unsigned int)pos >= node_color->nr)
	  fatal("Error in color asignation"); 
     if (node_color->data[pos].cp == NULL)
	  color = NOCOLOR;
     else
	  color = node_color->data[pos].cp->id;

     return color;
}

static void update_saturation_degree(const struct graph_adj *g, pqueue_t *pq, int node,
				     const struct node_ptr_array *node_color)
{
     int r, color;
     struct arc_array adjs;
     unsigned int i, nr;
	  
     color = get_color(node_color, node);
     assert(color != NOCOLOR);
     adjs = get_adjacent_list(g, node);
     nr = adjs.nr;
     for (i = 0; i < nr; i++) {
	  r = increase_key(g, pq, adjs.data[i].to, color);
	  if (r == -1)
	       fatal("Error in update saturation degree\n");
     }
}

static int greatest_saturation_node(const struct graph_adj *g, pqueue_t *pq, 
				    const struct node_ptr_array *node_color)
{
     int r, color, node;
     saturation_t *ns;

     node = NONODE;
     ns = NULL;
     color = INFTY;
     r = extract_max(g, pq, &ns);

     if (r == -1)
	  fatal("No node without color");
     if (ns) {
	  color = get_color(node_color, ns->node);
	  node = ns->node;
     } else {
	  fatal("Error in get the greatest saturation node");
     }
     if (color != NOCOLOR)
	  fatal("Error in node to coloring");
#ifdef PRGDEBUG
     printf("Node %d; Num. of adjacent %d; Degree in %d\n", node, ns->n_adj_colors, g->degree[node].din);
#endif
     free_saturation_node(ns); 
     return node;
}

static bool *get_ady_used_color(const struct graph_adj *g,
				const struct node_ptr_array *node_color, int node)
{
     struct arc_array adjs;
     bool *color_used;
     int color;
     size_t alloc;
     unsigned int i, nr;
     
     alloc = g->n_nodes*sizeof(bool);
     color_used = xmalloc(alloc);
     memset(color_used, false, alloc);
     adjs = get_adjacent_list(g, node);
     nr = adjs.nr;
     for (i = 0; i < nr; i++) {
	  color = get_color(node_color, adjs.data[i].to);
	  assert(color < g->n_nodes);
	  if (color != NOCOLOR) 
	       color_used[color] = true;
     }
     return color_used;
}

static void get_free_colors(const struct graph_adj *g,
			    const struct node_ptr_array *solution,
			    int node, struct int_array *available_colors,
			    color_ptr_array_t *partitions)
{
     int cn;
     bool *color_used;
     struct color *ctmp;
     unsigned int i, n;
     
     n = partitions->nr;
     color_used = NULL;
     assert(g->n_nodes > node);
     assert(g->n_nodes >= (int)n+1); /* Number of colors not used is n+1 */
     color_used = get_ady_used_color(g, solution, node);
     cn = get_color(solution, node);
     if (cn != NOCOLOR) 
	  fatal("A adjacent node are using the same color");
     for (i = 0; i < n; i++) {
	  ctmp = partitions->data[i];
	  assert(ctmp->id == (int)i);
	  if (!color_used[i]) {
	       available_colors->data[available_colors->nr++] = i;
	  }
     }
     /* New color */
     if (available_colors->nr == 0)
	  available_colors->data[available_colors->nr++] = i;
     free(color_used);
}

static inline double avg_similarity_entities(struct color_entry *e)
{
     unsigned int n;
     double sim;
     
     n = e->entities.nr;
     if (n > 1) {
	  sim = e->sim_entities / ((n*n-n)/2.0);    
     } else if (n == 1) {
	  sim = e->sim_entities;
     } else {
	  if (eq_double(e->sim_entities, 0.0))
	       fatal("error in the similarity of the set of entities\n");
	  sim = e->sim_entities;
     }
     return sim;
}

static inline double get_sum_with_all_elements(const struct hash_iset *set,
					       unsigned int e,
					       const struct matrix *ematrix)
{
     unsigned int i, n, k, alloc;
     int *data;
     double total;

     k = 0;
     i = 0;
     n = set->nr;
     data = set->htable;
     alloc = set->alloc;
     total = 0.0;
     while ((i < n) && (k < alloc)) {
	  if (data[k] != NOELEM) {
	       total += similarity(ematrix, data[k], e);
	       i++;
	  }
	  k++;
     }
     assert(i == n);
     return total;
}

static inline sim_tuple_t agregate_similarity(double current_sim,
					      const struct hash_iset *set,
					      unsigned int e,
					      const struct matrix *ematrix)
{
     unsigned int n;
     double sim, total;
     sim_tuple_t t;
   
     sim = 0.0;
     total = 0.0;
     n = set->nr;
     if (lookup_hash_iset(set, e)) {
	  total = current_sim;
	  if (n == 1)
	       sim = current_sim;
	  else 
	       sim = current_sim / ((n*n-n)/2.0);
     } else {
	  total = get_sum_with_all_elements(set, e, ematrix);
	  /* If n == 1 the current_sim = similarity(a, a), then
	     the current_sim will be not added to  
	     partitions with n > 1 elements */  
	  if (n > 1)
	       total += current_sim;
	  n++;
	  sim = total / ((n*n-n)/2.0); 
     }
     t.total = total;
     t.sim = sim;

     return t;
}

static inline sim_tuple_t agregate_similarity_two_nodes(double current_sim,
							const struct hash_iset *set,
							unsigned int e1, unsigned int e2,
							const struct matrix *ematrix)
{
     unsigned int n, e;
     double sim, total;
     sim_tuple_t t;
     bool in_e1, in_e2;
   
     sim = 0.0;
     total = 0.0;
     n = set->nr;
     
     in_e1 = lookup_hash_iset(set, e1);
     in_e2 = lookup_hash_iset(set, e2);

     if (in_e1 && in_e2) {
	  total = current_sim;
	  if (n == 1) {
	       sim = current_sim;
	  } else {
	       sim = current_sim / ((n*n-n)/2.0);
	  }
     } else {
	  if ( (in_e1 && !in_e2) || (!in_e1 && in_e2) ) {
	       if (in_e1 && !in_e2)
		    e = e2;
	       else
		    e = e1;
	       total = get_sum_with_all_elements(set, e, ematrix);
	       /* If n == 1 the current_sim = similarity(a, a), then
		  the current_sim will be not added to  
		  partitions with n > 1 elements */  
	       if (n > 1) {
		    total += current_sim;
	       } 
	       n++;
	       sim = total / ((n*n-n)/2.0);
	  } else {
	       total = get_sum_with_all_elements(set, e1, ematrix) +
		    get_sum_with_all_elements(set, e2, ematrix) + current_sim;
	       if (in_e1 || in_e2)
		    fatal("Error in node search");
	       n += 2;
	       sim = total / ((n*n-n)/2.0);
	  }
     }
     t.total = total;
     t.sim = sim;
     return t;
}

static inline struct color_entry *get_color_entry(const struct hash_map *ce,
						  unsigned int key)
{
     struct hash_item *hentry;
     
     hentry = hmap_find_member_int(ce, key);
     if (hentry == NULL)
	  return NULL;
     return hash_entry(hentry, struct color_entry, entry);
}

static inline double sim_avg_type(double current_sim, const struct hash_iset *set)
{
     unsigned int n;
     double sim;
     
     n = set->nr;
     if (n == 1)
	  sim = current_sim;
     else 
	  sim = current_sim / ((n*n-n)/2.0);

     return sim;
}

static info_partition_t density_with_new_node(const color_ptr_array_t *partitions,
					      const struct node *new_node,
					      const struct entity_array *ea,
					      int color)
{

     struct color *cptr;
     double bt;
     info_partition_t ip;
     sim_tuple_t t1, t2;
     int type1, type2, n_new_types;
     struct color_entry *ctmp;
     double sim_type1, sim_type2, not_modf_cDensity, sim_type_avg1, sim_type_avg2;

     cptr = partitions->data[color];
     assert(cptr->id == color);
     type1 = new_node->e1.type;
     type2 = new_node->e2.type;
     n_new_types = 0;
     if (type1 == type2) {
	  ip.eq_type = true;
	  /* Case 1: Adjacent entities in the edge have the same type */
	  ctmp = get_color_entry(&cptr->ce, type1);
	  if (ctmp == NULL) {
	       sim_type1 = 0.0;
	       sim_type_avg1 = 0.0;
	       t1.total = similarity(&ea->data[type1].ematrix,
				     new_node->e1.pos,
				     new_node->e2.pos);
	       t1.sim =  t1.total;
	       n_new_types += 1;
	  } else {
	       sim_type1 = ctmp->sim_entities;
	       sim_type_avg1 = sim_avg_type(sim_type1, &ctmp->entities); 
	       t1 = agregate_similarity_two_nodes(sim_type1, &ctmp->entities,
						  new_node->e1.pos, new_node->e2.pos,
						  &ea->data[type1].ematrix);
	  }
	  t2.total = NOSIM;
	  t2.sim = NOSIM;
	  sim_type_avg2 = NOSIM;
     } else {
	  ip.eq_type = false;
	  /* Case 2: Adjacent entities in the edge has different types */
	  ctmp = get_color_entry(&cptr->ce, type1);
	  if (ctmp == NULL) {
	       sim_type1 = 0;
	       sim_type_avg1 = 0;
	       t1.total = similarity(&ea->data[type1].ematrix,
				     new_node->e1.pos,
				     new_node->e1.pos);
	       t1.sim = t1.total;
	       n_new_types += 1;
	  } else {
	       sim_type1 = ctmp->sim_entities;
	       sim_type_avg1 = sim_avg_type(sim_type1, &ctmp->entities);
	       t1 = agregate_similarity(sim_type1, &ctmp->entities,
					new_node->e1.pos,
					&ea->data[type1].ematrix);
	  }
	  ctmp = get_color_entry(&cptr->ce, type2);
	  if (ctmp == NULL) {
	       sim_type_avg2 = 0;
	       t2.total = similarity(&ea->data[type2].ematrix,
				     new_node->e2.pos,
				     new_node->e2.pos);
	       t2.sim = t2.total;
	       n_new_types += 1;
	  } else {
	       sim_type2 = ctmp->sim_entities;
	       sim_type_avg2 = sim_avg_type(sim_type2, &ctmp->entities);
	       t2 = agregate_similarity(sim_type2, &ctmp->entities,
					new_node->e2.pos,
					&ea->data[type2].ematrix);
	  }
     }
     ip.is_new_color = false;
     ip.color = color;
     ip.sim_entity_1.type = type1;
     ip.sim_entity_1.sim = t1.total; 
     ip.sim_entity_2.type = type2;
     ip.sim_entity_2.sim = t2.total;
     ip.sim_between = cptr->sim_between + new_node->sim;
     bt = (double)ip.sim_between/(cptr->id_nodes.nr + 1.0);
     
     if  ( cptr->id_nodes.nr > 1 ) {
	  if (type1 == type2) {
	       not_modf_cDensity = ((1.0 + cptr->ce.fill) * cptr->cDensity)
		    - sim_type_avg1 - cptr->sim_between/cptr->id_nodes.nr;

	       ip.cDensity = (double)((not_modf_cDensity + t1.sim + bt)  /
				      (n_new_types + 1.0 + cptr->ce.fill));
	  } else {
	       not_modf_cDensity = ((1.0 + cptr->ce.fill) * cptr->cDensity)
		    - sim_type_avg1 - sim_type_avg2 - cptr->sim_between/cptr->id_nodes.nr;

	       ip.cDensity = (double)((not_modf_cDensity + t1.sim + t2.sim + bt)  /
				      (n_new_types + 1.0 + cptr->ce.fill));
	  }
     } else {
	  if (type1 == type2) {
	       if (n_new_types > 0)
		    ip.cDensity = (t1.sim + bt + 1.0*n_new_types) /
			 (n_new_types + 1.0 + cptr->ce.fill);
	       else
		     ip.cDensity = (t1.sim + bt)/(n_new_types + 1.0 + cptr->ce.fill);
	  } else {
	       if (n_new_types > 0)
		    ip.cDensity = (t1.sim + t2.sim + bt + 1.0*n_new_types) /
			 (n_new_types + 1.0 + cptr->ce.fill);
	       else
		    ip.cDensity = (t1.sim + t2.sim + bt) /
			 (n_new_types + 1.0 + cptr->ce.fill);
	  }
     }
     if (ip.cDensity >= ONE)
	  fatal("Error in the computation of the density of a partition");
     
     return ip;
}

static info_partition_t get_best_color(const struct node_ptr_array *vn,
				       clusters_t *c, int new_node,
				       const struct int_array *free_colors,
				       const struct entity_array *ea)
{
     int i, n, n_colors, curr_color;
     double nc_best, nc_new, nc_current;
     struct node *nptr;
     info_partition_t ip_aux, ip_best;

     nc_new = 0;
     init_info_partition(&ip_best);
     n = free_colors->nr;
     nc_best = INFTY;
     n_colors = c->partitions.nr;
     nc_current = c->nc;
     nptr = &vn->data[new_node];
     for (i = 0; i < n; i++) {
	  curr_color = free_colors->data[i];
	  assert((curr_color >= 0) && (curr_color <= n_colors));
	  DEBUG("Color to evaluate %d ", curr_color);
	  if (n_colors == curr_color) {
	       /* We need to use a new color, then we have a partition with a element */
	       /* The density of a partition with a element is epsilon */
	       nc_new = nc_current + ((curr_color+1) - EPSILON);
	       ip_aux.is_new_color = true;
	       ip_aux.color = curr_color;
	       ip_aux.sim_between = nptr->sim;
	       ip_aux.cDensity = EPSILON;
	       ip_aux.nc = nc_new; 
	  } else {
	       /* We coloring the new node with a used color */
	       ip_aux = density_with_new_node(&c->partitions, nptr, ea, curr_color);
	       nc_new = nc_current + c->partitions.data[curr_color]->cDensity - ip_aux.cDensity;
	       ip_aux.nc = nc_new;
	  }
	  if (nc_best > nc_new) {
	       nc_best = nc_new;
	       ip_best = ip_aux;
	  }
     }
     return ip_best;
}

static void set_colors(clusters_t *c, info_partition_t *ip, struct node *nptr,
		       const struct entity_array *ea)
{
     struct color *cptr = NULL;
     struct color_entry *c_entry;
     
     if (ip->is_new_color) {
	  cptr = new_color(ip->color, ip->sim_between, nptr->id, nptr->e1, nptr->e2, ea);
	  cptr->cDensity = ip->cDensity;
	  ARRAY_PUSH(c->partitions, cptr);
     } else {
	  cptr = c->partitions.data[ip->color];
	  cptr->sim_between = ip->sim_between;
	  cptr->cDensity = ip->cDensity;
	  if (nptr->e1.type == nptr->e2.type) {
	       assert(ip->sim_entity_2.sim == NOSIM);
               /* Update entity of type 1 */
	       c_entry = get_color_entry(&cptr->ce, nptr->e1.type);
	       if (c_entry != NULL) {
		    c_entry->sim_entities = ip->sim_entity_1.sim;
		    insert_hash_iset(&c_entry->entities, nptr->e1.pos);
		    insert_hash_iset(&c_entry->entities, nptr->e2.pos);
	       } else {
		    c_entry = xcalloc(1, sizeof(struct color_entry));
		    c_entry->sim_entities =  ip->sim_entity_1.sim;
		    init_hash_iset(&c_entry->entities);
		    insert_hash_iset(&c_entry->entities, nptr->e1.pos);
		    insert_hash_iset(&c_entry->entities, nptr->e2.pos);
		    hmap_add_int(&cptr->ce, &c_entry->entry, nptr->e1.type);
	       }	       
	  } else {
               /* Update entity of type 1 */
	       c_entry = get_color_entry(&cptr->ce, nptr->e1.type);
	       if (c_entry != NULL) {
		    c_entry->sim_entities = ip->sim_entity_1.sim;
		    insert_hash_iset(&c_entry->entities, nptr->e1.pos);
	       } else {
		    c_entry = xcalloc(1, sizeof(struct color_entry));
		    c_entry->sim_entities =  ip->sim_entity_1.sim;
		    init_hash_iset(&c_entry->entities);
		    insert_hash_iset(&c_entry->entities, nptr->e1.pos);
		    hmap_add_int(&cptr->ce, &c_entry->entry, nptr->e1.type);
	       }
	       /* Update entity of type 2 */
	       c_entry = get_color_entry(&cptr->ce, nptr->e2.type);
	       if (c_entry != NULL) {
		    c_entry->sim_entities = ip->sim_entity_2.sim;
		    insert_hash_iset(&c_entry->entities, nptr->e2.pos);
	       } else  {
		    c_entry = xcalloc(1, sizeof(struct color_entry));
		    c_entry->sim_entities =  ip->sim_entity_2.sim;
		    init_hash_iset(&c_entry->entities);
		    insert_hash_iset(&c_entry->entities, nptr->e2.pos);
		    hmap_add_int(&cptr->ce, &c_entry->entry, nptr->e2.type);
	       }
	  }
	  ARRAY_PUSH(cptr->id_nodes, nptr->id);
     }
     nptr->cp = cptr;
     c->nc = ip->nc;
}

/**
 * Partiion graph solver based on the coloring algorithm called DSATUR
 */ 
static void coloring(const struct entity_array *entities, const struct graph_adj *g,
		     const struct node_ptr_array *nodes, clusters_t *c)
{
     struct color *cptr;
     struct node *nptr;
     int n, colored_nodes, new_node;
     pqueue_t pq_saturation;
     struct int_array free_colors;
     info_partition_t ip;

     init_struct_array(free_colors);
     colored_nodes = 0;
     ALLOC_GROW(free_colors.data, (unsigned int)g->n_nodes, free_colors.alloc);
     pq_init(&pq_saturation);
     init_saturation_pq(g, &pq_saturation);
     assert(c->partitions.nr == 0);
          
     /* We color a first node */
     new_node = greatest_saturation_node(g, &pq_saturation, nodes);

     if (new_node == NONODE)
	  fatal("Error getting the greatest saturation node");
     nptr = &nodes->data[new_node];
     assert(new_node == nptr->id);
     cptr = new_color(0, nptr->sim, nptr->id, nptr->e1, nptr->e2, entities);
     cptr->cDensity = EPSILON;
     ARRAY_PUSH(c->partitions, cptr);
     c->nc = 1.0 - cDensity(cptr, entities); /* First color */
     assert(cDensity(cptr, entities) ==
	    cDensity(c->partitions.data[c->partitions.nr-1], entities));
     nptr->cp = cptr; 
     colored_nodes++;
     if (pq_saturation.size != 0) {
	  update_saturation_degree(g, &pq_saturation, new_node, nodes);
     }

     /* We color all the nodes */
     n = g->n_nodes;
     while (colored_nodes < n) {
	  if ( (colored_nodes % COLOR_MARK) == 0)
	       printf("++++ Number of nodes colored so far: %d ++++\n", colored_nodes);
	       
	  new_node = greatest_saturation_node(g, &pq_saturation, nodes);

	  if (new_node == NONODE)
	       fatal("Error getting the greatest saturation node");
	  nptr = &nodes->data[new_node];
	  assert(new_node == nptr->id);
	  free_colors.nr = 0;
	  get_free_colors(g, nodes, new_node, &free_colors, &c->partitions);
#ifdef PRGDEBUG
	  printf("Free colors: ");
	  print_int_array(free_colors);
#endif
	  ip = get_best_color(nodes, c, new_node, &free_colors, entities);
	  if (ip.nc == NOSIM)
	       fatal("Best color is NULL");

	  set_colors(c, &ip, nptr, entities);
	  colored_nodes++;
	  if (pq_saturation.size != 0) {
	       update_saturation_degree(g, &pq_saturation, new_node, nodes);
	  }
#ifdef PRGDEBUG	  
	  if (!eq_double(color_density(c, entities), c->nc)) 
	       error("Different in the color density value %.3f %.3f\n",
		     color_density(c, entities), c->nc);
	  else
	       printf("Equal cDensity %.3f -- nc %.3f\n", color_density(c, entities), c->nc);
#endif
     }
     if (pq_saturation.size != 0)
	  fatal("Incomplete coloration\n");
     pq_delete(&pq_saturation);
     free(free_colors.data);
}

static double get_density_average(clusters_t *c)
{
     unsigned i, n;
     double r;

     r = 0.0;
     n = c->partitions.nr;
     for (i = 0; i < n; i++) {
	  r += c->partitions.data[i]->cDensity;
     }
     return r/n;
}

/*************************************
 *************************************
 **
 **  Get output files
 **
 ************************************
 ************************************/

static char *print_clustering(struct node_ptr_array *color_nodes,
				struct color_ptr_array *partitions,
				const struct entity_array *entities,
				const char *name)
{
     FILE *f;
     unsigned i, j, n, m;
     char *output1, *output2, *message, *s1, *s2;
     struct stat st;
     struct node edge;
     struct color *cluster;
     int id_node;

    if (asprintf(&output1, "%s-Clusters", name) == -1)
	  fatal("Error in output directory");
    
     if (stat(output1, &st) == -1)
	  mkdir(output1, 0700);
     else
	  fatal("The folder %s exists. Please, delete the folder and rerun SemEP", output1);

     if (asprintf(&message, "Cluster directory: %s\n", output1) == -1)
	  fatal("Error in output message");

     printf("Number of partitions: %u\n", partitions->nr);
     n = partitions->nr;
     for (i = 0; i < n; i++) {
	  cluster = partitions->data[i];
	  if (asprintf(&output2, "%s/cluster-%u.txt", output1, i) == -1)
	       fatal("Error in cluster file");
          f = fopen(output2, "w");
	  if (!f)
	       fatal("No descriptor file specified, abort\n");
	  m = cluster->id_nodes.nr;
	  for (j = 0; j < m; j++) {
	       id_node = cluster->id_nodes.data[j];
	       edge = color_nodes->data[id_node];
	       s1 = entities->data[edge.e1.type].vertices.data[edge.e1.pos];
	       s2 = entities->data[edge.e2.type].vertices.data[edge.e2.pos];
	       fprintf(f ,"%s\t%s\t%s\t%.4f\n", s1, s2, edge.relation, edge.sim);
	  }
	  fclose(f);
	  free(output2);
     }
     free(output1);
     return message;
}

/*************************************
 *************************************
 **
 **  Freeing memory
 **
 ************************************
 ************************************/

void free_entity_array(struct entity_array *e)
{
     unsigned int i, n;

     n = e->nr;
     for (i = 0; i < n; i++) {
	  free_string_array(&e->data[i].vertices);
	  free_double_matrix(e->data[i].ematrix.data, 0, 0);
     }
     free(e->data);
     e->nr = 0;
     e->alloc = 0;
}

void free_node_ptr_array(struct node_ptr_array *nptr)
{
     unsigned int i, n;

     n = nptr->nr;
     for (i = 0; i < n; i++) {
	  free(nptr->data[i].relation);
	  nptr->data[i].cp = NULL;
     }
     free(nptr->data);
     nptr->nr = 0;
     nptr->alloc = 0;
}

static void free_color(struct color *c)
{
     struct hash_item *hentry;
     struct hlist_node *n;
     struct color_entry *item;
    
     free(c->id_nodes.data);
     c->id_nodes.nr = 0;
     c->id_nodes.alloc = 0;
     hmap_for_each_safe_int(hentry, n, &c->ce) {
	  item = hash_entry(hentry, struct color_entry, entry);
	  if (item) {
	       free_hash_iset(&item->entities);
	       free(item);
	  }
     }
     hmap_destroy(&c->ce);
}

void free_clusters(clusters_t *c)
{
     unsigned int i, n;

     n = c->partitions.nr;
     for (i = 0; i < n; i++) {
	  free_color(c->partitions.data[i]);
	  free(c->partitions.data[i]);
     }
     free(c->partitions.data);
     c->partitions.nr = 0;
     c->partitions.alloc = 0;
}

/*************************************
 *************************************
 **
 **  semEP solver main function
 **
 ************************************
 ************************************/

clusters_t *semEP_solver(const struct entity_array *entities, const char *graph_name,
			 bool rel_constraint, struct node_ptr_array *color_nodes)
{
     clusters_t *c;
     double density;
     struct graph_adj gc;
     clock_t ti, tf;
     char *message;

     c = (clusters_t *)xcalloc(1, sizeof(clusters_t));
     init_clusters(c);     
     density = 0.0;
     ti = clock();
     init_graph_adj(&gc, color_nodes->nr);
     build_graph_to_coloring_matrix(entities, rel_constraint, &gc, color_nodes);
     tf = clock();
     printf("Time to build the graph to coloring: %.4f secs\n",
	    (double)(tf-ti)/CLOCKS_PER_SEC);
     printf("Graph to Coloring - Num. of Nodes: %d; Num. of Edges: %ld\n",
	    gc.n_nodes, gc.n_arcs/2);
#ifdef PRGDEBUG
     print_graph_adj(&gc);
#endif
     ti = clock();
     if (gc.n_nodes != 0) {
	  coloring(entities, &gc, color_nodes, c);
	  density = get_density_average(c);
#ifdef PRGDEBUG
	  printf("--- cDensity %.3f -- nc %.3f\n", color_density(c, entities), c->nc);
	  if (!eq_double(color_density(c, entities), c->nc)) 
	       fatal("Error in the color density value %.3f %.3f\n",
		     color_density(c, entities), c->nc);
	  print_coloring(c);
#endif	  
     } else {
	  fatal("Graph to coloring has no nodes");
     }
     tf = clock();
     printf("Coloring solver time: %.4f secs\n", (double)(tf-ti)/CLOCKS_PER_SEC);
     message = print_clustering(color_nodes, &c->partitions, entities, graph_name);
     printf("%s", message);
     free_graph_adj(&gc);
     free(message);
     printf("Average density of the partitions: %.4f \n", density);
     
     return c;
}
