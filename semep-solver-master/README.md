# SemEP 

## 1.  About

SemEP is an edge partitioning approach that combines
a data mining framework for link prediction, semantic knowledge
(similarities) from ontologies, and an algorithmic approach
to partition the edges of a heterogeneous graph.

For more information about SemEP see:

Guillermo Palma, Maria-Esther Vidal and Louiqa Raschid.
*Drug-Target Interaction Prediction Using Semantic Similarity and Edge
Partitioning*. Proceedings of the 12th International Semantic Web
Conference (ISWC 2013). Italy. 2014. [(PDF)](http://ldc.usb.ve/~gpalma/papers/semEP-ISWC14.pdf)

## 2. Content

* AUTHORS: list of contributors of the SemEP project.
* doc: Documentation about SemEP. 
* LICENSE: GPL version 2.
* Makefile: builds SemEP.
* README.md: this file.
* src: source code.
* datasets: datasets to test SemEP
** datasets/yamanishi: Drug-Target dataset proposed by Yamanishi [1].
** datasets/iDrug-Target: Drug-Target dataset proposed by Xuan Xiao et al. [2].
** datasets/pauwels: Drug side-effect dataset proposed by Pauwels et al. [3].
** datasets/datasets/biblio: artificial bibliographic graph.
** datasets/datasets/one_type_bp: artificial bibliographic graph with one type of vertices.
** datasets/datasets/testing: artificial graph for testing. 
* VERSION: SemEP version.

[1] K. Bleakley and Y. Yamanishi.
*Supervised prediction of drug-target interactions
using bipartite local models*.
Bioinformatics, 25(18):2397-2403, 2009.

[2] Xiao, X., Min, J.-L., Lin, W.-Z., Liu, Z., Cheng, X., and Chou, K.-C
"iDrug-Target: predicting the interactions between drug compounds 
and target proteins in cellular networking via benchmark dataset 
optimization approach"
Journal of Biomolecular Structure and Dynamics 33 (2015), 1–13.

[3] Pauwels, E., Stoven, V., and Yamanishi, Y.,
"Predicting drug side-effect profiles: a chemical fragment-based approach",
BMC Bioinformatics, 12:169, 2011.

## 3. License

GNU GENERAL PUBLIC LICENSE Version 2.

## 4. Requirements

* GNU Compiler Collection (GCC) or Clang.
* GNU make (make).

SemEP has been tested on FreeBSD, GNU/Linux, and OS X.

## 5. Installation

Clean and generate necessary files:

`$>make clean`

`$>make`

The result is that the executable file 'semEP' will be created.

## 6. Usage

SemEP has several mandatory command line arguments, and it has two optional command line arguments.
Mandatory means that without specifying this argument, the program will not work.

SemEP command synopsis:

`semEP [-c] [-p] -n <number of vertex type> <vertices type 1> <similarity matrix type 1> <threshold type 1> ... <vertices type n>  <similarity matrix type n> <threshold type n> <graph>`

where mandatory arguments are:
* <number of vertex types>: Number of vertex types 
* <vertices type 1> file with the list of vertices of type 1.
* <similarity matrix type 1> similarity matrix file with the similarities between the vertices of type 1.
* <threshold type 1>: threshold of similarity between the type 1 vertices of the graph. This value must be equal to or greater than zero.
* <vertices type n> file with the list of vertices of type n.
* <similarity matrix type 1> similarity matrix file with the similarities between the vertices of type n.
* <threshold type n>: threshold of similarity between the type n vertices of the graph. This value must be equal to or greater than zero
* <graph>: file with the bipartite graph to study.

where optional arguments are:

* [-c]: apply edge restriction. When you use this option, in the output clustering,
  two edges never are in the same cluster if they have different relationships.   
* [-p]: get the predicted links from the clusters. 

## 7. Running some samples

### 7.1 Computing the partitions of a drug-target bipartite graph on Yamanishi's enzyme dataset:

>./semEP -n 2 datasets/yamanishi/e/e_drugs.txt datasets/yamanishi/e/e_matrix_drugs.txt 0.2174 datasets/yamanishi/e/e_targets.txt datasets/yamanishi/e/e_matrix_targets.txt 0.0198 datasets/yamanishi/e/e_drug-target_graph.txt

### 7.2 Computing the partitions of a drug-target bipartite graph and get the predicted links on Yamanishi's nuclear receptor dataset:

>./semEP -p -n 2 datasets/yamanishi/nr/nr_drugs.txt datasets/yamanishi/nr/nr_matrix_drugs.txt 0.3061 datasets/yamanishi/nr/nr_targets.txt datasets/yamanishi/nr/nr_matrix_targets.txt 0.1614 datasets/yamanishi/nr/nr_drug-target_graph.txt

### 7.3 Computing the partitions using relationship constraints, of a drug-target bipartite graph on Xiao's enzyme dataset:

>./semEP -c -n 2 datasets/iDrug-Target/e/e_drugs.txt datasets/iDrug-Target/e/e_matrix_drugs.txt 0.25 datasets/iDrug-Target/e/e_targets.txt datasets/iDrug-Target/e/e_matrix_targets.txt 0.03 datasets/iDrug-Target/e/e_drug-target_graph.txt

### 7.4 Computing the partitions using relationship constraints, and get the predicted links of a drug-target bipartite graph on Xiao's ion channel dataset:

>./semEP -p -c -n 2 datasets/iDrug-Target/ic/ic_drugs.txt datasets/iDrug-Target/ic/ic_matrix_drugs.txt 0.4 datasets/iDrug-Target/ic/ic_targets.txt datasets/iDrug-Target/ic/ic_matrix_targets.txt 0.01 datasets/iDrug-Target/ic/ic_drug-target_graph.txt

### 7.5 Computing the partitions of a artificial bibliographic graph with 4 types of vertices, and get the predicted links:

./semEP -p -n 4 datasets/biblio/authors.txt datasets/biblio/matrix_authors.txt 0.4 datasets/biblio/conf.txt datasets/biblio/matrix_conf.txt 0.3 datasets/biblio/paper.txt datasets/biblio/matrix_paper.txt 0.1 datasets/biblio/term.txt datasets/biblio/matrix_term.txt 0.1 datasets/biblio/biblio_graph.txt

## 8 SemEP input

### 8.1. The file format of the graph

The input graph is an undirected graph. The first line of the graph file contains
the number of edges. The following lines in the file show the edges of the graph. 
Each line corresponds to an edge data and contains the first node,
the second node, the relationship type, and the cost of the edge.

	[number of edges n]
	[first-node edge-1][TAB][second-node edge-1][TAB][relationship-1][TAB][edge-cost-1]
	...
	...
	...
	[first-node edge-n][TAB][second-node edge-n][TAB][relationship-n][TAB][edge-cost-n]

### 8.2. The file format of the similarity matrix

Files with the matrices must contain the similarities between the vertices of the same type.
The file format is as follows:

	[number of rows and columns]
	[sim vertex-1 vertex-1][SPC]...[SPC][sim vertex-1 vertex-n]
	...
	...
	...
	[sim vertex-n vertex-1][SPC]...[SPC][sim vertex-n vertex-n]

### 8.3. The file format of the list vertices of the same type

The file format is as follows:

	[number of vertices n]
	[vertex 1]
	...
	...
	...
	[vertex n]

Where *[vertex x]* is the identifier of a vertex in the graph. The order of each vertex
corresponds to the position of the vertex in the corresponding similarity matrix.

## 9. SemEP output

When SemEP is executed without the option [-p], then SemEP produces as output
a directory with suffix *-Clusters* that contains the clusters with the
edge partitioning of the graph. Each cluster corresponds to a file on the directory.

Additionally, when SemEP is executed with the option [-p], SemEP creates a file with
suffix *-Predictions.txt* that contains the SemEP predictions for each cluster
and the probability of each prediction. 

## 10. Contact

I hope you find SemEP a useful tool. Please, let me know
any comment, problem, bug, or suggestion.

Guillermo Palma
[palma at l3s dot de](mailto:palma@l3s.de)
