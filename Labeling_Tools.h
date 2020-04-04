/*
 * Labeling_algorithms.h
 *
 *  Created on: 01 feb 2016
 *      Author: Mattia D'Emidio
 */
#include <omp.h>
#include <xmmintrin.h>
#include "CustomDataTypes.h"
#include "Labeling.h"
#include "boost/heap/fibonacci_heap.hpp"
#include <tuple>
#include <unistd.h>
#include <stdlib.h>

#include <vector>

#ifndef LABELING_TOOLS_H
#define LABELING_TOOLS_H

struct heap_data
{
	custom_node node;
	custom_weight weight;
	heap_data(custom_node n,custom_weight w)
	{
		node = n;
		weight = w;
	}

	bool operator<(heap_data const & rhs) const
	{
		return this->weight < rhs.weight;
	}
};
class Labeling_Tools {

private:

	NetworKit::Graph* graph;


	Labeling* index;


	//first element: order - second element: reverse_order
	std::pair<std::vector<custom_node>,std::vector<custom_node>> keeper;

	void check_closure(UpdateData*,custom_node,custom_node,custom_weight,bool);



	void clean_target(UpdateData*, custom_node, custom_node,std::vector<LabelEntry> &,custom_weight, bool);
	custom_node hasIndex(std::vector<LabelEntry> &, custom_node);

	custom_node x,y;
	custom_weight weight;
	NetworKit::GraphEvent::Type typ;
	custom_weight temp_weight;

	custom_node min_affected_x, max_affected_x, min_affected_y, max_affected_y;


	std::vector<bool> affected_vector_x;
	std::vector<bool> affected_vector_y;
	std::vector<custom_node> affected_x,affected_y;
	//we need to check whether a node u toward v (or from v) has an hub providing distance smaller than that going through bfs_i

	void setAffectedX(custom_node);
	void setAffectedY(custom_node);


	boost::heap::fibonacci_heap<heap_data>* Q;
	boost::heap::fibonacci_heap<heap_data>::handle_type* H;
	std::vector<short int> S;
	std::queue<custom_node> T;

	boost::heap::fibonacci_heap<heap_data>* prio_que;
	boost::heap::fibonacci_heap<heap_data>::handle_type* handles_prio_que;
	std::vector<short int> status_prio_que;
	std::queue<custom_node> trace_prio_que;


	std::vector<std::pair<custom_node,custom_weight>> queue;
	int que_h, que_t;
	std::vector<bool> marked;


	std::vector<custom_node> added_per_visit;

	std::vector<std::unordered_map<custom_node,custom_weight>> distances_of_affected_x_toward_affected_y;
	std::vector<std::unordered_map<custom_node,custom_weight>> distances_of_affected_y_from_affected_x;



	std::vector<custom_weight> root_label;




	void dynamicSetUp();

	void unweighted_build();
	void weighted_build();


	void decremental_algorithm(UpdateData*,bool);
	void incremental_algorithm(UpdateData*,bool);
	void aff_detection();


	int removal_of_outdated();
	void evaluation(UpdateData*);
	void former_restore(UpdateData*);
	int restore_cover(custom_node, std::vector<bool>&, bool, int);
	void restore(UpdateData*);

	void restoration(UpdateData*);

	bool dynamic;



	custom_node index_to_node(custom_node);
	custom_node node_to_index(custom_node);

	void resume_shortest_paths_over_affected(custom_node,custom_node,custom_weight,UpdateData*,bool,std::vector<bool>&);
	void evaluate(custom_node,custom_node,UpdateData*,bool);



	bool isIndexInSet(custom_node,std::vector<custom_node>&);

	//OPERATIONS ON STANDARD FIFO QUEUE
	std::pair<custom_node,custom_weight>* unweighted_get_min();
	bool unweighted_empty();
	void unweighted_pop();
	void unweighted_reset();
	void unweighted_insert(custom_node,custom_weight);


	//OPERATIONS ON AUXILIARY FIBONACCI HEAP
	void aux_insert(custom_node, custom_weight);
	void aux_dec_key(custom_node, custom_weight);
	void aux_erase(custom_node);
	custom_weight aux_min_dist();
	custom_node aux_min_node();
	custom_weight aux_key(custom_node);
	void aux_pop();
	bool aux_empty();
	void aux_reset();


	//OPERATIONS ON GENERIC FIBONACCI HEAP
	void weighted_insert(custom_node, custom_weight);
	void weighted_decrease_key(custom_node, custom_weight);
	custom_weight weighted_get_min_dist();
	custom_node weighted_get_min_node();
	custom_weight weighted_get_key(custom_node);
	void weighted_pop();
	bool weighted_empty();
	void weighted_reset();

	void affected_via_distances();

	void unweighted_trivial_affected(custom_node);
	void weighted_trivial_affected(custom_node);

	void unweighted_affected_via_distances(custom_node,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&);
	void weighted_affected_via_distances(custom_node,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&);
	void unweighted_affected_via_hubs(custom_node);
	void weighted_affected_via_hubs(custom_node);



	void printPaths(std::set<std::vector<custom_node> >&);
//	void clean_distance(std::vector<custom_weight>& , std::queue<custom_node> & );



	void custom_bfs(custom_node, std::vector<custom_weight>&, std::queue<custom_node>&, bool);
	void custom_dijkstra(custom_node, std::vector<custom_weight>&, std::queue<custom_node>&, bool);



//	void refresh(UpdateData*);

	void clean();


	void resume_shortest_paths(custom_node,custom_node,custom_weight,UpdateData*,bool,bool);


//	custom_node hasIndex(custom_node, custom_node, bool);
	void unweighted_relax(custom_node , custom_weight,bool);
	void weighted_relax(custom_node , custom_weight,bool);
	void unweighted_affected_relax(custom_node , custom_weight,bool,std::vector<bool>&);
	void weighted_affected_relax(custom_node , custom_weight,bool,std::vector<bool>&);

	bool test_affected_via_distances(custom_node, custom_node, std::vector<bool>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&,
			std::vector<custom_weight>&,std::vector<custom_weight>&);
	bool test_affected_via_hubs(custom_node, custom_node ,std::vector<bool>&,custom_weight);
	void getPaths(custom_node,custom_node, bool, std::set<std::vector<custom_node> >&);


	//funzione per creazione Labeling : ritorna l'ultimo nodo inserito nell'array che contiene l'ordinametno dei nodi
	custom_node lastNode();
	//funzione per aggiungere elemento a keeper 



public:

    void add_node_to_keeper(custom_node node, int index);
	std::vector<UpdateData> handle_affected_comparisons(NetworKit::GraphEvent*);
	double preprocessing_time;
	void update(NetworKit::GraphEvent*, UpdateData*,bool=true);
	Labeling_Tools(NetworKit::Graph*,Labeling*,std::pair<std::vector<custom_node>,std::vector<custom_node>>, bool);
    Labeling_Tools(NetworKit::Graph*,Labeling*,std::pair<std::vector<custom_node>,std::vector<custom_node>>);
	virtual ~Labeling_Tools();
    void weighted_build_RXL();

};

#endif /* LABELING_ALGORITHMS_H_ */
