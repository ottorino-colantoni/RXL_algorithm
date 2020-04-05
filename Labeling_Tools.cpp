/*
 * Labeling_algorithms.cpp
 *
 *  Created on: 01 feb 2016
 *      Author: Mattia D'Emidio
 */

#include "Labeling_Tools.h"
#include "mytimer.h"
#include "callgrind.h"


namespace {

// Return the number of threads that would be executed in parallel regions
int get_max_threads() {
#ifdef _OPENMP
	return omp_get_max_threads();
#else
	return 1;
#endif
}

int get_current_procs(){
#ifdef _OPENMP
	return omp_get_num_procs();
#else
	return 1;
#endif
}

// Set the number of threads that would be executed in parallel regions
void set_num_threads(int num_threads) {
#ifdef _OPENMP
	omp_set_num_threads(num_threads);
#else
	if (num_threads != 1) {
		throw std::runtime_error("compile with -fopenmp");
	}
#endif
}

// Return my thread ID
int get_thread_id() {
#ifdef _OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}

template<typename T>
struct parallel_vector {
	parallel_vector(size_t size_limit)
	: v(get_max_threads(), std::vector<T>(size_limit)),
	  n(get_max_threads(), 0) {}

	void push_back(const T &x) {
		int id = get_thread_id();
		v[id][n[id]++] = x;
	}

	void clear() {
		for(int i = 0;i<get_max_threads();i++)
			n[i] = 0;
	}

	std::vector<std::vector<T> > v;
	std::vector<size_t> n;
};
}  // namespace


custom_node Labeling_Tools::lastNode(){

	return this->index->order.back();
}





void Labeling_Tools::getPaths(custom_node n1,custom_node n2, bool forward_visit, std::set<std::vector<custom_node> >& paths){

	if(!forward_visit)
		std::swap(n1,n2);

	NetworKit::BFS* dij = new NetworKit::BFS(*graph,n1,true,false,n2);
	dij->run();
	paths = dij->getPaths(n2, true);
	assert(paths.size() == dij->getPaths(n2, true).size());
	delete dij;




}


void Labeling_Tools::printPaths(std::set<std::vector<custom_node> > & paths){

	int counter = 1;
	for(std::set<std::vector<custom_node> >::iterator it = paths.begin();it!=paths.end();it++){
		std::cout<<"#"<<counter<<" ";
		for(size_t t=0;t<(*it).size();t++){
			std::cout<<(*it)[t]<<" ";
		}
		std::cout<<"\n";
		counter++;
	}

}

Labeling_Tools::Labeling_Tools(NetworKit::Graph* g,Labeling* l,std::pair<std::vector<custom_node>,std::vector<custom_node>> order_keeper, bool dyn){

	CALLGRIND_STOP_INSTRUMENTATION;
	keeper = order_keeper;
	graph = g;
	index = l;
	preprocessing_time = std::numeric_limits<double>::max();
	dynamic = dyn ? true:false;


	boost::timer::auto_cpu_timer t;

	mytimer time_counter;
	time_counter.restart();

	INFO("Setting Up Data Structures");


	index->order = keeper.first;
	index->reverse_order = keeper.second;

	added_per_visit.resize(graph->upperNodeIdBound(),0);

	index->in_labels.clear();
	index->in_labels.resize(graph->upperNodeIdBound());

	if(graph->isDirected()){
		index->out_labels.clear();
		index->out_labels.resize(graph->upperNodeIdBound());
	}
	graph->parallelForNodes([&] (custom_node v){
		index->in_labels[v].push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));
		if(graph->isDirected())
			index->out_labels[v].push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));

	});


	if(!graph->isWeighted()){
		queue.resize(graph->upperNodeIdBound());
		marked.resize(graph->upperNodeIdBound(),false);
		que_h = 0, que_t = 0;
		unweighted_build();
	}
	else{
		prio_que = new boost::heap::fibonacci_heap<heap_data>();
		handles_prio_que = new boost::heap::fibonacci_heap<heap_data>::handle_type[graph->upperNodeIdBound()];
		status_prio_que.resize(graph->upperNodeIdBound(),0);
		weighted_build();
	}

	preprocessing_time =  time_counter.elapsed();
	//labeling_size = index->getLabelingSize();

	int sum = 0;

#pragma omp parallel for reduction(+:sum)
	for(int t=0;t<added_per_visit.size();++t)
		sum+=added_per_visit[t];
#pragma omp flush

	double total = 0;
	for(int t=0;;t++){
		total+=(double)added_per_visit[t]/(double)sum;
		if(total>=0.9){
			std::cout<<"90% of entries filled by "<<100*(double)t/(double)graph->numberOfNodes()<<" % of the visits "<<std::endl;
			break;
		}

	}
	if(dynamic)
		dynamicSetUp();

	assert(affected_x.empty() && affected_y.empty());

}


Labeling_Tools::Labeling_Tools(NetworKit::Graph* g,Labeling* l,std::pair<std::vector<custom_node>,std::vector<custom_node>> order_keeper){

    CALLGRIND_STOP_INSTRUMENTATION;
    keeper = order_keeper;
    graph = g;
    index = l;
    preprocessing_time = std::numeric_limits<double>::max();
    dynamic = false;

    boost::timer::auto_cpu_timer t;

    mytimer time_counter;
    time_counter.restart();

    INFO("Setting Up Data Structures");

    index->order = keeper.first;
    index->reverse_order = keeper.second;

    added_per_visit.resize(graph->upperNodeIdBound(),0);

    index->in_labels.clear();
    index->in_labels.resize(graph->upperNodeIdBound());

    graph->parallelForNodes([&] (custom_node v){
        index->in_labels[v].push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));
    });

    prio_que = new boost::heap::fibonacci_heap<heap_data>();
    handles_prio_que = new boost::heap::fibonacci_heap<heap_data>::handle_type[graph->upperNodeIdBound()];
    status_prio_que.resize(graph->upperNodeIdBound(),0);

}




void Labeling_Tools::unweighted_build() {

	ProgressStream builder_(graph->numberOfNodes());

	if(graph->isDirected())
		builder_.label() << "Building UNWEIGHTED DIRECTED labeling for " <<graph->numberOfNodes()<< " vertices";
	else
		builder_.label() << "Building UNWEIGHTED UNDIRECTED labeling for " <<graph->numberOfNodes()<< " vertices";



	std::vector<custom_time> crr_time(graph->upperNodeIdBound(), NULL_TIME), nxt_time(graph->upperNodeIdBound(), NULL_TIME);
	parallel_vector<std::pair<custom_node, LabelEntry> > pdiff_labels(graph->upperNodeIdBound()), pdiff_labels_bwd(graph->upperNodeIdBound());
	parallel_vector<custom_node> pdiff_nxt_que(graph->upperNodeIdBound()), pdiff_touched_ys(graph->upperNodeIdBound());

	for(custom_node order_of_source = 0; order_of_source<graph->upperNodeIdBound();order_of_source++){

		custom_node s = index_to_node(order_of_source);
		if(!graph->hasNode(s))
			continue;

		++builder_;


		//FORWARD
		std::vector<custom_node> crr_que;

		crr_que.push_back(s);
		crr_time[s] = 0;

		for (custom_weight explored_dist = 0; !crr_que.empty(); ++explored_dist) {

			size_t crr_que_size = crr_que.size();

#pragma omp parallel for schedule(guided, 1)

			for(size_t que_i = 0; que_i < crr_que_size; que_i++){

				custom_node v = crr_que[que_i]; // node in the queue, initially s
				custom_time t = crr_time[v];
				// Prune
				const std::vector<LabelEntry> &s1 = graph->isDirected() ? index->out_labels[s] : index->in_labels[s];
				const std::vector<LabelEntry> &s2 = index->in_labels[v];

				for (size_t i1 = 0, i2 = 0;;) {

					if (s1[i1].d > explored_dist && s1[i1].d < NULL_WEIGHT)	++i1;
					else if (s2[i2].d > explored_dist && s2[i2].d < NULL_WEIGHT) ++i2;
					else if (s1[i1].v < s2[i2].v) ++i1;
					else if (s1[i1].v > s2[i2].v) ++i2;
					//if none of the above conditions holds, we have reached a suitable (common) HUB that can be used to compute the current value of distance from s to v
					else{
						//pick vertex of entry i1 of label s1 (at s)
						custom_node v = s1[i1].v;
						assert(s1[i1].v == s2[i2].v);

						//if equal to NULL_NODE, end of label reached, break
						if (v == NULL_NODE)
							break;
						//sum of the distances given by the hub

						assert(s1[i1].d < NULL_WEIGHT && s2[i2].d < NULL_WEIGHT);
						custom_weight sum = (custom_weight)s1[i1].d + (custom_weight)s2[i2].d;
						if (sum <= explored_dist)
							goto prune; // pruning can be performed, by simply skipping the rest of the routine
						//otherwise, increase the two counters
						++i1; ++i2;
					}
				}

				assert(order_of_source<=node_to_index(v));
				pdiff_labels.push_back(std::make_pair(v, LabelEntry(order_of_source,explored_dist)));
				added_per_visit[order_of_source]++;


				graph->forNeighborsOf(v, [&](custom_node tv) {
					custom_time tt = std::max(t, EDGE_TIME_NULL);
					if (tt < crr_time[tv] && tt < nxt_time[tv]){
						for (;;) {
							custom_time prv_tt = nxt_time[tv];
							if (prv_tt <= tt)
								break;


							if (__sync_bool_compare_and_swap(&nxt_time[tv], prv_tt, tt)) {
								if (prv_tt == NULL_TIME) {
									pdiff_nxt_que.push_back(tv);
									if (crr_time[tv] == NULL_TIME)
										pdiff_touched_ys.push_back(tv);
								}
								break;
							}
						}
					}
				});


				prune:{}

			}
#pragma omp flush
			crr_que.clear();

			for(int i = 0;i<get_max_threads();i++){ // each thread does something
				for(size_t j = 0;j<pdiff_nxt_que.n[i];j++) {
					custom_node v = pdiff_nxt_que.v[i][j];
					crr_time[v] = nxt_time[v];
					nxt_time[v] = NULL_TIME;
					crr_que.push_back(v);
				}


			}
			pdiff_nxt_que.clear();

#pragma omp flush
		}

		for(int i = 0;i<get_max_threads();i++){
			for(size_t j = 0;j<pdiff_touched_ys.n[i];j++){
				custom_node v = pdiff_touched_ys.v[i][j];
				crr_time[v] = NULL_TIME;
			}
		}

		pdiff_touched_ys.clear();

		if(!graph->isDirected()){

			for(int i = 0;i<get_max_threads();i++){
				for(size_t j = 0;j<pdiff_labels.n[i];j++) {
					index->in_labels[pdiff_labels.v[i][j].first].back() = pdiff_labels.v[i][j].second;
					index->in_labels[pdiff_labels.v[i][j].first].push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));
				}
			}
			pdiff_labels.clear();
			continue;
		}



		//BACKWARD VISIT
		crr_que.push_back(s);
		crr_time[s] = 0;

		for (custom_weight explored_dist = 0; !crr_que.empty(); ++explored_dist) {

			size_t crr_que_size = crr_que.size();

#pragma omp parallel for schedule(guided, 1)

			for(size_t que_i = 0; que_i < crr_que_size;que_i++){

				custom_node v = crr_que[que_i]; // node in the queue, initially s
				custom_time t = crr_time[v];
				// Prune
				const std::vector<LabelEntry> &s1 = index->out_labels[v];
				const std::vector<LabelEntry> &s2 = index->in_labels[s];

				for (size_t i1 = 0, i2 = 0;;) {

					if (s1[i1].d > explored_dist && s1[i1].d < NULL_WEIGHT)	++i1;
					else if (s2[i2].d > explored_dist && s2[i2].d < NULL_WEIGHT) ++i2;
					else if (s1[i1].v < s2[i2].v) ++i1;
					else if (s1[i1].v > s2[i2].v) ++i2;

					else{
						assert(s1[i1].v == s2[i2].v);

						custom_node v = s1[i1].v;
						if (v == NULL_NODE)
							break;
						assert(s1[i1].d < NULL_WEIGHT && s2[i2].d < NULL_WEIGHT);
						custom_weight sum = (custom_weight)s1[i1].d + (custom_weight)s2[i2].d;
						if (sum <= explored_dist)
							goto prune_bwd; // pruning can be performed, by simply skipping the rest of the routine
						++i1; ++i2;
					}
				}
				assert(order_of_source<=node_to_index(v));
				pdiff_labels_bwd.push_back(std::make_pair(v, LabelEntry(order_of_source,explored_dist)));
				added_per_visit[order_of_source]++;
				graph->forInNeighborsOf(v, [&](custom_node tv) {

					custom_time tt = std::max(t, EDGE_TIME_NULL);
					if (tt < crr_time[tv] && tt < nxt_time[tv]){
						for (;;) {
							custom_time prv_tt = nxt_time[tv];
							if (prv_tt <= tt)
								break;

							if (__sync_bool_compare_and_swap(&nxt_time[tv], prv_tt, tt)) {
								if (prv_tt == NULL_TIME) {
									pdiff_nxt_que.push_back(tv);
									if (crr_time[tv] == NULL_TIME)	pdiff_touched_ys.push_back(tv);
								}
								break;
							}
						}
					}
				});


				prune_bwd:{}

			}
#pragma omp flush

			crr_que.clear();

			for(int i = 0;i<get_max_threads();i++){
				for(size_t j = 0;j<pdiff_nxt_que.n[i];j++) {
					custom_node v = pdiff_nxt_que.v[i][j];
					crr_time[v] = nxt_time[v];
					nxt_time[v] = NULL_TIME;
					crr_que.push_back(v);
				}
			}
			pdiff_nxt_que.clear();
#pragma omp flush

		}

		for(int i = 0;i<get_max_threads();i++){
			for(size_t j = 0;j<pdiff_touched_ys.n[i];j++){
				custom_node v = pdiff_touched_ys.v[i][j];
				crr_time[v] = NULL_TIME;
			}
		}
		pdiff_touched_ys.clear();



		for(int i = 0;i<get_max_threads();i++){
			for(size_t j = 0;j<pdiff_labels.n[i];j++) {
				index->in_labels[pdiff_labels.v[i][j].first].back() = pdiff_labels.v[i][j].second;
				index->in_labels[pdiff_labels.v[i][j].first].push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));
			}
			for(size_t j = 0;j<pdiff_labels_bwd.n[i];j++) {
				index->out_labels[pdiff_labels_bwd.v[i][j].first].back() = pdiff_labels_bwd.v[i][j].second;
				index->out_labels[pdiff_labels_bwd.v[i][j].first].push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));
			}
		}

		pdiff_labels.clear();
		pdiff_labels_bwd.clear();
	}
}

void Labeling_Tools::aux_insert(custom_node v,custom_weight d){
	assert(S[v]==0);
	H[v] = Q->emplace(heap_data(v,-d));
	S[v] = 1;
	T.push(v);
	assert((*H[v]).node==v);
	assert(aux_key(v)==d);
	//std::cout<<"INSERT\n";
}
bool Labeling_Tools::aux_empty(){
	return Q->empty();
}


custom_node Labeling_Tools::aux_min_node(){
	return Q->top().node;
}

custom_weight Labeling_Tools::aux_min_dist(){
	return -Q->top().weight;
}
custom_weight Labeling_Tools::aux_key(custom_node n){
	assert(S[n]>=1);
	return -(*H[n]).weight;
}
void Labeling_Tools::aux_pop(){
	Q->pop();
}



void Labeling_Tools::aux_erase(custom_node v){
#ifndef NDEBUG
	assert(S[v]==1);
	assert((*H[v]).node==v);
#endif
	Q->erase(H[v]);
	S[v]=0;

}


void Labeling_Tools::aux_dec_key(custom_node v,custom_weight d){
#ifndef NDEBUG
	assert(S[v]==1);
	assert(d<aux_key(v));
	assert((*H[v]).node==v);
#endif

	(*H[v]).weight = -d;
	Q->increase(H[v]);
	assert(aux_key(v)==d);
}
void Labeling_Tools::aux_reset(){
	while(!T.empty()) {
		assert(S[T.front()]>=0);
		S[T.front()] = 0;
		T.pop();
	}
#ifndef NDEBUG
	for(size_t t = 0;t<S.size();t++)
		assert(S[t]==0);
	assert(T.empty());
#endif
}

void Labeling_Tools::weighted_insert(custom_node v,custom_weight d){
	assert(status_prio_que[v]==0);
	handles_prio_que[v] = prio_que->emplace(heap_data(v,-d));
	status_prio_que[v] = 1;
	trace_prio_que.push(v);
	assert((*handles_prio_que[v]).node==v);
	assert(weighted_get_key(v)==d);
}
bool Labeling_Tools::weighted_empty(){
	return prio_que->empty();
}


custom_node Labeling_Tools::weighted_get_min_node(){
	return prio_que->top().node;
}
custom_weight Labeling_Tools::weighted_get_min_dist(){
	return -prio_que->top().weight;
}

custom_weight Labeling_Tools::weighted_get_key(custom_node n){
	assert(status_prio_que[n]>=1);
	return -(*handles_prio_que[n]).weight;
}
void Labeling_Tools::weighted_pop(){
	prio_que->pop();
}

void Labeling_Tools::weighted_decrease_key(custom_node v,custom_weight d){
#ifndef NDEBUG
	assert(status_prio_que[v]==1);
	assert(d<weighted_get_key(v));
	assert((*handles_prio_que[v]).node==v);

#endif
	(*handles_prio_que[v]).weight = -d;
	prio_que->increase(handles_prio_que[v]);
	assert(weighted_get_key(v)==d);
}
void Labeling_Tools::weighted_reset(){
	while(!trace_prio_que.empty()) {
		assert(status_prio_que[trace_prio_que.front()]>0);
		status_prio_que[trace_prio_que.front()] = 0;
		trace_prio_que.pop();
	}
	prio_que->clear();

#ifndef NDEBUG
	for(size_t t = 0;t<status_prio_que.size();t++)
		assert(status_prio_que[t]==0);
	assert(trace_prio_que.empty());
	assert(prio_que->empty());
#endif

}


void Labeling_Tools::unweighted_relax(custom_node current, custom_weight d, bool out){

	if(out){
		graph->forNeighborsOf(current, [&](custom_node n){
			if(!marked[n])
				unweighted_insert(n, d+1);

		});
	}
	else{
		graph->forInNeighborsOf(current, [&](custom_node n){
			if(!marked[n])
				unweighted_insert(n, d+1);

		});
	}
}

void Labeling_Tools::weighted_relax(custom_node current,custom_weight current_distance, bool out){
	if(out){
		graph->forNeighborsOf(current,[&](custom_node dn) {
			if(status_prio_que[dn]==0)
				weighted_insert(dn,current_distance + graph->weight(current, dn));
			else{
				if(status_prio_que[dn]==1 && weighted_get_key(dn) > current_distance + graph->weight(current,dn))
					weighted_decrease_key(dn,current_distance + graph->weight(current, dn));
			}
		});
	}
	else{
		graph->forInNeighborsOf(current,[&](custom_node up) {
			if(status_prio_que[up]==0)
				weighted_insert(up,current_distance + graph->weight(up, current));
			else{
				if(status_prio_que[up]==1 && weighted_get_key(up) > current_distance + graph->weight(up, current))
					weighted_decrease_key(up,current_distance + graph->weight(up, current));
			}
		});
	}

}


void Labeling_Tools::unweighted_affected_relax(custom_node current, custom_weight d, bool out,std::vector<bool>& A){

	if(out){
		graph->forNeighborsOf(current, [&](custom_node n){

#ifndef NDEBUG
			if(S[n]>=1)
				assert(A[node_to_index(n)]);
#endif
			//RELAX is restricted to affected nodes, that have not been already extracted from the main queue Q (if extracted, shortest path has been already found)
			//i.e. nodes having S[n]==2 have already been processed w.r.t. the root of the visit -- filter them out
			if(S[n]<=1 && A[node_to_index(n)] && !marked[n])
				unweighted_insert(n, d+1);

		});
	}
	else{
		graph->forInNeighborsOf(current, [&](custom_node n){
#ifndef NDEBUG
			if(S[n]>=1)
				assert(A[node_to_index(n)]);
#endif
			//RELAX is restricted to affected nodes, that have not been already extracted from the main queue Q (if extracted, shortest path has been already found)
			//i.e. nodes having S[n]==2 have already been processed w.r.t. the root of the visit -- filter them out
			if(S[n]<=1 && A[node_to_index(n)] && !marked[n])
				unweighted_insert(n, d+1);
		});
	}
}


void Labeling_Tools::weighted_affected_relax(custom_node current,custom_weight current_distance, bool out,std::vector<bool>& A){
	if(out){
		graph->forNeighborsOf(current,[&](custom_node dn) {
#ifndef NDEBUG
			if(S[dn]>=1)
				assert(A[node_to_index(dn)]);
#endif
			//RELAX is restricted to affected nodes, that have not been already extracted from the main queue Q (if extracted, shortest path has been already found)
			//i.e. nodes having S[dn]==2 have already been processed w.r.t. the root of the visit -- filter them out
			if(S[dn]<=1 && A[node_to_index(dn)]){
				if(status_prio_que[dn]==0)
					weighted_insert(dn,current_distance + graph->weight(current, dn));

				else{
					if(status_prio_que[dn]==1 && weighted_get_key(dn) > current_distance + graph->weight(current,dn)){
						weighted_decrease_key(dn,current_distance + graph->weight(current, dn));

					}
				}
			}
		});
	}
	else{

		graph->forInNeighborsOf(current,[&](custom_node up) {
#ifndef NDEBUG
			if(S[up]>=1)
				assert(A[node_to_index(up)]);
#endif
			//RELAX is restricted to affected nodes, that have not been already extracted from the main queue Q (if extracted, shortest path has been already found)
			//i.e. nodes having S[up]==2 have already been processed w.r.t. the root of the visit -- filter them out
			if(S[up]<=1 && A[node_to_index(up)]){
				if(status_prio_que[up]==0 )
					weighted_insert(up,current_distance + graph->weight(up, current));

				else{
					if(status_prio_que[up]==1 && weighted_get_key(up) > current_distance + graph->weight(up, current))
						weighted_decrease_key(up,current_distance + graph->weight(up, current));
				}

			}
		});
	}

}

void Labeling_Tools::custom_dijkstra(custom_node root,std::vector<custom_weight>& distance, std::queue<custom_node>& visited, bool forward){


	assert(weighted_empty());
	weighted_insert(root,0);
	distance[root] = 0;
	visited.push(root);

	custom_node current;
	custom_weight current_distance;


	while(!weighted_empty()){

		current = weighted_get_min_node();
		current_distance =  weighted_get_min_dist();
		status_prio_que[current]=2;
		weighted_pop();
		distance[current] = current_distance;
		visited.push(current);


		weighted_relax(current,current_distance,forward);
	}

	weighted_reset();
};


void Labeling_Tools::weighted_build(){

	ProgressStream builder_(graph->numberOfNodes());
	std::pair<custom_node,custom_weight> encoded = std::make_pair(NULL_NODE,NULL_WEIGHT);
	if(graph->isDirected())
		builder_.label() << "Building WEIGHTED DIRECTED labeling for " <<graph->numberOfNodes()<< " vertices";
	else
		builder_.label() << "Building WEIGHTED UNDIRECTED labeling for " <<graph->numberOfNodes()<< " vertices";

	for(custom_node order_of_source = 0; order_of_source<this->graph->upperNodeIdBound();order_of_source++){

		custom_node s =index_to_node(order_of_source);
		if(!this->graph->hasNode(s))
		 continue;

		++builder_;

#ifndef NDEBUG
		for(size_t t = 0;t<marked.size();t++)
			assert(status_prio_que[t]==0);
		assert(trace_prio_que.empty());
#endif

		//FORWARD VISIT
		assert(weighted_empty());
		weighted_insert(s,0);

		custom_node current;
		custom_weight current_distance;
		custom_weight labeling_distance;

		while(!weighted_empty()){
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			status_prio_que[current]=2;
			weighted_pop();

			if(node_to_index(current)<order_of_source)
				continue;
			index->query(s,current,current_distance,encoded);
			labeling_distance = encoded.second;

			if(labeling_distance<=current_distance)
				continue;

			index->in_labels[current].back().v = order_of_source;
			index->in_labels[current].back().d = current_distance;
			index->in_labels[current].push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));

			added_per_visit[order_of_source]++;

			weighted_relax(current,current_distance,true);

		}

		weighted_reset();

		if(!graph->isDirected())
			continue;

#ifndef NDEBUG
		for(size_t t = 0;t<marked.size();t++)
			assert(status_prio_que[t]==0);
		assert(trace_prio_que.empty());
#endif
		assert(weighted_empty());
		weighted_insert(s,0);

		while(!weighted_empty()){
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			status_prio_que[current]=2;
			weighted_pop();

			if(node_to_index(current)<order_of_source)
				continue;

			index->query(current,s,current_distance,encoded);
			labeling_distance = encoded.second;

			if(labeling_distance<=current_distance)
				continue;

			index->out_labels[current].back().v = order_of_source;
			index->out_labels[current].back().d = current_distance;
			index->out_labels[current].push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));

			added_per_visit[order_of_source]++;

			weighted_relax(current,current_distance,false);


		}
		weighted_reset();


	}
}





void Labeling_Tools::weighted_build_RXL(){

    ProgressStream builder_(graph->numberOfNodes());
    std::pair<custom_node,custom_weight> encoded = std::make_pair(NULL_NODE,NULL_WEIGHT);
    if(graph->isDirected())
        builder_.label() << "Building WEIGHTED DIRECTED labeling for " <<graph->numberOfNodes()<< " vertices";
    else
        builder_.label() << "Building WEIGHTED UNDIRECTED labeling for " <<graph->numberOfNodes()<< " vertices";


    custom_node s = lastNode();
    if(!this->graph->hasNode(s))
        return;	//continue;

    ++builder_;

#ifndef NDEBUG
    for(size_t t = 0;t<marked.size();t++)
        assert(status_prio_que[t]==0);
    assert(trace_prio_que.empty());
#endif

    //FORWARD VISIT
    assert(weighted_empty());
    weighted_insert(s,0);

    custom_node current;
    custom_weight current_distance;
    custom_weight labeling_distance;

    while(!weighted_empty()){
        current = weighted_get_min_node();
        current_distance =  weighted_get_min_dist();
        status_prio_que[current]=2;
        weighted_pop();

        //if(node_to_index(current)<order_of_source)
        //	continue;
		
        index->query(s,current,current_distance,encoded);
        labeling_distance = encoded.second;

        if(labeling_distance<=current_distance){
            continue;
	}
        index->in_labels[current].back().v = s;		//order_of_source;
        index->in_labels[current].back().d = current_distance;
        index->in_labels[current].push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));

        //added_per_visit[keeper.first.size()-1]++;

        weighted_relax(current,current_distance,true);

    }

    weighted_reset();
}



























void Labeling_Tools::unweighted_reset(){
	assert(!graph->isWeighted());

	for(custom_node i = 0;i<que_t;i++)
		marked[queue[i].first] = false;

	assert((double)que_h/(double)graph->numberOfNodes()<=1);
	que_h=0;
	que_t=0;
}

void Labeling_Tools::unweighted_insert(custom_node n, custom_weight d) {
	assert(!marked[n]);

	queue[que_t] = std::make_pair(n,d);
	que_t++;
	marked[n] = true;

}
std::pair<custom_node,custom_weight>* Labeling_Tools::unweighted_get_min(){
	return &queue[que_h];
}
void Labeling_Tools::unweighted_pop(){
	++que_h;
}

bool Labeling_Tools::unweighted_empty(){
	return que_h == que_t;
}

void Labeling_Tools::custom_bfs(custom_node root,std::vector<custom_weight>& distance, std::queue<custom_node>& visited, bool forward){


	assert((custom_node)queue.size() >= graph->upperNodeIdBound());
	assert(!marked[root]);
	assert(que_h==0 && que_t==0);

	unweighted_insert(root,0);

	std::pair<custom_node,custom_weight>* current_element;
	custom_node current;
	custom_weight current_distance;

	while (!unweighted_empty()) {
		current_element = unweighted_get_min();
		current = (*current_element).first;
		current_distance = (*current_element).second;
		unweighted_pop();

		distance[current] = current_distance;
		visited.push(current);


		unweighted_relax(current,current_distance,forward);


	}
	unweighted_reset();

};





void Labeling_Tools::dynamicSetUp(){


	INFO("Initializing data structures for dynamic algorithm");
	root_label.resize(graph->upperNodeIdBound(),NULL_WEIGHT);

	affected_vector_x.resize(graph->upperNodeIdBound(),false);
	affected_vector_y.resize(graph->upperNodeIdBound(),false);


	Q = new boost::heap::fibonacci_heap<heap_data>();
	H = new boost::heap::fibonacci_heap<heap_data>::handle_type[graph->upperNodeIdBound()];
	S.resize(graph->upperNodeIdBound(),0);
	assert(T.empty());

	INFO("Done!");

}


void Labeling_Tools::update(NetworKit::GraphEvent* e,UpdateData* data,bool new_vers){


	x = e->u;
	y = e->v;
	weight = e->w;
	temp_weight = graph->weight(x,y);
	typ = e->type;


	if(typ == NetworKit::GraphEvent::EDGE_REMOVAL || (typ == NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE && weight>0)){
		assert(graph->hasEdge(x,y));
		assert(weight==-1 || weight > graph->weight(x,y));
		CALLGRIND_START_INSTRUMENTATION;
		decremental_algorithm(data,new_vers);
		CALLGRIND_STOP_INSTRUMENTATION;
		std::cout<<"Removed: "<<data->removed<<" "<<" Added: "<<data->added<<"\n";
	}
	else{


#ifndef NDEBUG
		assert(typ == NetworKit::GraphEvent::EDGE_ADDITION || (typ == NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE && weight<0));
		if(typ == NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE){
			assert(-weight < graph->weight(x,y));
			assert(-weight > 0);
		}
		else
			assert(!graph->hasEdge(x,y));

#endif
		if(typ == NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE)
			weight=-weight;
		CALLGRIND_START_INSTRUMENTATION;
		incremental_algorithm(data,new_vers);
		CALLGRIND_STOP_INSTRUMENTATION;

		std::cout<<"Removed: "<<data->removed<<" "<<" Added: "<<data->added<<"\n";


	}


}

void Labeling_Tools::decremental_algorithm(UpdateData* data, bool new_vers) {

	mytimer time_counter;
	time_counter.restart();
	if(typ == NetworKit::GraphEvent::EDGE_REMOVAL)
		graph->removeEdge(x,y);
	else
		graph->setWeight(x,y,weight);

	aff_detection();




	data->affected_x = affected_x.size();
	data->affected_y = affected_y.size();
	data->affected = affected_x.size() + affected_y.size();
	double expon = (double)4/3;
	long long int tx = round(std::pow(graph->upperNodeIdBound(),expon));
	long long int val = (long long int)data->affected_x*(long long int)data->affected_y;

	std::cout<<"("<<affected_x.size()<<","<<affected_y.size()<<")\n";
	std::cout<<val<<" <=> "<<tx<<"\n";


	INFO("Updating");

	data->removed = removal_of_outdated();


	if(new_vers)
		restore(data);
	else
		former_restore(data);


	data->time =  time_counter.elapsed();


}



void Labeling_Tools::affected_via_distances(){

	INFO("Computing Distances");

	boost::timer::auto_cpu_timer t;
	std::vector<custom_weight> d_from_x, d_from_y, d_toward_x, d_toward_y, d_prime_toward_y, d_prime_from_x;
	std::queue<custom_node> visited_from_x, visited_from_y, visited_toward_x, visited_toward_y, visited_prime_toward_y, visited_prime_from_x;
	std::vector<custom_weight> d_x, d_y, d_prime_x, d_prime_y;
	std::queue<custom_node> visited_x, visited_y, visited_prime_x, visited_prime_y;

	if(graph->isDirected()){
		d_from_x.resize(graph->upperNodeIdBound(),NULL_WEIGHT);
		d_from_y.resize(graph->upperNodeIdBound(),NULL_WEIGHT);
		d_toward_x.resize(graph->upperNodeIdBound(),NULL_WEIGHT);
		d_toward_y.resize(graph->upperNodeIdBound(),NULL_WEIGHT);
		d_prime_from_x.resize(graph->upperNodeIdBound(),NULL_WEIGHT);
		d_prime_toward_y.resize(graph->upperNodeIdBound(),NULL_WEIGHT);

	}
	else{
		d_x.resize(graph->upperNodeIdBound(),NULL_WEIGHT);
		d_y.resize(graph->upperNodeIdBound(),NULL_WEIGHT);
		d_prime_x.resize(graph->upperNodeIdBound(),NULL_WEIGHT);
		d_prime_y.resize(graph->upperNodeIdBound(),NULL_WEIGHT);
	}
	if(graph->isDirected()){

		if(!graph->isWeighted()){
			custom_bfs(x,d_toward_x,visited_toward_x,false);
			custom_bfs(y,d_toward_y,visited_toward_y,false);
			custom_bfs(x,d_from_x,visited_from_x,true);
			custom_bfs(y,d_from_y,visited_from_y,true);
		}
		else{
			custom_dijkstra(x,d_toward_x,visited_toward_x,false);
			custom_dijkstra(y,d_toward_y,visited_toward_y,false);
			custom_dijkstra(x,d_from_x,visited_from_x,true);
			custom_dijkstra(y,d_from_y,visited_from_y,true);
		}
		if(typ == NetworKit::GraphEvent::EDGE_REMOVAL){
			graph->removeEdge(x,y);
		}
		else
			graph->setWeight(x,y,weight);
		if(!graph->isWeighted()){
			custom_bfs(x,d_prime_from_x,visited_prime_from_x,true);
			custom_bfs(y,d_prime_toward_y,visited_prime_toward_y,false);
		}
		else{
			custom_dijkstra(x,d_prime_from_x,visited_prime_from_x,true);
			custom_dijkstra(y,d_prime_toward_y,visited_prime_toward_y,false);
		}


	}
	else{

		if(!graph->isWeighted()){
			custom_bfs(x,d_x,visited_x,true);
			custom_bfs(y,d_y,visited_y,true);
		}
		else{
			custom_dijkstra(x,d_x,visited_x,true);
			custom_dijkstra(y,d_y,visited_y,true);
		}
		if(typ == NetworKit::GraphEvent::EDGE_REMOVAL){
			graph->removeEdge(x,y);
		}
		else
			graph->setWeight(x,y,weight);
		if(!graph->isWeighted()){
			custom_bfs(x,d_prime_x,visited_prime_x,true);
			custom_bfs(y,d_prime_y,visited_prime_y,true);
		}
		else{
			custom_dijkstra(x,d_prime_x,visited_prime_x,true);
			custom_dijkstra(y,d_prime_y,visited_prime_y,true);
		}

	}
	if(!graph->isWeighted()){
		unweighted_affected_via_distances(x,d_x,d_y,d_prime_x,d_prime_y,d_toward_x,d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y);
		unweighted_affected_via_distances(y,d_x,d_y,d_prime_x,d_prime_y,d_toward_x,d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y);
	}
	else{
		weighted_affected_via_distances(x,d_x,d_y,d_prime_x,d_prime_y,d_toward_x,d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y);
		weighted_affected_via_distances(y,d_x,d_y,d_prime_x,d_prime_y,d_toward_x,d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y);

	}



}
void Labeling_Tools::aff_detection(){


	min_affected_x = NULL_NODE;
	max_affected_x = 0;
	min_affected_y = NULL_NODE;
	max_affected_y = 0;

	INFO("Computation of Affected");


	if(!graph->isWeighted()){
		unweighted_affected_via_hubs(x);
		unweighted_affected_via_hubs(y);
	}
	else{
		weighted_affected_via_hubs(x);
		weighted_affected_via_hubs(y);
	}


}

void Labeling_Tools::unweighted_trivial_affected(custom_node node){

	boost::timer::auto_cpu_timer t;

	assert((custom_node)queue.size() >= graph->upperNodeIdBound());
	assert(!marked[node]);
	assert(que_h==0 && que_t==0);

	unweighted_insert(node,0);


	custom_node current, current_index;
	custom_weight dist;
	std::pair<custom_node,custom_weight>* current_element;

	if(node==x && graph->isDirected()){

		assert(affected_x.empty());
		assert(affected_y.empty());

		//forward_visit = false;
		while (!unweighted_empty()) {
			current_element = unweighted_get_min();
			current = (*current_element).first;
			current_index = node_to_index(current);
			dist = (*current_element).second;
			unweighted_pop();

			if(index->query(current,y)==dist+1){
				setAffectedX(current_index);
				unweighted_relax(current,dist,false);
			}
		}

	}
	else if(node==x && !graph->isDirected()){
		assert(affected_x.empty());
		assert(affected_y.empty());

		//forward_visit = false;
		while (!unweighted_empty()) {
			current_element = unweighted_get_min();
			current = (*current_element).first;
			current_index = node_to_index(current);
			dist = (*current_element).second;
			unweighted_pop();

			if(index->query(current,y)==dist+1){
				setAffectedX(current_index);
				unweighted_relax(current,dist,true);

			}

		}

	}
	else{
		assert(node==y);
		assert(!affected_x.empty());

		while (!unweighted_empty()) {

			current_element = unweighted_get_min();
			current = (*current_element).first;
			current_index = node_to_index(current);
			dist = (*current_element).second;
			unweighted_pop();

			if(index->query(x,current)==dist+1){
				setAffectedY(current_index);
				unweighted_relax(current,dist,true);

			}


		}

	}
	unweighted_reset();


}


void Labeling_Tools::weighted_trivial_affected(custom_node node){

	boost::timer::auto_cpu_timer t;
	assert(weighted_empty());
	weighted_insert(node,0);

	custom_node current, current_index;
	custom_weight current_distance;



	if(node==x && graph->isDirected()){

		assert(affected_x.empty());
		assert(affected_y.empty());


		//forward_visit = false;
		while (!weighted_empty()) {
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			current_index = node_to_index(current);
			status_prio_que[current]=2;
			weighted_pop();

			if(index->query(current,y)==current_distance+temp_weight){
				setAffectedX(current_index);
				weighted_relax(current,current_distance,false);

			}

		}

	}
	else if(node==x && !graph->isDirected()){
		assert(affected_x.empty());
		assert(affected_y.empty());

		while (!weighted_empty()) {
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			current_index = node_to_index(current);
			status_prio_que[current]=2;
			weighted_pop();

			if(index->query(current,y)==current_distance+temp_weight){
				setAffectedX(current_index);
				weighted_relax(current,current_distance,true);
			}

		}

	}
	else{
		assert(node==y);
		assert(affected_x.size()>=0);

		while (!weighted_empty()) {
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			current_index = node_to_index(current);
			status_prio_que[current]=2;
			weighted_pop();
			if(index->query(x,current)==current_distance+temp_weight){
				setAffectedY(current_index);
				weighted_relax(current,current_distance,true);

			}


		}

	}
	weighted_reset();


}

void Labeling_Tools::unweighted_affected_via_hubs(custom_node node){

	boost::timer::auto_cpu_timer t;

	assert((custom_node)queue.size() >= graph->upperNodeIdBound());
	assert(!marked[node]);
	assert(que_h==0 && que_t==0);

	unweighted_insert(node,0);

	custom_node current, current_index;
	custom_weight dist;
	std::pair<custom_node,custom_weight>* current_element;

	if(node==x && graph->isDirected()){

		assert(affected_x.empty());
		assert(affected_y.empty());

		while (!unweighted_empty()) {

			current_element = unweighted_get_min();
			current = (*current_element).first;
			current_index = node_to_index(current);
			dist = (*current_element).second;
			unweighted_pop();

			if(!test_affected_via_hubs(current,y,affected_vector_x,dist+1))
				continue;

			setAffectedX(current_index);
			unweighted_relax(current,dist,false);



		}

	}
	else if(node==x && !graph->isDirected()){
		assert(affected_x.empty());
		assert(affected_y.empty());


		//forward_visit = false;
		while (!unweighted_empty()) {

			current_element = unweighted_get_min();
			current = (*current_element).first;
			current_index = node_to_index(current);
			dist = (*current_element).second;
			unweighted_pop();
			if(!test_affected_via_hubs(current,y,affected_vector_x,dist+1))
				continue;

			setAffectedX(current_index);

			unweighted_relax(current,dist,true);


		}

	}

	else{
		assert(node==y);
		assert(!affected_x.empty());
		//		neighbor = x;

		//		forward_visit = true;
		while (!unweighted_empty()) {

			current_element = unweighted_get_min();
			current = (*current_element).first;
			current_index = node_to_index(current);
			dist = (*current_element).second;
			unweighted_pop();

			if(!test_affected_via_hubs(x,current,affected_vector_y,dist+1))
				continue;

			setAffectedY(current_index);
			unweighted_relax(current,dist,true);



		}

	}



	unweighted_reset();

}


void Labeling_Tools::weighted_affected_via_hubs(custom_node node){



	boost::timer::auto_cpu_timer t;

	assert(weighted_empty());
	weighted_insert(node,0);

	custom_node current, current_index;
	custom_weight current_distance;

	if(node==x && graph->isDirected()){

		assert(affected_x.empty());
		assert(affected_y.empty());

		//forward_visit = false;
		while (!weighted_empty()) {
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			current_index = node_to_index(current);
			status_prio_que[current]=2;
			weighted_pop();

			if(!test_affected_via_hubs(current,y,affected_vector_x,current_distance+temp_weight))
				continue;

			setAffectedX(current_index);
			weighted_relax(current,current_distance,false);
		}

	}
	else if(node==x && !graph->isDirected()){
		assert(affected_x.empty());
		assert(affected_y.empty());
		while (!weighted_empty()) {
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			current_index = node_to_index(current);
			status_prio_que[current]=2;
			weighted_pop();

			if(!test_affected_via_hubs(current,y,affected_vector_x,current_distance+temp_weight))
				continue;

			setAffectedX(current_index);
			weighted_relax(current,current_distance,true);


		}

	}

	else{
		assert(node==y);
		assert(affected_x.size()>=0);

		while (!weighted_empty()) {
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			current_index = node_to_index(current);
			status_prio_que[current]=2;
			weighted_pop();

			if(!test_affected_via_hubs(x,current,affected_vector_y,current_distance+temp_weight))
				continue;

			setAffectedY(current_index);
			weighted_relax(current,current_distance,true);

		}

	}

	weighted_reset();

}



void Labeling_Tools::unweighted_affected_via_distances(custom_node node,
		std::vector<custom_weight>& d_x,
		std::vector<custom_weight>& d_y,
		std::vector<custom_weight>& d_prime_x,
		std::vector<custom_weight>& d_prime_y,
		std::vector<custom_weight>& d_toward_x,
		std::vector<custom_weight>& d_toward_y,
		std::vector<custom_weight>& d_from_x,
		std::vector<custom_weight>& d_from_y,

		std::vector<custom_weight>& d_prime_from_x,
		std::vector<custom_weight>& d_prime_toward_y
){

	//forward refers to the direction of the visit
	boost::timer::auto_cpu_timer t;
	custom_node current, current_index;
	std::pair<custom_node,custom_weight>* current_element;

	assert((custom_node)queue.size() >= graph->upperNodeIdBound());
	assert(!marked[node]);
	assert(que_h==0 && que_t==0);

	unweighted_insert(node,0);

	if(node==x && graph->isDirected()){

		assert(affected_x.empty());
		assert(affected_y.empty());
		while (!unweighted_empty()) {

			current_element = unweighted_get_min();
			current = (*current_element).first;
			current_index = node_to_index(current);
			unweighted_pop();


			if(d_prime_toward_y[current] != d_toward_y[current] || (d_toward_y[current] == d_toward_x[current] + 1 &&
					test_affected_via_distances(current,y,affected_vector_x,
							d_x,d_y,d_prime_x,d_prime_y,d_toward_x,
							d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y))){
#ifndef NDEBUG
				if(d_toward_y[current] != d_toward_x[current] + 1){
					assert(false);
				}
#endif
				setAffectedX(current_index);
				unweighted_relax(current,0,false);


			}
		}
	}
	else if(node==x && !graph->isDirected()){
		assert(affected_x.empty());
		assert(affected_y.empty());

		while (!unweighted_empty()) {

			current_element = unweighted_get_min();
			current = (*current_element).first;
			current_index = node_to_index(current);
			unweighted_pop();

			if(d_y[current] != d_prime_y[current] || (d_y[current] == d_x[current] + 1 &&
					test_affected_via_distances(current,y,affected_vector_x,d_x,d_y,d_prime_x,d_prime_y,d_toward_x,
							d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y))){
#ifndef NDEBUG
				if(d_y[current] != d_x[current] + 1){
					assert(false);
				}
#endif
				setAffectedX(current_index);
				unweighted_relax(current,0,true);


			}


		}


	}
	else if(node==y && !graph->isDirected()){
		assert(node==y);
		assert(!affected_x.empty());

		while (!unweighted_empty()) {

			current_element = unweighted_get_min();
			current = (*current_element).first;
			current_index = node_to_index(current);
			unweighted_pop();

			if(d_x[current] != d_prime_x[current] || (d_x[current] == d_y[current] + 1 &&
					test_affected_via_distances(x,current,affected_vector_y,d_x,d_y,d_prime_x,d_prime_y,d_toward_x,
							d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y))){

#ifndef NDEBUG
				if(d_x[current] != d_y[current] + 1){
					assert(false);
				}
#endif

				setAffectedY(current_index);

				unweighted_relax(current,0,true);


			}


		}
	}
	else{
		assert(node==y);
		assert(!affected_x.empty());
		while (!unweighted_empty()) {

			current_element = unweighted_get_min();
			current = (*current_element).first;
			current_index = node_to_index(current);
			unweighted_pop();

			if(d_prime_from_x[current] != d_from_x[current] || (d_from_x[current] == d_from_y[current] + 1 &&
					test_affected_via_distances(x,current,affected_vector_y,d_x,d_y,d_prime_x,d_prime_y,d_toward_x,
							d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y)


			)){

#ifndef NDEBUG
				if(d_from_x[current] != d_from_y[current] + 1){
					assert(false);
				}
#endif

				setAffectedY(current_index);

				unweighted_relax(current,0,true);

			}


		}
	}

	unweighted_reset();

}



void Labeling_Tools::weighted_affected_via_distances(
		custom_node node,
		std::vector<custom_weight>& d_x,
		std::vector<custom_weight>& d_y,
		std::vector<custom_weight>& d_prime_x,
		std::vector<custom_weight>& d_prime_y,
		std::vector<custom_weight>& d_toward_x,
		std::vector<custom_weight>& d_toward_y,
		std::vector<custom_weight>& d_from_x,
		std::vector<custom_weight>& d_from_y,
		std::vector<custom_weight>& d_prime_from_x,
		std::vector<custom_weight>& d_prime_toward_y) {


	boost::timer::auto_cpu_timer t;
	assert(weighted_empty());
	weighted_insert(node,0);

	custom_node current, current_index;
	custom_weight current_distance;

	if(node==x && graph->isDirected()){

		assert(affected_x.empty());
		assert(affected_y.empty());
		//forward_visit = false;
		while (!weighted_empty()) {
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			current_index = node_to_index(current);
			status_prio_que[current]=2;
			weighted_pop();

			if(d_prime_toward_y[current] != d_toward_y[current] ||
					(d_toward_y[current] == d_toward_x[current] + temp_weight &&
							test_affected_via_distances(current,y,affected_vector_x,d_x,d_y,d_prime_x,d_prime_y,d_toward_x,
									d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y))){
#ifndef NDEBUG
				if(d_toward_y[current] != d_toward_x[current] + temp_weight){
					assert(false);
				}
#endif
				setAffectedX(current_index);


				weighted_relax(current,current_distance,false);


			}
		}
	}
	else if(node==x && !graph->isDirected()){
		assert(affected_x.empty());
		assert(affected_y.empty());
		while (!weighted_empty()) {
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			current_index = node_to_index(current);
			status_prio_que[current]=2;
			weighted_pop();

			if(d_y[current] != d_prime_y[current] || (d_y[current] == d_x[current] + temp_weight
					&& test_affected_via_distances(current,y,affected_vector_x,d_x,d_y,d_prime_x,d_prime_y,d_toward_x,
							d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y))){
#ifndef NDEBUG
				if(d_y[current] != d_x[current] + temp_weight){
					assert(false);
				}
#endif
				setAffectedX(current_index);


				weighted_relax(current,current_distance,true);
			}

		}


	}
	else if(node==y && !graph->isDirected()){
		assert(node==y);
		assert(affected_x.size()>=0);
		while (!weighted_empty()) {
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			current_index = node_to_index(current);
			status_prio_que[current]=2;
			weighted_pop();

			if(d_x[current] != d_prime_x[current] || (d_x[current] == d_y[current] + temp_weight &&
					test_affected_via_distances(x,current,affected_vector_y,d_x,d_y,d_prime_x,d_prime_y,d_toward_x,
							d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y))){

#ifndef NDEBUG
				if(d_x[current] != d_y[current] + temp_weight){
					assert(false);
				}
#endif
				setAffectedY(current_index);
				weighted_relax(current,current_distance,true);
			}


		}
	}
	else{
		assert(node==y);
		assert(affected_x.size()>=0);
		while (!weighted_empty()) {
			current = weighted_get_min_node();
			current_distance =  weighted_get_min_dist();
			current_index = node_to_index(current);
			status_prio_que[current]=2;
			weighted_pop();

			if(d_prime_from_x[current] != d_from_x[current] || (d_from_x[current] == d_from_y[current] + temp_weight &&
					test_affected_via_distances(x,current,affected_vector_y,d_x,d_y,d_prime_x,d_prime_y,d_toward_x,
							d_toward_y,d_from_x,d_from_y,d_prime_from_x,d_prime_toward_y))){

#ifndef NDEBUG
				if(d_from_x[current] != d_from_y[current] + temp_weight){
					assert(false);
				}
#endif

				setAffectedY(current_index);

				weighted_relax(current,current_distance,true);
			}


		}
	}

	weighted_reset();

}

void Labeling_Tools::incremental_algorithm(UpdateData* data, bool minimal) {
	boost::timer::auto_cpu_timer t;

	mytimer time_counter;
	time_counter.restart();

	data->added = 0;
	data->removed = 0;

	std::cout<<"Updating labeling... ";

	assert(weight>0);
	min_affected_x = NULL_NODE;
	max_affected_x = 0;
	min_affected_y = NULL_NODE;
	max_affected_y = 0;


	if(typ == NetworKit::GraphEvent::EDGE_ADDITION)
		graph->addEdge(x,y,weight);
	else
		graph->setWeight(x,y,weight);

	if(!graph->isDirected()){

		std::vector<LabelEntry> &idx_x = index->in_labels[x];
		std::vector<LabelEntry> &idx_y = index->in_labels[y];

		int index_x = 0, index_y = 0, previous_vertex = -1;

		for (;;) {

			LabelEntry l_x = idx_x[index_x];
			LabelEntry l_y = idx_y[index_y];

			if (index_x > 0 && (int)l_x.v == previous_vertex)
				++index_x;
			else if (index_y > 0 && (int)l_y.v == previous_vertex)
				++index_y;
			else if (l_x.v < l_y.v){
				// u -> v
				resume_shortest_paths(previous_vertex = l_x.v, y, l_x.d + graph->weight(x, y),data,true,minimal);
				++index_x;

			}
			else if (l_x.v > l_y.v){  // v -> u
				resume_shortest_paths(previous_vertex = l_y.v, x, l_y.d + graph->weight(x, y),data,true,minimal);
				++index_y;
			}
			else{
				// u <-> v - both
				assert(l_x.v==l_y.v);
				if (l_x.v == NULL_NODE)
					break;

				if (l_x.d + graph->weight(x, y) < l_y.d)
					resume_shortest_paths(l_x.v, y, l_x.d + graph->weight(x, y),data,true,minimal);


				if (l_y.d + graph->weight(x, y) < l_x.d)
					resume_shortest_paths(l_y.v, x, l_y.d + graph->weight(x, y),data,true,minimal);


				++index_x;
				++index_y;
				previous_vertex = l_x.v;

			}
		}

	}
	else{


		std::vector<LabelEntry> &idx_in_x = index->in_labels[x];
		std::vector<LabelEntry> &idx_in_y = index->in_labels[y];
		std::vector<LabelEntry> &idx_out_x = index->out_labels[x];
		std::vector<LabelEntry> &idx_out_y = index->out_labels[y];


		int index_in_x = 0 , index_in_y = 0;
		int prv_in = -1;

		int index_out_x = 0, index_out_y = 0;
		int prv_out = -1;


		for (;;) {

			LabelEntry l_in_x = idx_in_x[index_in_x];
			LabelEntry l_in_y = idx_in_y[index_in_y];
			LabelEntry l_out_x = idx_out_x[index_out_x];
			LabelEntry l_out_y = idx_out_y[index_out_y];

			bool min_side = std::min(l_in_x.v,l_in_y.v)<=std::min(l_out_y.v,l_out_x.v) ? true : false;
			if(min_side){
				if (l_in_x.v == l_in_y.v){
					if (l_in_y.v == NULL_NODE)
						break;
					assert(l_in_y.v == l_in_x.v);

					if (l_in_x.d + graph->weight(x, y) < l_in_y.d)
						resume_shortest_paths(prv_in = l_in_x.v, y, l_in_x.d + graph->weight(x, y), data,true,minimal);

					++index_in_x;
					++index_in_y;
				}
				else if (l_in_x.v < l_in_y.v){
					if(index_in_x > 0 && (int)l_in_x.v == prv_in)
						++index_in_x;
					else{
						resume_shortest_paths(prv_in = l_in_x.v, y, l_in_x.d + graph->weight(x, y), data,true,minimal);
						++index_in_x;
					}
				}
				else{
					assert(l_in_x.v > l_in_y.v);
					++index_in_y;
				}
				continue;
			}
			else{
				if (l_out_y.v == l_out_x.v){
					if (l_out_y.v == NULL_NODE)
						break;

					assert(l_out_y.v == l_out_x.v);
					if (l_out_y.d + graph->weight(x, y) < l_out_x.d)
						resume_shortest_paths(prv_out = l_out_x.v, x, l_out_y.d + graph->weight(x, y), data,false,minimal);

					++index_out_x;
					++index_out_y;
				}
				else if(l_out_y.v < l_out_x.v){
					if(index_out_y > 0 && (int)l_out_y.v == prv_out)
						++index_out_y;

					else{
						resume_shortest_paths(prv_out = l_out_y.v, x, l_out_y.d + graph->weight(x, y), data,false,minimal);
						++index_out_y;
					}
				}
				else{
					assert(l_out_y.v > l_out_x.v);
					++index_out_x;
				}
				continue;
			}

		}


	}


	std::cout<<"done!\n"<<std::flush;
	data->affected_x = std::count(affected_vector_x.begin(),affected_vector_x.end(),true);
	data->affected_y = std::count(affected_vector_y.begin(),affected_vector_y.end(),true);
	data->affected = data->affected_x+data->affected_y;

	data->time = time_counter.elapsed();
}
//
//void Labeling_Tools::post(custom_node n, UpdateData* data, bool x_side){
//
//	std::vector<LabelEntry>& label = graph->isDirected() && x_side ? index->out_labels[n] : index->in_labels[n];
//
//
//	std::vector<bool>& AFF_mine = x_side ? affected_vector_x : affected_vector_y;
//	std::vector<bool>& AFF_others = x_side ? affected_vector_y : affected_vector_x;
//
//	for(int t = 0;;++t) {
//
//
//		custom_node ind_x = label[t].v;
//
//		if(ind_x==NULL_NODE)
//			break;
//
//		assert(ind_x!=NULL_NODE);
//
//		if(AFF_mine[ind_x])
//			continue;
//
//
//
//
//		custom_weight d = label[t].d;
//
////#ifndef NDEBUG
////		custom_weight comparison_d;
////		std::vector<custom_node> hubset = x_side ? index->hubs(n,index_to_node(ind_x),comparison_d) :
////						index->hubs(index_to_node(ind_x),n, comparison_d);
////		for(auto&el:hubset)
////			assert(el<=ind_x && el<=node_to_index(n));
////		//no hub can be larger than the two endpoints indices ind_x and node_to_index(n)
////
////#endif
//		//IF DIRECTED AND X_SIDE then label is out_label (query from n to index_to_node(ind_x))
//		//HENCE hub_label must be set to in_labels
//
//		//IF DIRECTED AND !X_SIDE then label is in_label(query from index_to_node(ind_x) to n)
//		//HENCE hub_label must be set to out_labels
//		//IF UNDIRECTED -  NEVERMIND - set to in_labels;
//
//		size_t index_n = 0;
//		size_t index_h = 0;
//
//		std::vector<LabelEntry> & target_label = (!graph->isDirected() || x_side)
//				? index->in_labels[index_to_node(ind_x)] : index->out_labels[index_to_node(ind_x)];
//
//		for(;;){
//			LabelEntry &l_n = label[index_n];
//			LabelEntry &l_h = target_label[index_h];
//			if(l_n.v==NULL_NODE || l_h.v == NULL_NODE || l_n.v >= ind_x)
//				break;
//
//			//no label entry of h=index_to_node(ind_x) can be larger in terms of vertex ordering than ind_x
//			assert(l_h.v <= ind_x);
//
//			if(l_n.v < l_h.v || l_n.d > d)
//				++index_n;
//			else if(l_n.v > l_h.v  || l_h.d > d)
//				++index_h;
//			else{
//				assert(l_n.v == l_h.v);
//
//				if((custom_weight)l_n.d + (custom_weight)l_h.d < d){
//					label.erase(label.begin()+t);
//					data->removed++;
//					t--;
//					break;
//				}
//				else if((custom_weight)l_n.d + (custom_weight)l_h.d == d){
//					if(l_n.v < ind_x){
//						label.erase(label.begin()+t);
//						data->removed++;
//						t--;
//
//						break;
//
//					}
//				}
//				++index_n;
//				++index_h;
//
//
//			}
//		}
//
//
//	}
//
//
//
//}


bool Labeling_Tools::isIndexInSet(custom_node index, std::vector<custom_node>& hubset){

	size_t i = 0;

	assert(hubset.size()>0);


	if(hubset.size()==1)
		return hubset[0] == index;

	if(hubset[0] > index || hubset[hubset.size()-1] < index)
		return false;

	size_t mid, left = 0 ;
	size_t right = hubset.size(); // one position passed the right end

	while (left < right) {
		mid = left + round((right - left)/2);

		assert(hubset[mid]!=NULL_NODE);

		if (index > hubset[mid])
			left = mid+1;

		else if (index < hubset[mid])
			right = mid;
		else
			return true;

	}

	return false;
}


bool Labeling_Tools::test_affected_via_hubs(custom_node n1, custom_node n2, std::vector<bool>& isA, custom_weight d){

	custom_weight comparison_d;
	std::vector<custom_node> hubset = index->hubs(n1,n2,comparison_d);

	if(d!=comparison_d){
		assert(d>comparison_d);
		return false;
	}

	assert(comparison_d!=NULL_NODE);
	assert(!hubset.empty());

	if(hubset.size()>1){

		for(custom_node & el : hubset){
			if(isA[el]){
				//hubset.clear();
				return true;
			}
		}
		//hubset.clear();
		return false;
	}
	else{
		//qui invece sono "costretto" a considerare hubset[0]==node_to_index(n1) || hubset[0]==node_to_index(n2) perchè non so se esiste o meno un altro cammino di pari peso
		//ossia l'hub n1 o n2 copre un cammino, non so se questo cammino include l'arco o meno, per "sicurezza" lo devo considerare affected
		if(hubset[0]==node_to_index(n1) || hubset[0]==node_to_index(n2) || isA[hubset[0]]){
			hubset.clear();
			return true;
		}
		//hubset.clear();
		return false;
	}





}
//

bool Labeling_Tools::test_affected_via_distances(
		custom_node n1, custom_node n2, std::vector<bool>& isA,
		std::vector<custom_weight>& d_x,
		std::vector<custom_weight>& d_y,
		std::vector<custom_weight>& d_prime_x,
		std::vector<custom_weight>& d_prime_y,
		std::vector<custom_weight>& d_toward_x,
		std::vector<custom_weight>& d_toward_y,
		std::vector<custom_weight>& d_from_x,
		std::vector<custom_weight>& d_from_y,

		std::vector<custom_weight>& d_prime_from_x,
		std::vector<custom_weight>& d_prime_toward_y


){



	std::vector<custom_node> hubset = index->hubs(n1,n2);


#ifndef NDEBUG

	if(graph->isDirected()){
		assert((n2==y && d_prime_toward_y[n1] == d_toward_y[n1]) || (n1==x && d_prime_from_x[n2] == d_from_x[n2]));
		assert(!hubset.empty());
		custom_weight comparison_d;
		std::vector<custom_node> h = index->hubs(n1,n2,comparison_d);
		assert(comparison_d==index->query(n1,n2));

		if(n2==y){
			if(d_toward_y[n1] != d_toward_x[n1] + temp_weight)
				assert(false);
		}
		if(n1==x){
			if(d_from_x[n2] != d_from_y[n2] + temp_weight)
				assert(false);
		}
	}
	else{
		assert((n2==y && d_prime_y[n1] == d_y[n1]) || (n1==x && d_prime_x[n2] == d_x[n2]));

		assert(!hubset.empty());
		custom_weight comparison_d;
		std::vector<custom_node> h = index->hubs(n1,n2,comparison_d);
		assert(comparison_d==index->query(n1,n2));

		if(n2==y){
			if(d_y[n1] != d_x[n1] + temp_weight)
				assert(false);
		}
		if(n1==x){
			if(d_x[n2] != d_y[n2] + temp_weight)
				assert(false);
		}
	}

#endif


	if(hubset.size()>1){


#ifndef NDEBUG
		bool found_n1 = false;
		bool found_n2 = false;

		for(custom_node & el : hubset){
			if(el==node_to_index(n1))
				found_n1 = true;
			if(el==node_to_index(n2))
				found_n2 = true;
		}

		assert(!found_n2 && !found_n1);

#endif
		//E' GIUSTO PRENDERE TRUE SE E' VERO PER UN SOLO HUB NON DEVE VALERE PER TUTTI
		for(custom_node & el : hubset){
			if(isA[el]){
				//hubset.clear();
				return true;
			}
		}
		//hubset.clear();

		return false;
	}
	else{

		//se la distanza non cambia e l'hub è n1 oppure n2, allora sicuramente non è affected (resta valida la label entry)

#ifndef NDEBUG
		if(isA[hubset[0]])
			assert(hubset[0]!=node_to_index(n1) && hubset[0]!=node_to_index(n2));

#endif
		if(isA[hubset[0]]){
			//hubset.clear();
			return true;
		}
		//hubset.clear();
		return false;
	}


}




int Labeling_Tools::removal_of_outdated() {


	INFO("Removal of Outdated Labels");

	boost::timer::auto_cpu_timer t;

	int removed_entries = 0;

	assert(std::count(affected_vector_x.begin(),affected_vector_x.end(),true)==affected_x.size());
	assert(std::count(affected_vector_y.begin(),affected_vector_y.end(),true)==affected_y.size());



	for(size_t t = 0;t<affected_x.size();t++){

		int local_remove = 0;
		custom_node fn_index = affected_x[t];
		custom_node first_node = index_to_node(fn_index);
		int bias = affected_vector_y[fn_index] ? 1 : 0;

		std::vector<LabelEntry> &l1 = graph->isDirected() ? index->out_labels[first_node] : index->in_labels[first_node];
		assert(l1[l1.size()-1].v==NULL_NODE);
		assert(l1[l1.size()-1].d==NULL_WEIGHT);

		int prv = -1;

		for(int t = l1.size()-2;t>=0;t--){
			custom_node idx = l1[t].v;
			assert(idx!=NULL_NODE);

			if (affected_vector_y[idx] && idx!=fn_index) {
				if(prv!=-1 && prv == (int) idx){
					prv = (int) idx;
					l1.erase(l1.begin()+t);
				}
				else{
					prv = (int) idx;
					l1.erase(l1.begin()+t);
					local_remove++;
				}
			}
			else{
				prv = (int) idx;
			}
			if(local_remove==affected_y.size()-bias)
				break;

		}
		removed_entries+=local_remove;


	}


	for(size_t t = 0;t<affected_y.size();t++){

		custom_node fn_index = affected_y[t];
		custom_node first_node = index_to_node(fn_index);
		int local_remove = 0;
		int bias = affected_vector_x[fn_index] ? 1 : 0;
		std::vector<LabelEntry> &l1 = index->in_labels[first_node];
		assert(l1[l1.size()-1].v==NULL_NODE);
		assert(l1[l1.size()-1].d==NULL_WEIGHT);
		int prv = -1;


		for(int t = l1.size()-2;t>=0;t--){
			custom_node idx = l1[t].v;
			assert(idx!=NULL_NODE);
			if (affected_vector_x[idx] && idx!=fn_index) {
				if(prv!=-1 && prv== (int) idx){
					prv = (int) idx;
					l1.erase(l1.begin()+t);
				}
				else{
					prv = (int) idx;
					l1.erase(l1.begin()+t);
					local_remove++;
				}
			}
			else {
				prv = (int) idx;
			}

			//se fn_index affected sicuramente non lo rimuovo, posso killare prima la search
			if(local_remove==affected_x.size()-bias){
				//INFO("BREAK");
				break;
			}

		}

		removed_entries+=local_remove;

	}

	return removed_entries;
}








int Labeling_Tools::restore_cover(custom_node index_root, std::vector<bool>& aff_vect, bool forward_visit, int num) {


	int amount_of_added = 0;

	custom_node root = index_to_node(index_root);

	std::vector<LabelEntry>& idx_r = forward_visit && graph->isDirected() ? index->out_labels[root] : index->in_labels[root];

	if(!graph->isWeighted()){
		assert(unweighted_empty());
		unweighted_insert(root,0);
		assert(!unweighted_empty());

	}
	else{
		assert(weighted_empty());
		weighted_insert(root,0);
		assert(!weighted_empty());

	}

	custom_node v;
	custom_weight v_dist;
	custom_node v_index;

	int seen = 0;


	while (graph->isWeighted() ? !weighted_empty() : !unweighted_empty())   {
		if(seen==num)
			break;


		if(!graph->isWeighted()){
			v = (*unweighted_get_min()).first;
			v_dist = (*unweighted_get_min()).second;
			v_index = node_to_index(v);
			unweighted_pop();
		}
		else{
			v = weighted_get_min_node();
			v_dist = weighted_get_min_dist();
			v_index = node_to_index(v);
			weighted_pop();
			status_prio_que[v]=2;

		}


		if(aff_vect[v_index])
			seen++;

		if(v_index < index_root)
			continue;


		std::vector<LabelEntry>& idx_v =  (!graph->isDirected() || forward_visit) ? index->in_labels[v] : index->out_labels[v];

		if(!aff_vect[v_index])
			if(hasIndex(idx_v,index_root)==NULL_INDEX)
				goto r_c_prune;

		if(aff_vect[v_index]){

			custom_weight provided_dist = NULL_WEIGHT;
			size_t t_r = 0;
			size_t t_v = 0;
			for (;;) {

				LabelEntry &l_r = idx_r[t_r];
				LabelEntry &l_v = idx_v[t_v];
				if(l_r.v==NULL_NODE || l_v.v==NULL_NODE)
					break;
				if(l_r.v>=index_root || l_v.v>=index_root)
					break;
				if(l_r.v>=v_index || l_v.v>=v_index)
					break;
				assert(l_r.d < NULL_WEIGHT && l_v.d < NULL_WEIGHT);
				if (l_r.v < l_v.v){
					t_r++;
				}
				else if(l_r.v > l_v.v){
					t_v++;
				}

				else{

					assert(l_r.v == l_v.v);
					assert(l_r.v<=index_root && l_r.v<=v_index);

					if(v_dist>l_r.d+l_v.d){
						goto r_c_prune;
					}
					else if(v_dist==l_r.d+l_v.d){

						if(l_r.v>index_root){
							std::cout<<l_r.v<<" "<<index_root<<"\n";
						}
						assert(l_r.v<=index_root);
						if(l_r.v<index_root)
							goto r_c_prune;

					}

					else{
						provided_dist=std::min(provided_dist, l_r.d+l_v.d);
						t_v++;
						t_r++;

					}


				}
			}

			assert(provided_dist>v_dist);
			size_t suited_pos = 0;
			for (;idx_v[suited_pos].v<=index_root; ++suited_pos) {
				LabelEntry &l = idx_v[suited_pos];
				if (l.v == index_root)
					break;
			}
			for (int j = (int)idx_v.size() - 1; j - 1 >= (int)suited_pos; --j)
				idx_v[j] = idx_v[j - 1];

			idx_v.push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));
			idx_v[suited_pos] = LabelEntry(index_root,v_dist);
			amount_of_added++;

		}

		if(!graph->isWeighted())
			unweighted_relax(v,v_dist,forward_visit);

		else
			weighted_relax(v,v_dist,forward_visit);




		r_c_prune:{


		};


	}

	if(!graph->isWeighted())
		unweighted_reset();
	else
		weighted_reset();


	assert(graph->isWeighted() ? weighted_empty() : unweighted_empty());
	return amount_of_added;

}

void Labeling_Tools::former_restore(UpdateData* data) {


	boost::timer::auto_cpu_timer t;

	bool directed = graph->isDirected();


	if(affected_x.empty() || affected_y.empty()){

#ifndef NDEBUG
		assert(graph->isWeighted());
#endif
		INFO("Restore unnecessary, removed/weight-changed non-shortest-path edge");
		return;

	}


	ProgressStream updater_(affected_x.size()+affected_y.size());
	updater_.label()<< "Restoring label entries, performing " <<affected_x.size()+affected_y.size()<< " visits ("<<affected_x.size()<<","<<affected_y.size()<<") ";


	custom_node min_t = std::min(min_affected_x,min_affected_y);
	custom_node max_t = std::max(max_affected_x,max_affected_y);

	for(custom_node root_index = min_t;root_index<=max_t;root_index++){


		if(!affected_vector_x[root_index] && !affected_vector_y[root_index])
			continue;
		if(affected_vector_x[root_index]){
			std::vector<bool>& from_vector = affected_vector_x;
			std::vector<bool>& to_vector = affected_vector_y;

			data->added += restore_cover(root_index,to_vector,true,affected_y.size());
			++updater_;
			affected_vector_x[root_index]=false;
		}
		if(affected_vector_y[root_index]){


			std::vector<bool>& from_vector = affected_vector_y;
			std::vector<bool>& to_vector = affected_vector_x;

			data->added += restore_cover(root_index,to_vector,false,affected_x.size());
			++updater_;
			affected_vector_y[root_index]=false;
		}



	}

//	for(custom_node root_index = min_t;root_index<=max_t;root_index++){
//		if(affected_vector_x[root_index])
//			affected_vector_x[root_index]=false;
//		if(affected_vector_y[root_index])
//			affected_vector_y[root_index]=false;
//	}

	affected_x.clear();
	affected_y.clear();



}

void Labeling_Tools::restore(UpdateData* data) {


	bool directed = graph->isDirected();


	if(affected_x.empty() || affected_y.empty()){

#ifndef NDEBUG
		assert(graph->isWeighted());
#endif
		INFO("Restore unnecessary, removed/weight-changed non-shortest-path edge");
		return;

	}






	evaluation(data);
	restoration(data);


}

void Labeling_Tools::restoration(UpdateData*data){

	INFO("Restore of New Labels");
	boost::timer::auto_cpu_timer t;



	ProgressStream updater_(affected_x.size()+affected_y.size());

	updater_.label()<< "Restoring label entries, performing (" <<affected_x.size()<<"+"<<affected_y.size()<< ") phases "  ;

	//double num_visits = 0;

	custom_node min_t = std::min(min_affected_x,min_affected_y);
	custom_node max_t = std::max(max_affected_x,max_affected_y);

	for(custom_node root_index = min_t;root_index<=max_t;root_index++){


		if(!affected_vector_x[root_index] && !affected_vector_y[root_index])
			continue;

		if(affected_vector_y[root_index]){

			for(size_t k = 0;k<distances_of_affected_x_toward_affected_y.size();k++){
				if(affected_x[k]>max_affected_x)
					break;
				if(affected_x[k]<root_index || distances_of_affected_x_toward_affected_y[k].empty())
					continue;
				std::unordered_map<custom_node,custom_weight>::iterator p = distances_of_affected_x_toward_affected_y[k].find(root_index);
				if(p==distances_of_affected_x_toward_affected_y[k].end())
					continue;

				assert(affected_x[k]>=root_index);
				assert(p->first==root_index);
				assert(p->second<NULL_WEIGHT);
				aux_insert(index_to_node(affected_x[k]),p->second);
			}

			while(!aux_empty()){
				custom_node cu = aux_min_node();
				custom_weight wg = aux_min_dist();
				assert((*H[cu]).node==cu);
				assert((*H[cu]).weight==-wg);
				aux_pop();
				assert(node_to_index(cu)>=root_index);
				resume_shortest_paths_over_affected(root_index, cu, wg, data,!graph->isDirected(),affected_vector_x);
				S[cu]=2;
			}

			assert(aux_empty());
			aux_reset();
			++updater_;


			affected_vector_y[root_index]=false;

		}

		if(affected_vector_x[root_index]){



			for(size_t k = 0;k<distances_of_affected_y_from_affected_x.size();k++){
				if(affected_y[k]>max_affected_y)
					break;
				if(affected_y[k]<root_index || distances_of_affected_y_from_affected_x[k].empty())
					continue;
				std::unordered_map<custom_node,custom_weight>::iterator p = distances_of_affected_y_from_affected_x[k].find(root_index);
				if(p==distances_of_affected_y_from_affected_x[k].end())
					continue;
				assert(affected_y[k]>=root_index);
				assert(p->first==root_index);
				assert(p->second<NULL_WEIGHT);
				aux_insert(index_to_node(affected_y[k]),p->second);

			}
			while(!aux_empty()){
				custom_node cu = aux_min_node();
				custom_weight wg = aux_min_dist();
				assert((*H[cu]).node==cu);
				assert((*H[cu]).weight==-wg);
				aux_pop();
				assert(node_to_index(cu)>=root_index);

				resume_shortest_paths_over_affected(root_index, cu, wg, data, true,affected_vector_y);
				S[cu]=2;

			}
			assert(aux_empty());
			aux_reset();

			++updater_;
			affected_vector_x[root_index]=false;
		}



	}
	affected_x.clear();
	affected_y.clear();

}
void Labeling_Tools::evaluation(UpdateData*data){

	INFO("Evaluation of Status of Labeling");
	boost::timer::auto_cpu_timer t;

	distances_of_affected_x_toward_affected_y.clear();
	distances_of_affected_y_from_affected_x.clear();

	distances_of_affected_x_toward_affected_y.resize(affected_x.size());
	distances_of_affected_y_from_affected_x.resize(affected_y.size());

	ProgressStream progress(affected_x.size()+affected_y.size());
	progress.label()<< "Evaluating remaining label entries onto " <<affected_x.size()+affected_y.size()<< " nodes ("<<affected_x.size()<<","<<affected_y.size()<<") ";

	size_t tracker = 0;

	while(true){
		if(tracker == affected_x.size())
			break;
		evaluate(tracker++,affected_x[tracker], data, true);
		++progress;
	}

	tracker = 0;

	while(true){
		if(tracker == affected_y.size())
			break;
		evaluate(tracker++,affected_y[tracker], data,false);
		++progress;

	}

#ifndef NDEBUG
	assert(distances_of_affected_x_toward_affected_y.size()==affected_x.size());
	for(auto&el : distances_of_affected_x_toward_affected_y)
		assert(el.size()<=affected_y.size());
	assert(distances_of_affected_y_from_affected_x.size()==affected_y.size());
	for(auto&el : distances_of_affected_y_from_affected_x)
		assert(el.size()<=affected_x.size());
#endif

}
void Labeling_Tools::evaluate(custom_node id, custom_node v_index,UpdateData* data, bool afct_x) {



	custom_node vertex = index_to_node(v_index);
	std::unordered_map<custom_node,custom_weight>::iterator p;

	if(afct_x){


		graph->forNeighborsOf(vertex, [&](custom_node n) {
			custom_node idx = node_to_index(n);

#ifndef NDEBUG
			if(affected_vector_x[idx]){
				const std::vector<LabelEntry>& label_ = graph->isDirected() ? index->out_labels[n] : index->in_labels[n];
				for (size_t i1 = 0;;i1++) {

					if (label_[i1].v >= v_index || label_[i1].v == NULL_NODE)
						break;

					if(affected_vector_y[label_[i1].v])
						assert(label_[i1].v==node_to_index(n));
				}
			}
#endif

			if(!affected_vector_x[idx]){
				const std::vector<LabelEntry>& label_ = graph->isDirected() ? index->out_labels[n] : index->in_labels[n];

				for (size_t i1 = 0;;i1++) {
					//when sentinel is reached, continue
					if (label_[i1].v >= v_index || label_[i1].v == NULL_NODE)
						break;

					if(affected_vector_y[label_[i1].v]){
						p = distances_of_affected_x_toward_affected_y[id].find(label_[i1].v);
						if(p!=distances_of_affected_x_toward_affected_y[id].end())
							p->second = std::min(p->second,label_[i1].d + graph->weight(vertex,n));
						else
							distances_of_affected_x_toward_affected_y[id][label_[i1].v] = label_[i1].d + graph->weight(vertex,n);
					}

				}
			}

			else{
				if(v_index>idx){
					if(affected_vector_y[idx]){
						p = distances_of_affected_x_toward_affected_y[id].find(idx);
						if(p!=distances_of_affected_x_toward_affected_y[id].end())
							p->second = std::min(p->second,0 + graph->weight(vertex,n));
						else
							distances_of_affected_x_toward_affected_y[id][idx] = 0 + graph->weight(vertex,n);
					}
				}
			}

		});
	}
	else{
		if(!graph->isDirected()){

			graph->forNeighborsOf(vertex, [&](custom_node n) {
				custom_node idx = node_to_index(n);

#ifndef NDEBUG

				if(affected_vector_y[idx]){
					const std::vector<LabelEntry> &label_ = index->in_labels[n];
					for (size_t i1 = 0;;i1++) {
						//when sentinel is reached, continue
						if (label_[i1].v >= v_index || label_[i1].v == NULL_NODE)
							break;

						if(affected_vector_x[label_[i1].v])
							assert(label_[i1].v==idx);
					}
				}

#endif
				if(!affected_vector_y[idx]){
					const std::vector<LabelEntry> &label_ = index->in_labels[n];
					for (size_t i1 = 0;;i1++) {
						//when sentinel is reached, continue
						if (label_[i1].v >= v_index || label_[i1].v == NULL_NODE)
							break;

						if(affected_vector_x[label_[i1].v]){
							p = distances_of_affected_y_from_affected_x[id].find(label_[i1].v);
							if(p!=distances_of_affected_y_from_affected_x[id].end())
								p->second = std::min(p->second,label_[i1].d + graph->weight(vertex,n));
							else
								distances_of_affected_y_from_affected_x[id][label_[i1].v] = label_[i1].d + graph->weight(vertex,n);
						}

					}
				}
				else{
					if(v_index>idx){
						if(affected_vector_x[idx]){
							p = distances_of_affected_y_from_affected_x[id].find(idx);
							if(p!=distances_of_affected_y_from_affected_x[id].end())
								p->second = std::min(p->second,0 + graph->weight(vertex,n));
							else
								distances_of_affected_y_from_affected_x[id][idx] = 0 + graph->weight(vertex,n);
						}
					}
				}
			});

		}
		else{

			graph->forInNeighborsOf(vertex, [&](custom_node n) {
				custom_node idx = node_to_index(n);

#ifndef NDEBUG

				if(affected_vector_y[idx]){
					const std::vector<LabelEntry> &label_ = index->in_labels[n];
					for (size_t i1 = 0;;i1++) {
						if (label_[i1].v >= v_index || label_[i1].v == NULL_NODE)
							break;

						if(affected_vector_x[label_[i1].v])
							assert(label_[i1].v==node_to_index(n));
					}
				}

#endif
				if(!affected_vector_y[idx]){

					const std::vector<LabelEntry> & label_ = index->in_labels[n];


					for (size_t i1 = 0;;i1++) {
						if (label_[i1].v >= v_index || label_[i1].v == NULL_NODE)
							break;
						if(affected_vector_x[label_[i1].v]){
							p = distances_of_affected_y_from_affected_x[id].find(label_[i1].v);
							if(p!=distances_of_affected_y_from_affected_x[id].end())
								p->second = std::min(p->second,label_[i1].d + graph->weight(n,vertex));
							else
								distances_of_affected_y_from_affected_x[id][label_[i1].v] = label_[i1].d + graph->weight(n,vertex);
						}

					}
				}
				else{

					if(v_index>idx){
						if(affected_vector_x[idx]){
							p = distances_of_affected_y_from_affected_x[id].find(idx);
							if(p!=distances_of_affected_y_from_affected_x[id].end())
								p->second = std::min(p->second,0 + graph->weight(n,vertex));
							else
								distances_of_affected_y_from_affected_x[id][idx] = 0 + graph->weight(n,vertex);

						}
					}
				}

			});



		}



	}

}





//bool Labeling_Tools::isAffectedX(custom_node i){
//	return affected_vector_x[i];
//}
//bool Labeling_Tools::isAffectedY(custom_node i){
//	return affected_vector_y[i];
//}
void Labeling_Tools::setAffectedX(custom_node i){
	assert(!affected_vector_x[i]);
	affected_x.push_back(i);
	affected_vector_x[i] = true;
	min_affected_x=std::min(min_affected_x,i);
	max_affected_x=std::max(max_affected_x,i);

}
void Labeling_Tools::setAffectedY(custom_node i){
	assert(!affected_vector_y[i]);
	affected_y.push_back(i);
	affected_vector_y[i] = true;
	min_affected_y=std::min(min_affected_y,i);
	max_affected_y=std::max(max_affected_y,i);
}



void Labeling_Tools::resume_shortest_paths_over_affected(custom_node bfs_i,custom_node vertex, custom_weight distance, UpdateData* data, bool forward, std::vector<bool>& A) {

#ifndef NDEBUG
	assert((graph->isWeighted() && status_prio_que.size()>= graph->upperNodeIdBound()) || (!graph->isWeighted() && (custom_node)queue.size() >= graph->upperNodeIdBound()));
	if(S[vertex]!=1){
		std::cout<<vertex<<" "<<S[vertex]<<"\n";
		assert(false);
	}
#endif

	custom_node root = index_to_node(bfs_i);

	std::vector<LabelEntry> &idx_r = graph->isDirected() && forward ? index->out_labels[root] : index->in_labels[root];



	for (int i = (int)idx_r.size() - 1; i >= 0; --i) {
		//when sentinel is reached, continue
		if (idx_r[i].v == NULL_NODE)
			continue;
		root_label[idx_r[i].v] = idx_r[i].d;
	}

	if(!graph->isWeighted()){

		assert((custom_node)queue.size() >= graph->upperNodeIdBound());
		assert(!marked[vertex]);
		assert(que_h==0 && que_t==0);
		unweighted_insert(vertex,distance);


		while (!unweighted_empty()) {

			custom_node v = (*unweighted_get_min()).first;
			custom_weight d = (*unweighted_get_min()).second;
			unweighted_pop();


			{
				//PATH FROM v to bfs_i equal or shorter than that already in the global QUEUE
				if(S[v]==2 || (S[v]==1 && aux_key(v)<d)){
					assert(v!=vertex);
					goto partial_prune_;
				}
				if(S[v]==1 && aux_key(v)>d){
					assert(v!=vertex);
					aux_erase(v);
				}
				if(bfs_i > node_to_index(v))
					goto partial_prune_;


				std::vector<LabelEntry> &idx_v = forward ? index->in_labels[v] : index->out_labels[v];


				size_t i = 0;
				for (; idx_v[i].v <= bfs_i; ++i) {
					LabelEntry &l = idx_v[i];
					if (root_label[l.v] != NULL_WEIGHT && ((root_label[l.v] + l.d) <= d)){ //if query through l.v returns smaller value, prune on
						goto partial_prune_;
					}
					if (l.v == bfs_i)
						break;
				}

				assert(bfs_i<=node_to_index(v));

				if(bfs_i==idx_v[i].v){
					idx_v[i] = LabelEntry(bfs_i,d);
					//replacement
				}
				else{
					for (int j = (int)idx_v.size() - 1; j - 1 >= (int)i; --j)
						idx_v[j] = idx_v[j - 1];
					idx_v.push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));
					idx_v[i] = LabelEntry(bfs_i,d);
					data->added++;

				}
			}
			unweighted_affected_relax(v,d,forward,A);
			partial_prune_:{};
		}
		unweighted_reset();



	}
	else{
		assert(weighted_empty());
		weighted_insert(vertex,distance);

		while (!weighted_empty()) {

			custom_node v = weighted_get_min_node();
			custom_weight d = weighted_get_min_dist();
			status_prio_que[v]=2;

			weighted_pop();


			{

				//PATH FROM v to bfs_i equal or shorter than that already in the global QUEUE
				if(S[v]==2 || (S[v]==1 && aux_key(v)<d)){
					assert(v!=vertex);
					goto w_partial_prune_;
				}
				if(S[v]==1 && aux_key(v)>d){
					assert(v!=vertex);
					aux_erase(v);
				}
				if(bfs_i > node_to_index(v))
					goto w_partial_prune_;

				std::vector<LabelEntry> &idx_v = forward ? index->in_labels[v] : index->out_labels[v];

				size_t i = 0;
				for (; idx_v[i].v <= bfs_i; ++i) {
					LabelEntry &l = idx_v[i];

					if (root_label[l.v] != NULL_WEIGHT && ((root_label[l.v] + l.d) <= d))
						//if query through l.v returns smaller/equal value, prune on
						goto w_partial_prune_;

					if (l.v == bfs_i)
						break;
				}

				assert(bfs_i<=node_to_index(v));

				if(bfs_i==idx_v[i].v)
					idx_v[i] = LabelEntry(bfs_i,d);
				else{
					for (int j = (int)idx_v.size() - 1; j - 1 >= (int)i; --j)
						idx_v[j] = idx_v[j - 1];
					idx_v.push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));
					idx_v[i] = LabelEntry(bfs_i,d);
					data->added++;

				}



			}
			weighted_affected_relax(v,d,forward,A);
			w_partial_prune_:{};
		}
		weighted_reset();

	}




	for (int i = (int)idx_r.size() - 1; i >= 0; --i) {
		if (idx_r[i].v == NULL_NODE)
			continue;
		root_label[idx_r[i].v] = NULL_WEIGHT;
	}

}

void Labeling_Tools::check_closure(UpdateData* data,
		custom_node root_index, custom_node v, custom_weight dist,bool forward){


	assert(aux_empty());
	custom_node v_index = node_to_index(v);
	assert(root_index<v_index);
	aux_insert(v,0);



	while(!aux_empty()){
		custom_node cu = aux_min_node();
		custom_weight wg = aux_min_dist();
		assert((*H[cu]).node==cu);
		assert((*H[cu]).weight==-wg);
		aux_pop();
		assert(node_to_index(cu)>=v_index);

		S[cu]=2;

		std::vector<LabelEntry> &idx_cu = graph->isDirected() && forward ? index->out_labels[cu] : index->in_labels[cu];

		custom_node seek_index_v = hasIndex(idx_cu,v_index);
		if(seek_index_v == NULL_INDEX)
			continue;

		custom_node seek_index_r = hasIndex(idx_cu,root_index);


		if(seek_index_r != NULL_INDEX){
			if(idx_cu[seek_index_r].d+dist<=idx_cu[seek_index_v].d){
				idx_cu.erase(idx_cu.begin()+seek_index_v);
				data->removed++;
			}
			else
				continue;

		}
		if(!graph->isDirected() || !forward){
			graph->forNeighborsOf(cu, [&](custom_node tv) {
				if(S[tv]==0 && node_to_index(tv)>v_index)
					aux_insert(tv,wg+graph->weight(cu, tv));
				else if(S[tv]==1 && aux_key(tv)>wg+graph->weight(cu, tv)){
					assert(node_to_index(tv)>v_index);
					aux_dec_key(tv,wg+graph->weight(cu, tv));
				}

			});
		}

		else{
			graph->forInNeighborsOf(cu, [&](custom_node tv) {
				if(S[tv]==0 && node_to_index(tv)>v_index)
					aux_insert(tv,wg+graph->weight(tv, cu));
				else if(S[tv]==1 && aux_key(tv)>wg+graph->weight(tv, cu)){
					assert(node_to_index(tv)>v_index);
					aux_dec_key(tv,wg+graph->weight(tv, cu));
				}
			});

		}



	}


	aux_reset();










}

void Labeling_Tools::resume_shortest_paths(custom_node bfs_i,custom_node vertex, custom_weight distance, UpdateData* data, bool forward, bool minimal) {

	assert((graph->isWeighted() && status_prio_que.size()>= graph->upperNodeIdBound()) || (!graph->isWeighted() && (custom_node)queue.size() >= graph->upperNodeIdBound()));






	custom_node root = index_to_node(bfs_i);

	std::vector<LabelEntry> &idx_r = graph->isDirected() && forward ? index->out_labels[root] : index->in_labels[root];



	for (int i = (int)idx_r.size() - 1; i >= 0; --i) {
		if (idx_r[i].v == NULL_NODE)
			continue;
		root_label[idx_r[i].v] = idx_r[i].d;

	}




	if(!graph->isWeighted()){

		assert((custom_node)queue.size() >= graph->upperNodeIdBound());
		assert(!marked[vertex]);
		assert(que_h==0 && que_t==0);
		unweighted_insert(vertex,distance);

		while (!unweighted_empty()) {

			custom_node v = (*unweighted_get_min()).first;
			custom_weight d = (*unweighted_get_min()).second;

			unweighted_pop();

			{
				if(bfs_i > node_to_index(v))
					continue;

				assert(forward || graph->isDirected());
				std::vector<LabelEntry> &idx_v = !graph->isDirected() || forward ? index->in_labels[v] : index->out_labels[v];


				size_t i = 0;
				for (; idx_v[i].v <= bfs_i; ++i) {
					LabelEntry &l = idx_v[i];
					if (root_label[l.v] != NULL_WEIGHT && ((root_label[l.v] + l.d) <= d)){ //if query through l.v returns smaller value, prune on
						goto partial_prune;
					}


					if (l.v == bfs_i)
						break;
				}
				assert(bfs_i<=node_to_index(v));

				//distance from root to v has been shortened here

				if(bfs_i==idx_v[i].v)
					idx_v[i] = LabelEntry(bfs_i,d);

				else{
					for (int j = (int)idx_v.size() - 1; j - 1 >= (int)i; --j)
						idx_v[j] = idx_v[j - 1];
					idx_v.push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));
					idx_v[i] = LabelEntry(bfs_i,d);
					data->added++;
				}
				assert(idx_v[i].v==bfs_i);
				assert(idx_v[i].d==d);

				if(minimal){
					clean_target(data,bfs_i,i,idx_v,d,forward);
					check_closure(data,bfs_i,v,d,forward);

				}

			}

			unweighted_relax(v,d,forward);
			partial_prune:{};
		}

		unweighted_reset();
	}
	else{
		assert(weighted_empty());
		weighted_insert(vertex,distance);

		while (!weighted_empty()) {

			custom_node v = weighted_get_min_node();
			custom_weight d = weighted_get_min_dist();
			status_prio_que[v]=2;

			weighted_pop();
			{
				if(bfs_i > node_to_index(v))
					continue;
				assert(forward || graph->isDirected());

				std::vector<LabelEntry> &idx_v = !graph->isDirected() || forward ? index->in_labels[v] : index->out_labels[v];


				size_t i = 0;
				for (; idx_v[i].v <= bfs_i; ++i) {
					LabelEntry &l = idx_v[i];

					if (root_label[l.v] != NULL_WEIGHT && ((root_label[l.v] + l.d) <= d)){ //if query through l.v returns smaller/equal value, prune on
						goto w_partial_prune;
					}


					if (l.v == bfs_i)
						break;
				}
				assert(bfs_i<=node_to_index(v));

				//either adds or overwrite

				if(bfs_i==idx_v[i].v)
					idx_v[i] = LabelEntry(bfs_i,d);

				else{
					for (int j = (int)idx_v.size() - 1; j - 1 >= (int)i; --j)
						idx_v[j] = idx_v[j - 1];

					idx_v.push_back(LabelEntry(NULL_NODE,NULL_WEIGHT));
					idx_v[i] = LabelEntry(bfs_i,d);
					data->added++;
				}
				assert(idx_v[i].v==bfs_i);
				assert(idx_v[i].d==d);


				if(minimal){
					clean_target(data,bfs_i,i,idx_v,d,forward);
					check_closure(data,bfs_i,v,d,forward);
				}





			}

			weighted_relax(v,d,forward);
			w_partial_prune:{};
		}

		weighted_reset();


	}
	for (int i = (int)idx_r.size() - 1; i >= 0; --i) {
		if (idx_r[i].v == NULL_NODE)
			continue;
		root_label[idx_r[i].v] = NULL_WEIGHT;
	}




}

void Labeling_Tools::clean_target(UpdateData* data, custom_node bfs_i, custom_node indice, std::vector<LabelEntry> &idx_v,custom_weight d, bool forward){

	size_t t = indice;
	++t;
	for (;;) {

		LabelEntry &l = idx_v[t];
		if(l.v==NULL_NODE)
			break;
		custom_node hub_node = index_to_node(l.v);

		std::vector<LabelEntry> &hub_label = graph->isDirected() && forward ? index->out_labels[hub_node] : index->in_labels[hub_node];



		assert(hub_label[hasIndex(hub_label,l.v)].d==0);

		custom_node seek_index = hasIndex(hub_label,bfs_i);
		if(seek_index!=NULL_INDEX){
			if(hub_label[seek_index].d+d<=l.d){
				idx_v.erase(idx_v.begin()+t);
				data->removed++;
				continue;
			}
		}

		++t;


	}
}
void Labeling_Tools::clean(){



	for(custom_node &el : affected_x){
		assert(affected_vector_x[el]);
		affected_vector_x[el] = false;
	}
	for(custom_node &el : affected_y){
		assert(affected_vector_y[el]);
		affected_vector_y[el] = false;
	}
	affected_x.clear();
	affected_y.clear();

#ifndef NDEBUG

	assert(std::count(affected_vector_x.begin(),affected_vector_x.end(),true)==0);
	assert(std::count(affected_vector_y.begin(),affected_vector_y.end(),true)==0);



	assert(affected_x.empty());
	assert(affected_y.empty());

#endif




};





Labeling_Tools::~Labeling_Tools() {
	if(graph->isWeighted()){
		delete prio_que;
		delete[] handles_prio_que;
	}
	if(dynamic){
		root_label.clear();
		affected_vector_x.clear();
		affected_vector_y.clear();
		delete Q;
		delete[] H;
		S.clear();
	}
}

//
//

custom_node Labeling_Tools::hasIndex(std::vector<LabelEntry> & label_vect, custom_node index_v){


	size_t i = 0;
	assert(label_vect[label_vect.size()-1].v==NULL_NODE);

	if(label_vect.size()==1)
		return NULL_INDEX;
	if(label_vect.size()==2){
		if(label_vect[0].v == index_v)
			return 0;
		else
			return NULL_INDEX;
	}
	if(label_vect[0].v > index_v || label_vect[label_vect.size()-2].v < index_v)
		return NULL_INDEX;

	size_t mid, left = 0 ;
	size_t right = label_vect.size()-1; // one position passed the right end

	while (left < right) {
		mid = left + round((right - left)/2);

		assert(label_vect[mid].v!=NULL_NODE);

		if (index_v > label_vect[mid].v)
			left = mid+1;

		else if (index_v < label_vect[mid].v)
			right = mid;

		else
			return (custom_node)mid;

	}

	return NULL_INDEX;


}

custom_node Labeling_Tools::index_to_node(custom_node indice){
	return keeper.first[indice];
}

custom_node Labeling_Tools::node_to_index(custom_node node){
	return keeper.second[node];
}

std::vector<UpdateData> Labeling_Tools::handle_affected_comparisons(NetworKit::GraphEvent* e){
	std::vector<UpdateData> trace;
	x = e->u;
	y = e->v;
	weight = e->w;
	temp_weight = graph->weight(x,y);
	typ = e->type;






	if(!graph->isWeighted()){
		unweighted_trivial_affected(x);
		unweighted_trivial_affected(y);
	}
	else{
		weighted_trivial_affected(x);
		weighted_trivial_affected(y);
	}


	trace.push_back(UpdateData(0.0,affected_x.size(),affected_y.size(),affected_x.size()+affected_y.size(),0,0));
	clean();
	assert(affected_x.empty() && affected_y.empty());


	if(!graph->isWeighted()){
		unweighted_affected_via_hubs(x);
		unweighted_affected_via_hubs(y);
	}
	else{
		weighted_affected_via_hubs(x);
		weighted_affected_via_hubs(y);
	}

	trace.push_back(UpdateData(0.0,affected_x.size(),affected_y.size(),affected_x.size()+affected_y.size(),0,0));

	clean();
	assert(affected_x.empty() && affected_y.empty());





	affected_via_distances();
	trace.push_back(UpdateData(0.0,affected_x.size(),affected_y.size(),affected_x.size()+affected_y.size(),0,0));

	clean();
	assert(affected_x.empty() && affected_y.empty());

	return trace;


}

void Labeling_Tools::add_node_to_keeper(custom_node node, int counter){

    this->index->order.push_back(node);

    this->index->reverse_order[node]= counter;
	
}



