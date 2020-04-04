/*
 * Labeling.cpp
 *
 *  Created on: 01 feb 2016
 *      Author: Mattia D'Emidio
 */

#include "Labeling.h"

Labeling::Labeling(bool d){
	dir = d;
};

Labeling::~Labeling(){};

custom_node Labeling::index_to_node(custom_node indice){
	return order[indice];
}

custom_node Labeling::node_to_index(custom_node node){
	return reverse_order[node];
}



void Labeling::printInLabels(){

     for(int i =0 ;i < this->in_labels.size() ;i++){

		for(int j=0; j<in_labels[i].size()-1;j++){

			std::cout << "labels of : "<< i <<"nodo :"<< this->in_labels[i][j].v<< "distanza: "<<  this->in_labels[i][j].d << "\n";	
		}
	}

}

void Labeling::query(custom_node v, custom_node w, custom_weight explored_dist,std::pair<custom_node,custom_weight>& encoded) {


	assert(!dir || this->in_labels.size() == this->out_labels.size());
	assert(v<=this->in_labels.size() && w<=this->in_labels.size());


	const std::vector<LabelEntry> &s1 = dir ? this->out_labels[v] : this->in_labels[v];
	const std::vector<LabelEntry> &s2 = this->in_labels[w];

	custom_node max_index = std::max(node_to_index(v),node_to_index(w));

	for (size_t i1 = 0, i2 = 0;;) {
		if(s1[i1].v > max_index || s2[i2].v > max_index)
			break;

		else if (s1[i1].v < s2[i2].v)
			++i1;
		else if (s1[i1].v > s2[i2].v)
			++i2;
		else if (s1[i1].d > explored_dist && s1[i1].d < NULL_WEIGHT)
			++i1;
		else if (s2[i2].d > explored_dist && s2[i2].d < NULL_WEIGHT)
			++i2;
		else{
			custom_node vertex = s1[i1].v;
			assert(s1[i1].v == s2[i2].v);

			if (vertex == NULL_NODE)
				break;


			custom_weight sum = (custom_weight)s1[i1].d + (custom_weight)s2[i2].d;

			if (sum <= explored_dist){
#ifndef NDEBUG
				if(sum<query(v,w))
					assert(false);
#endif
				encoded.first = s1[i1].v;
				encoded.second = sum;
				return;
			}

		  for (++i1; s1[i1].v == vertex; ++i1){}
		  for (++i2; s2[i2].v == vertex; ++i2){}
		}
	}
	encoded.first = NULL_NODE;
	encoded.second = NULL_WEIGHT;

}

custom_weight Labeling::query(custom_node v, custom_node w) {


	  assert(!dir || this->in_labels.size() == this->out_labels.size());
	  assert(v<=this->in_labels.size() && w<=this->in_labels.size());

	  if(v==w){
		  return 0;
	  }
	  custom_weight distance = NULL_WEIGHT;


	  const std::vector<LabelEntry> &s1 = dir ? this->out_labels[v] : this->in_labels[v];
	  const std::vector<LabelEntry> &s2 = this->in_labels[w];


	  size_t i1 = 0, i2 = 0;
	  custom_node max_index = std::max(node_to_index(v),node_to_index(w));

	  while (i1 < s1.size() && i2 < s2.size()) {
		  if(s1[i1].v > max_index || s2[i2].v > max_index)
		  		break;

		  else if (s1[i1].v < s2[i2].v)
			  ++i1;
		  else if (s1[i1].v > s2[i2].v)
			  ++i2;
		  else{
			  assert(s1[i1].v == s2[i2].v);

			  custom_node vertex = s1[i1].v;

			  if (vertex == NULL_NODE)
				  break;


			  distance = std::min(distance, s1[i1].d + s2[i2].d);
			  for (++i1; s1[i1].v == vertex; ++i1){}
			  for (++i2; s2[i2].v == vertex; ++i2){}
		  }
	  }
	  //might be NULL_WEIGHT
	  return distance;

}


std::vector<custom_node> Labeling::hubs(custom_node v, custom_node w,custom_weight& distance) {


	  std::vector<custom_node> hubset;
	  assert(!dir || this->in_labels.size() == this->out_labels.size());
	  assert(v<=this->in_labels.size() && w<=this->in_labels.size());

	  if(v==w){
		  distance = 0;
		  hubset.push_back(node_to_index(v));
		  return hubset;
	  }
	  distance = NULL_WEIGHT;

	  const std::vector<LabelEntry> &s1 = dir ? this->out_labels[v] : this->in_labels[v];
	  const std::vector<LabelEntry> &s2 = this->in_labels[w];
	  size_t i1 = 0, i2 = 0;

	  custom_node max_index = std::max(node_to_index(v),node_to_index(w));

	  while (i1 < s1.size() && i2 < s2.size()) {
		  if(s1[i1].v > max_index || s2[i2].v > max_index)
			break;

		  else if (s1[i1].v < s2[i2].v)
			  ++i1;
		  else if (s1[i1].v > s2[i2].v)
			  ++i2;
		  else{
			  custom_node vertex = s1[i1].v;
			  assert(s1[i1].v == s2[i2].v);

			  if (vertex == NULL_NODE)
				  break;


			  if(s1[i1].d + s2[i2].d < distance){
				  distance = s1[i1].d + s2[i2].d;
				  assert(s1[i1].v == s2[i2].v);
				  hubset.clear();
				  hubset.push_back(s1[i1].v);
			  }
			  else if(s1[i1].d + s2[i2].d == distance){
				  assert(s1[i1].v == s2[i2].v);
				  hubset.push_back(s1[i1].v);
			  }

			  for (++i1; s1[i1].v == vertex; ++i1){}
			  for (++i2; s2[i2].v == vertex; ++i2){}
		  }
	  }
	  assert(v!=w || hubset.size()==1);
	  assert(std::is_sorted(hubset.begin(),hubset.end()));
	  return hubset;
}

std::vector<custom_node> Labeling::hubs(custom_node v, custom_node w) {

	  std::vector<custom_node> hubset;
	  assert(!dir || this->in_labels.size() == this->out_labels.size());
	  assert(v<=this->in_labels.size() && w<=this->in_labels.size());

	  if(v==w){
		  hubset.clear();
		  hubset.push_back(node_to_index(v));
		  return hubset;
	  }
	  hubset.clear();

	  custom_weight distance = NULL_WEIGHT;


	  const std::vector<LabelEntry> &s1 = dir ? this->out_labels[v] : this->in_labels[v];
	  const std::vector<LabelEntry> &s2 = this->in_labels[w];
      int previous = -1;

	  size_t i1 = 0, i2 = 0;
	  custom_node max_index = std::max(node_to_index(v),node_to_index(w));

	  while (i1 < s1.size() && i2 < s2.size()) {

		  if(s1[i1].v > max_index || s2[i2].v > max_index)
			  break;
		  else if (s1[i1].v < s2[i2].v)
			  ++i1;
		  else if (s1[i1].v > s2[i2].v)
			  ++i2;
		  else{
			  custom_node vertex = s1[i1].v;
			  assert(s1[i1].v == s2[i2].v);

			  if (vertex == NULL_NODE)
				  break;
			  previous = (int) s1[i1].v;

			  if(s1[i1].d + s2[i2].d < distance){
				  distance = s1[i1].d + s2[i2].d;
				  assert(s1[i1].v == s2[i2].v);
				  hubset.clear();
				  hubset.push_back(s1[i1].v);
			  }
			  else if(s1[i1].d + s2[i2].d == distance){
				  assert(s1[i1].v == s2[i2].v);
				  hubset.push_back(s1[i1].v);
			  }
			  for (++i1; s1[i1].v == vertex; ++i1){}
			  for (++i2; s2[i2].v == vertex; ++i2){}
		  }
	  }
	  assert(v!=w || hubset.size()==1);
	  return hubset;
}




custom_node Labeling::getNumberOfLabelEntries() {
	custom_node num_labels = 0;
	for(size_t t = 0;t<in_labels.size();t++)
		num_labels += in_labels[t].size()-1;

	if(!dir)
		if(out_labels.size()!=0)
			throw new std::runtime_error("non empty out labels");

	for(size_t t = 0;t<out_labels.size();t++)
		num_labels += out_labels[t].size()-1;

	return num_labels;
}

