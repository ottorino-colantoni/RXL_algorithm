/*
 * CustomDataTypes.h
 *
 *  Created on: 01 feb 2016
 *      Author: anonym
 */
#include <limits>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <iomanip>
#include <boost/container/vector.hpp>
#include <string>
#include "mersenneTwister.h"
#include "progressBar.h"
#include "UpdateData.h"
#include <NetworKit/io/KONECTGraphReader.h>
#include <NetworKit/io/SNAPGraphReader.h>
#include <NetworKit/dynamics/GraphEvent.h>
#include <NetworKit/dynamics/GraphUpdater.h>
#include "UpdateData.h"
#include <boost/program_options.hpp>
#include <NetworKit/distance/Dijkstra.h>
#include <NetworKit/distance/BFS.h>
#include <NetworKit/graph/Graph.h>
#include <boost/lexical_cast.hpp>

#ifndef CUSTOMDATATYPES_H_
#define CUSTOMDATATYPES_H_

typedef NetworKit::node custom_node;
typedef NetworKit::edgeweight custom_weight;
typedef unsigned int custom_time;
typedef NetworKit::count custom_count;
constexpr custom_node NULL_NODE = std::numeric_limits<custom_node>::max();
constexpr custom_weight NULL_WEIGHT = std::numeric_limits<custom_weight>::max();
constexpr custom_node NULL_INDEX = std::numeric_limits<custom_node>::max()-1000000000000;
constexpr custom_time EDGE_TIME_NULL = 0;
constexpr custom_time NULL_TIME = std::numeric_limits<custom_time>::max();
constexpr custom_count NULL_COUNT = std::numeric_limits<custom_count>::max();




struct DEDGE {
	unsigned int source		: 31;
	bool forward		:  1;
	unsigned int target		: 31;
	bool backward		:  1;
	//EdgeWeight
	int weight;

	// constructors
	DEDGE(){}
	DEDGE(unsigned int s, unsigned int t, unsigned long long w, bool f, bool b) :
		source(s), forward(f), target(t), backward(b), weight(w) {}
};

template<typename Primitive>
void readPrimitive(std::istream& in, Primitive& p){
  in.read((char*)&p, sizeof(p));
}

template<typename Primitive>
void writePrimitive(std::fstream& f, Primitive& p){
  f.write((char*)&p, sizeof(p));
}


#endif /* CUSTOMDATATYPES_H_ */
