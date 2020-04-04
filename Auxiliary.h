/*
 * Labeling.h
 *
 *  Created on: 01 feb 2016
 *      Author: Mattia D'Emidio
 */

#ifndef AUXILIARY_H_
#define AUXILIARY_H_

#include "CustomDataTypes.h"
#include </usr/include/boost/filesystem.hpp>
#include <dirent.h>
#include "Labeling.h"
#include "callgrind.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "CustomDataTypes.h"
class Auxiliary {

public:
	Auxiliary();
	virtual ~Auxiliary();
	static std::pair<std::vector<custom_node>,std::vector<custom_node>> compute_ordering(NetworKit::Graph*,int);

	static void test_match(NetworKit::Graph*,Labeling*,int);
	static std::string getTimeString();
	static void convert(std::string,unsigned int type=0,unsigned int pream_length=0, bool weighted = false, bool directed = false);
	static void real_world_parser(std::string, bool weighted, bool directed,custom_node,custom_node);
	static std::string exec(const char*);
	static void process_data(std::string);
	static void box_plot_column(std::string,int);

	static bool exists (const std::string&);
	static void init_mods(std::vector<NetworKit::GraphEvent>&, int, int,NetworKit::Graph*,double=0);
	static void read_mods(std::string, std::vector<NetworKit::GraphEvent>&);

	static void copy(NetworKit::Graph**,NetworKit::Graph**);

	static void store_graph(std::string,NetworKit::Graph**);
	static void store_mods(std::string,std::vector<NetworKit::GraphEvent>&);

	static void read(std::string,bool,NetworKit::Graph**);
	static void generate(NetworKit::Graph**,int,bool=false);

	static double avg(std::vector<double>&);
	static double avg(std::vector<int>&);
	static double stddev(std::vector<double>&);
	static double stddev(std::vector<int>&);
	static void SCC(NetworKit::Graph*);
	static bool isInteger(const std::string &);
	static void get_all(const boost::filesystem::path&, const std::string& , std::vector<boost::filesystem::path>&);
	static double getQuartile(int, std::vector<double>&);

};

#endif /* AUXILIARY_H_ */
