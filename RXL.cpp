//
// Created by f4b3r on 03/04/20.
//


#include "Labeling.h"
#include "Auxiliary.h"
#include "CustomDataTypes.h"
#include "Labeling_Tools.h"
#include "Dijkstra.h"
#include "SamPG.h"
#include "networkit/graph/Graph.hpp"
#include <stdio.h>
#include "mytimer.h"
#include <boost/program_options.hpp>
#include <omp.h>
#include "progressBar.h"
#include "InputOutput.h"



int main(int argc, char** argv){

    //declare supported options
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");

    desc.add_options()
            ("graph_location,g", po::value<std::string>(), "Input Graph File Location")
            ("num_samples,k", po::value<int>(), "Number of samples to use at the beginning")
            ("num_counters,c", po::value<int>(), "Number of counters to use")
            ("num_newsamples,n", po::value<int>(), "Number of new samples introduced at each iteration")
            ("max_numtrees,m",po::value<int>(), "Max number of trees used as samples")
            ("output_location,o", po::value<std::string>(), "Location where to save the output")
            ;


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    int num_samples = -1, num_counters = -1, num_newsamples = -1, max_numtrees = -1;
    std::string graph_location = "", output_location = "";

    if (vm.empty()){
        std::cout << desc << "\n";
        throw std::runtime_error("Empty options");
    }
    if (vm.count("num_samples"))
        num_samples = vm["num_samples"].as<int>();

    if (vm.count("num_counters"))
        num_counters = vm["num_counters"].as<int>();

    if (vm.count("num_newsamples"))
        num_newsamples = vm["num_newsamples"].as<int>();

    if (vm.count("max_numtrees"))
        max_numtrees = vm["max_numtrees"].as<int>();

    if (vm.count("graph_location"))
        graph_location = vm["graph_location"].as<std::string>();

    if (vm.count("output_location"))
        output_location = vm["output_location"].as<std::string>();



    if(num_samples<=0){
        std::cout << desc << "\n";
        throw std::runtime_error("Number of samples must be greater than 0");
    }

    if(graph_location == ""){
        std::cout << desc << "\n";
        throw std::runtime_error("Wrong graph_location");
    }

    if(num_counters <= 0 || num_counters > num_samples){
        std::cout << desc << "\n";
        throw std::runtime_error("Number of counters must fall in [1, num_samples-1]");
    }

    if(max_numtrees < num_samples || max_numtrees < 0){
        std::cout << desc << "\n";
        throw std::runtime_error("Max number of trees must fall in [1, graph size]");
    }


    NetworKit::Graph *graph;
    Auxiliary::read(graph_location, false, &graph);
    SamPG *spg = new SamPG(num_samples, num_counters, graph);
    spg->createForest();
	int max;

    std::pair <std::vector<custom_node>, std::vector<custom_node>> keeper;
    keeper.second.resize(graph->upperNodeIdBound());
    Labeling *labeling = new Labeling(graph->isDirected());
    Labeling_Tools *lt = new Labeling_Tools(graph, labeling, keeper);

    ProgressStream builder_(graph->numberOfNodes());
    builder_.label() << "Building WEIGHTED UNDIRECTED labeling for " <<graph->numberOfNodes()<< " vertices";

    for (int i = 0; i < graph->numberOfNodes(); i++) {
        max = spg->maxDescNode();
        lt->add_node_to_keeper(max, i);
        lt->weighted_build_RXL();
        if(!spg->isEnded()) {
            spg->updateForest(max);
            if (spg->getNumSamples()< max_numtrees) {
                spg->encreaseForest(num_newsamples, labeling);
            }
        }
        ++builder_;
    }


    InputOutput* io = new InputOutput();
    io->printLabelsOnFile(labeling, output_location);

}

