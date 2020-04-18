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
#include <vector>
#include "matplotlibcpp.h"




void plotResult(std::vector<std::vector<float>> data, std::vector<int> rounds){

    matplotlibcpp::plot(rounds, data[4]);
    matplotlibcpp::title("Create forest time trend");
    matplotlibcpp::save("./createForestTime.png");
    matplotlibcpp::show();

    matplotlibcpp::plot(rounds, data[5]);
    matplotlibcpp::title("Total time trend");
    matplotlibcpp::save("./totalTime.png");
    matplotlibcpp::show();

    matplotlibcpp::plot(rounds, data[6]);
    matplotlibcpp::title("Update time trend");
    matplotlibcpp::save("./updateTime.png");
    matplotlibcpp::show();

    matplotlibcpp::plot(rounds, data[7]);
    matplotlibcpp::title("Encrease forest average time trend");
    matplotlibcpp::save("./encreaseForestAVGtime.png");
    matplotlibcpp::show();

    matplotlibcpp::plot(rounds, data[8]);
    matplotlibcpp::title("Number of labels trend");
    matplotlibcpp::save("./numberOfLabels.png");
    matplotlibcpp::show();

}


void runRXL(std::string graph_location,int num_samples,int num_counters,int num_newsamples,int max_numtrees,std::string output_location){

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


void testRXL(std::string graph_location,std::vector<int> num_samples,std::vector<int> num_counters,std::vector<int> num_newsamples,std::vector<int> max_numtrees,std::string output_location){


        /*
         * data:
         * - 0 num_samples;
         * - 1 num_counters;
         * - 2 num_newsamples;
         * - 3 max_numtrees;
         * - 4 create forest time;
         * - 5 total time;
         * - 6 update avg time;
         * - 7 encrease avg time;
         * - 8 number of labels
         */

		NetworKit::Graph *graph;
		Auxiliary::read(graph_location, false, &graph);
		int round = 0;
		std::vector<std::vector<float>> data;
		data.resize(9);
		std::vector<int> rounds;

		for(int i=0;i<num_samples.size();i++){
			
			for(int j=0;j<num_counters.size();j++){
			
				for(int z=0;z<num_newsamples.size();z++){
				
					for(int k=0; k<max_numtrees.size();k++){
					
						data[0].push_back(num_samples[i]);
						data[1].push_back(num_counters[j]);
						data[2].push_back(num_newsamples[z]);
						data[3].push_back(max_numtrees[k]);
						
						int max;
						int avgup=0;
						int avgen=0;
						mytimer timer;
						mytimer timerup;
						float avgupdate=0;
						float avgencrease=0;
						std::pair <std::vector<custom_node>, std::vector<custom_node>> keeper;
						keeper.second.resize(graph->upperNodeIdBound());
						Labeling *labeling = new Labeling(graph->isDirected());
						Labeling_Tools *lt = new Labeling_Tools(graph, labeling, keeper);
						SamPG *spg = new SamPG(num_samples[i], num_counters[j], graph);
						timer.restart();
						spg->createForest();
						data[4].push_back(timer.elapsed());
						for (int i = 0; i < graph->numberOfNodes(); i++) {
							max = spg->maxDescNode();
							lt->add_node_to_keeper(max, i);
							lt->weighted_build_RXL();
							if(!spg->isEnded()) {
								timerup.restart();
								spg->updateForest(max);
								avgup++;
								avgupdate+=timerup.elapsed();
								if (spg->getNumSamples()< max_numtrees[k]) {
									timerup.restart();
									spg->encreaseForest(num_newsamples[z], labeling);
									avgencrease+= timerup.elapsed();
									avgen++;
								}
								
							}
						}
						data[5].push_back(timer.elapsed());
						data[6].push_back(avgupdate/avgup);
						data[7].push_back(avgencrease/avgen);
						data[8].push_back(labeling->getNumberOfLabelEntries());
						rounds.push_back(round);
						round++;
						delete spg, labeling;
					}
				
				}
			
			}
		
		}

		plotResult(data, rounds);

}



int main(int argc, char** argv){

    //declare supported options
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");

    desc.add_options()
            ("graph_location,g", po::value<std::string>(), "Input Graph File Location")
            ("num_samples,k", po::value<std::vector<int>>()->multitoken(), "Number of samples to use at the beginning")
            ("num_counters,c", po::value<std::vector<int>>()->multitoken(), "Number of counters to use")
            ("num_newsamples,n", po::value<std::vector<int>>()->multitoken(), "Number of new samples introduced at each iteration")
            ("max_numtrees,m",po::value<std::vector<int>>()->multitoken(), "Max number of trees used as samples")
            ("output_location,o", po::value<std::string>(), "Location where to save the output")
            ("test,t", po::value<int>(),"execution mode : {test:1 or run:0}")
            ;


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    std::vector<int> num_samples, num_counters , num_newsamples , max_numtrees;
    std::string graph_location = "", output_location = "";
	int test= -1;
    if (vm.empty()){
        std::cout << desc << "\n";
        throw std::runtime_error("Empty options");
    }
    if (vm.count("num_samples"))
        num_samples = vm["num_samples"].as<std::vector<int>>();

    if (vm.count("num_counters"))
        num_counters = vm["num_counters"].as<std::vector<int>>();

    if (vm.count("num_newsamples"))
        num_newsamples = vm["num_newsamples"].as<std::vector<int>>();

    if (vm.count("max_numtrees"))
        max_numtrees = vm["max_numtrees"].as<std::vector<int>>();

    if (vm.count("graph_location"))
        graph_location = vm["graph_location"].as<std::string>();

    if (vm.count("output_location"))
        output_location = vm["output_location"].as<std::string>();
        
    if (vm.count("test"))
    test = vm["test"].as<int>();


	for(int i=0; i<num_samples.size();i++){
		if(num_samples[i]<=0){
		    std::cout << desc << "\n";
		    throw std::runtime_error("Number of samples must be greater than 0");
		}
    }

    if(graph_location == ""){
        std::cout << desc << "\n";
        throw std::runtime_error("Wrong graph_location");
    }

	for(int i=0; i<num_counters.size();i++){
		if(num_counters[i] <= 0 || num_counters > num_samples){
		    std::cout << desc << "\n";
		    throw std::runtime_error("Number of counters must fall in [1, num_samples-1]");
		}
    }

	for(int i=0; i<max_numtrees.size();i++){
		if(max_numtrees[i] < 0){
		    std::cout << desc << "\n";
		    throw std::runtime_error("Max number of trees must fall in [1, graph size]");
		}
    }
    

	if(test == 0){
		runRXL(graph_location,num_samples[0],num_counters[0],num_newsamples[0],max_numtrees[0],output_location);
    }
    
    else if(test == 1){
        testRXL(graph_location, num_samples, num_counters, num_newsamples, max_numtrees, output_location);
    }
    
    else{throw std::runtime_error("Bad value (test) : choose 1 for RXL benchmark or 0 to print labels on file");}
}







