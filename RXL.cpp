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

    InputOutput* io = new InputOutput();
    io->printPlot(rounds, {data[4]}, "Create Forest time trend", "iteration", "time (s)", "createForestTime.png");
    io->printPlot(rounds, {data[5]}, "Total time trend", "iteration", "time (s)", "totalTime.png");
    io->printPlot(rounds, {data[6]}, "Update time trend", "iteration", "time (s)", "updateTime.png");
    io->printPlot(rounds, {data[7]}, "Encrease forest time trend", "iteration", "time (s)", "encreaseForestTime.png");
    io->printPlot(rounds, {data[8]}, "Number Of Labels trend", "iteration", "num. labels", "numOfLabels.png");

}


void randomQueryTest(Labeling *labels1, Labeling *labels2, int graph_size, std::vector<std::vector<float>> &data){

    mytimer timer;
    float avg1=0;
    float avg2=0;
    BestRandom *random = new BestRandom(0, graph_size);
    std::vector<int> attempts;
    std::vector<float> time1;
    std::vector<float> time2;
    int num_attempts = graph_size/10;
    attempts.resize(num_attempts);
    time1.resize(num_attempts);
    time2.resize(num_attempts);
    int value1;
    int value2;
    int dist1;
    int dist2;
    int corrects = 0;

    for (int i = 0; i < num_attempts; i++) {
        value1 = random->random();
        value2 = random->random();
        timer.restart();
        dist1 = labels1->query(value1,value2);
        time1[i] = timer.elapsed();
        avg1 += time1[i];
        timer.restart();
        dist2 = labels2->query(value1,value2);
        time2[i] = timer.elapsed();
        avg2 += time2[i];
        if(dist1 == dist2){
            corrects += 1;
        }
        attempts[i] = i;
    }

    avg1 /= num_attempts;
    avg2 /= num_attempts;

    std::vector<std::vector<float>> y{time1, time2};

    InputOutput* io = new InputOutput();
    io->printPlot(attempts, y, "Query time trend", "iteration", "time(s)", "queryTimeTrend.png");

    data[0].push_back(avg1);
    data[1].push_back(avg2);
    data[0].push_back(corrects*100/num_attempts);
    data[1].push_back(corrects*100/num_attempts);
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
		    }
            if (spg->getNumSamples()< max_numtrees) {
                spg->encreaseForest(num_newsamples, labeling);
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
                    }
				}
			
			}
		
		}

        InputOutput* io = new InputOutput();
		io->printTestResult(data, output_location);
}

void compare(std::string graph_location,int num_samples,int num_counters,int num_newsamples,int max_numtrees,std::string output_location){

    NetworKit::Graph *graph;
    Auxiliary::read(graph_location, false, &graph);
    mytimer timer;
    std::vector<std::vector<float>> data;
    data.resize(2);

    /* data structure
     *
     * - 0 PLL
     * - 1 RXL
     *
     * - 0 preprocessing time;
     * - 1 labels' number;
     * - 2 avg time;
     * - 3 correct queries percentage.
     */

    //Creating Labels with PLL
    timer.restart();
    Labeling *PLLlabels = new Labeling(graph->isDirected());
    std::pair <std::vector<custom_node>, std::vector<custom_node>> keeperPLL;
    keeperPLL = Auxiliary::compute_ordering(graph, 0);
    Labeling_Tools *ltPLL = new Labeling_Tools(graph, PLLlabels, keeperPLL, false);
    data[0].push_back(timer.elapsed());



    //Creating Labels with RXL
    timer.restart();
    SamPG *spg = new SamPG(num_samples, num_counters, graph);
    spg->createForest();
    int max;
    std::pair <std::vector<custom_node>, std::vector<custom_node>> keeperRXL;
    keeperRXL.second.resize(graph->upperNodeIdBound());
    Labeling *RXLlabels = new Labeling(graph->isDirected());
    Labeling_Tools *lt = new Labeling_Tools(graph, RXLlabels, keeperRXL);

    for (int i = 0; i < graph->numberOfNodes(); i++) {
        max = spg->maxDescNode();
        lt->add_node_to_keeper(max, i);
        lt->weighted_build_RXL();
        if(!spg->isEnded()) {
            spg->updateForest(max);
            if (spg->getNumSamples()< max_numtrees) {
                spg->encreaseForest(num_newsamples, RXLlabels);
            }
        }
    }
    data[1].push_back(timer.elapsed());
    data[0].push_back(PLLlabels->getNumberOfLabelEntries());
    data[1].push_back(RXLlabels->getNumberOfLabelEntries());
    randomQueryTest(PLLlabels, RXLlabels, graph->numberOfNodes(), data);

    InputOutput* io = new InputOutput();
    io->printLogCompare(data, graph_location, output_location);


}


void generateTestsPlots(std::string output_location){
    InputOutput *io = new InputOutput();
    std::vector<std::vector<float>> data;
    std::vector<int> rounds;
    io->readTestResult(data, rounds, output_location);
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
            ("exe,e", po::value<int>(),"execution mode : {run:0, test:1, compare:2, plot:3}")
            ;


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    std::vector<int> num_samples, num_counters , num_newsamples , max_numtrees;
    std::string graph_location = "", output_location = "";
	int exe= -1;
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
        
    if (vm.count("exe"))
    exe = vm["exe"].as<int>();

    if(exe != 3) {
        for (int i = 0; i < num_samples.size(); i++) {
            if (num_samples[i] <= 0) {
                std::cout << desc << "\n";
                throw std::runtime_error("Number of samples must be greater than 0");
            }
        }

        if (graph_location == "") {
            std::cout << desc << "\n";
            throw std::runtime_error("Wrong graph_location");
        }

        for (int i = 0; i < num_counters.size(); i++) {
            if (num_counters[i] <= 0 || num_counters > num_samples) {
                std::cout << desc << "\n";
                throw std::runtime_error("Number of counters must fall in [1, num_samples-1]");
            }
        }

        for (int i = 0; i < max_numtrees.size(); i++) {
            if (max_numtrees[i] < 0) {
                std::cout << desc << "\n";
                throw std::runtime_error("Max number of trees must fall in [1, graph size]");
            }
        }
    }

	if(exe == 0){
		runRXL(graph_location,num_samples[0],num_counters[0],num_newsamples[0],max_numtrees[0],output_location);
    }
    
    else if(exe == 1){
        testRXL(graph_location, num_samples, num_counters, num_newsamples, max_numtrees, output_location);
    }

    else if(exe == 2){
        compare(graph_location, num_samples[0], num_counters[0], num_newsamples[0], max_numtrees[0], output_location);
    }

    else if(exe == 3){
        generateTestsPlots(output_location);
    }
    else if(exe == 4){
    	InputOutput* io= new InputOutput();

   		io->changeGraphWeight(graph_location);
    	
    }
    
    else{throw std::runtime_error("Bad value (test) : choose 1 for RXL benchmark or 0 to print labels on file");}
}







