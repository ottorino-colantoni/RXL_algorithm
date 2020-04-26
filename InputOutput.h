#include <iostream>
#include <fstream>
#include "progressBar.h"
#include "Labeling.h"
#include "LabelEntry.h"
#include <cassert>
#include <vector>
#include "matplotlibcpp.h"
#include <array>


#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H
#define  LOGS_LOCATION "./LogFiles/"
#define  UWGRAPHS_LOCATION "./uwGraphs/"
#define  GRAPHS_LOCATION  "./Graphs/"

class InputOutput {

public:

    InputOutput();

    bool readLabelsFromFile(Labeling* labels, std::string file_name);

    bool printLabelsOnFile(Labeling* labels, std::string file_name);

    bool printLogCompare(std::vector<std::vector<float>> data, std::string graph_name, std::string file_name);

    bool printPlot(std::vector<int> xs, std::vector<std::vector<float>> ys, std::string title, std::string xlabel, std::string ylabel, std::string filename);

    bool printTestResult(std::vector<std::vector<float>> &data, std::string file_name);

    bool readTestResult(std::vector<std::vector<float>> &data, std::vector<int> &rounds, std::string file_location);

	bool changeGraphWeight(std::string graph_location);

    bool networkRepoGraph(std::string file_name);
};


#endif //RXL_ALGORITHM_INPUTOUTPUT_H
