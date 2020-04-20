
#include "InputOutput.h"
#include <string>
#include <cassert>


InputOutput::InputOutput() {}

bool InputOutput::printPlot(std::vector<int> xs, std::vector<std::vector<float>> ys, std::string title, std::string xlabel, std::string ylabel, std::string filename) {
    for (int i = 0; i < ys.size() ; i++) {
        assert(xs.size() == ys[i].size());
        matplotlibcpp::plot(xs, ys[i]);
    }
    matplotlibcpp::xlabel(xlabel);
    matplotlibcpp::ylabel(ylabel);
    matplotlibcpp::title(title);
    matplotlibcpp::save(LOGS_LOCATION + filename);
    matplotlibcpp::close();
}

bool InputOutput::printLabelsOnFile(Labeling* labels, std::string file_name) {

    std::ofstream myfile;
    try{
        myfile.open(LOGS_LOCATION + file_name, std::ios::trunc);
        int size = labels->in_labels.size();
        myfile<<size<<'\n';
        for (int i = 0; i < labels->in_labels.size() ; i++) {
            for (int j = 0; j <labels->in_labels[i].size() ; j++) {
                if(labels->in_labels[i][j].d != NULL_WEIGHT) {
                    myfile << i << ' ' << labels->index_to_node(labels->in_labels[i][j].v) << ' ' << labels->in_labels[i][j].d <<'\n';
                }
            }
        }
    }
    catch (const std::ofstream::failure& e){
        std::cout<< "Exception opening file";
        return false;
    }

    myfile.close();

    return true;

}

bool InputOutput::readLabelsFromFile(Labeling *labels, std::string file_name) {

    std::ifstream ifs(LOGS_LOCATION + file_name);

    if (!ifs)
        throw std::runtime_error("Error opening File ");

    int size, node, hub, distance;
    ifs >>size;
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    assert(size>0);
    labels->in_labels.resize(size);
    while(true){


        ifs >> node >> hub >> distance;
        ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if(ifs.eof())
            break;
        assert(node >=0 && hub >=0 && distance >=0);

       labels->in_labels[node].push_back(LabelEntry(hub, distance));

    }

    ifs.close();

}

bool InputOutput::printLogCompare(std::vector<std::vector<float>> data, std::string graph_name,std::string file_name){

    std::ofstream myfile;
    try {
        myfile.open(LOGS_LOCATION + file_name, std::ios::trunc);
        myfile <<"Graph name: "<<graph_name << "\n";
        for (int i = 0; i < data.size(); i++) {
            if (i == 0) {
                myfile << "PLL data" << "\n";

            } else if (i == 1) {
                myfile << "RXL data" << "\n";
            }
            myfile << "Preprocessing time: " << data[i][0] << "  Number of Labels: " << data[i][1]
                   << "  Query average time (sec): " << data[i][2] << "  Correctness percentage: " << data[i][3]
                   << "\n";
        }
    }
    catch (const std::ofstream::failure& e){
        std::cout<< "Exception opening file";
        return false;
    }

    myfile.close();

    return true;

}