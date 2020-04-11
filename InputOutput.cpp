
#include "InputOutput.h"
#include <string>

InputOutput::InputOutput() {}


bool InputOutput::printLabelsOnFile(Labeling* labels, std::string file_location) {

    std::ofstream myfile;
    try{
        myfile.open(file_location, std::ios::trunc);
        int size = 32000;
        myfile<<size<<'\n';
        for (int i = 0; i < labels->in_labels.size() ; i++) {
            for (int j = 0; j <labels->in_labels[i].size() ; j++) {
                if(labels->in_labels[i][j].d != NULL_WEIGHT) {
                    myfile << i << ' ' << labels->in_labels[i][j].v << ' ' << labels->in_labels[i][j].d <<'\n';
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

bool InputOutput::readLabelsFromFile(Labeling *labels, std::string file_location) {

    std::ifstream ifs(file_location);

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