
#include "InputOutput.h"
#include <string>

InputOutput::InputOutput() {}


bool InputOutput::printLabelsOnFile(Labeling* labels, std::string file_location) {

    std::ofstream myfile;
    myfile.exceptions(std::ofstream::badbit);
    try{
        myfile.open(file_location);
        myfile<<labels->in_labels.size()<<"\n";
        for (int i = 0; i < labels->in_labels.size() ; i++) {
            for (int j = 0; j <labels->in_labels[i].size() ; ++j) {
                if(labels->in_labels[i][j].d != NULL_WEIGHT) {
                    myfile << i << " " << labels->in_labels[i][j].v << " " << labels->in_labels[i][j].d << "\n";
                }
                else{
                    continue;
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

    int prec_node = -1;

    std::ifstream ifs(file_location);

    if (!ifs)
        throw std::runtime_error("Error opening File ");

    int size = -1, node = -1, hub = -1, distance = -1;
    ifs >> size;
    assert(size>0);
    labels->in_labels.resize(size);
    ProgressStream reader(size);

    while(true){


        ifs >> node >> hub >> distance;

        if(ifs.eof())
            break;

        assert(node >=0 && hub >=0 && distance >=0);

        labels->in_labels[node].push_back(LabelEntry(hub, distance));

        if(prec_node != node){
            ++reader;
            prec_node = node;
        }
    }

    ifs.close();


}