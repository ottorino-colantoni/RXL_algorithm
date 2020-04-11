#include <iostream>
#include <fstream>
#include "progressBar.h"
#include "Labeling.h"
#include "LabelEntry.h"
#include <cassert>
#include <vector>

#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

class InputOutput {

public:

    InputOutput();

    bool readLabelsFromFile(Labeling* labels, std::string file_location);

    bool printLabelsOnFile(Labeling* labels, std::string file_location);


};


#endif //RXL_ALGORITHM_INPUTOUTPUT_H
