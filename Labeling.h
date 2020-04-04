/*
 * Labeling.h
 *
 *  Created on: 01 feb 2016
 *      Author: Mattia D'Emidio
 */

#ifndef LABELING_H_
#define LABELING_H_

#include "CustomDataTypes.h"
#include "LabelEntry.h"
class Labeling {

public:


	//MAIN DATA STRUCTURES
	std::vector<std::vector<LabelEntry>> in_labels;
	std::vector<std::vector<LabelEntry>> out_labels;

	std::vector<custom_node> order;
	std::vector<custom_node> reverse_order;
	custom_node index_to_node(custom_node);
	custom_node node_to_index(custom_node);
	bool dir;
	Labeling(bool);

	virtual ~Labeling();


	custom_weight query(custom_node,custom_node);
	void query(custom_node,custom_node,custom_weight,std::pair<custom_node,custom_weight>&);

//	std::vector<custom_node> shared_vertices(custom_node,custom_node);

	std::vector<custom_node> hubs(custom_node,custom_node);
	std::vector<custom_node> hubs(custom_node,custom_node,custom_weight&);
//	std::vector<custom_node> hubs(custom_weight,custom_node v, custom_node w);

//	int hubs_size(custom_node,custom_node);



//	custom_node getLabelingSize(); //returns bytes
//	double getAverageLabelSizePerNode(); //return bytes
	custom_node getNumberOfLabelEntries();
//	double getAverageNumberOfLabelsPerNode();

	void printInLabels();

};

#endif /* LABELING_H_ */
