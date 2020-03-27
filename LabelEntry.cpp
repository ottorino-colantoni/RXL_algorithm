/*
 * LabelEntry.cpp
 *
 *  Created on: 05 feb 2016
 *      Author: anonym
 */

#include "LabelEntry.h"

LabelEntry::LabelEntry(custom_node nd,custom_weight wg){
	this->v = nd;
	this->d = wg;

}
LabelEntry::LabelEntry() {
	this->v = NULL_NODE;
	this->d = NULL_WEIGHT;

}

LabelEntry::~LabelEntry() {

}


