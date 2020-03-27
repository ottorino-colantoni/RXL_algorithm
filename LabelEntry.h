/*
 * LabelEntry.h
 *
 *  Created on: 05 feb 2016
 *      Author: anonym
 */

#ifndef LABELENTRY_H_
#define LABELENTRY_H_
#include "CustomDataTypes.h"



class LabelEntry {
public:
	LabelEntry(custom_node,custom_weight);
	LabelEntry();

	custom_node v;
	custom_weight d;
	bool operator<(LabelEntry other) const
	{
		return v > other.v || (v == other.v && d > other.d);
	}

	virtual ~LabelEntry();
};



#endif /* LABELENTRY_H_ */
