/*
 * UpdateData.h
 *
 *  Created on: 05 feb 2016
 *      Author: anonym
 */

#ifndef UPDATEDATA_H_
#define UPDATEDATA_H_


class UpdateData {
public:
	UpdateData();
	UpdateData(double,int,int,int,int,int);
	double time;
	int affected;
	int affected_x;
	int affected_y;
	int added;
	int removed;
	virtual ~UpdateData();
};



#endif /* UPDATEDATA_H_ */
