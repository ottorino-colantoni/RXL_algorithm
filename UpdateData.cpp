/*
 * UpdateData.cpp
 *
 *  Created on: 05 feb 2016
 *      Author: anonym
 */

#include "UpdateData.h"

UpdateData::UpdateData() {
	this->time = 0.0;
	this->affected = 0;
	this->affected_x = 0;
	this->affected_y = 0;
	this->added = 0;
	this->removed = 0;

}

UpdateData::UpdateData(double t, int a,int b,int c, int d,int e) {
	this->time = t;
	this->affected_x = a;
	this->affected_y = b;
	this->affected = c;

	this->added = d;
	this->removed = e;

}


UpdateData::~UpdateData() {}


