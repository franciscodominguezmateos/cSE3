/*
 * twist.cpp
 *
 *  Created on: 2 Dec 2020
 *      Author: Francisco Dominguez
 */
#include "pose.h"
Pose Twist::operator*(Pose p){
	return this->exp()*p;
}
Pose Twist::exp(){
	return Pose::exp(*this);
}





