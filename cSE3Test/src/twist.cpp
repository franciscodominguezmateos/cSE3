/*
 * twist.cpp
 *
 *  Created on: 2 Dec 2020
 *      Author: Francisco Dominguez
 */
#include "pose.h"
Pose Twist::exp(){
	Twist &self=*this;
	return Pose::exp(self);
}




