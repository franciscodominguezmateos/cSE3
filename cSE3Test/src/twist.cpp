/*
 * twist.cpp
 *
 *  Created on: 2 Dec 2020
 *      Author: Francisco Dominguez
 *  This file can or can not be used.
 *  If you don't put it in your project everything works fine but this function
 */
#include "pose.h"
Pose Twist::exp(){
	Twist &self=*this;
	return Pose::exp(self);
}




