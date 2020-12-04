/*
 * twist.cpp
 *
 *  Created on: 2 Dec 2020
 *      Author: Francisco Dominguez
 */
#include "pose.h"
//Left composition
/*Pose Twist::operator+(Pose p){
	Pose pSelf=this->exp();
	return pSelf+p;
}*/
Pose Twist::exp(){
	Twist &self=*this;
	return Pose::exp(self);
}




