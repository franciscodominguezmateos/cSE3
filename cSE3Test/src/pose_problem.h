/*
 * pose_problem.h
 *
 *  Created on: 4 Dec 2020
 *      Author: Francisco Dominguez
 */

#pragma once
#include "pose.h"
class PoseProblem{
public:
	virtual int size()     =0;
	virtual Mat gradEls()  =0;
	virtual double  Els()  =0;
	virtual Pose &getPose()=0;
	virtual Mat res(int i) =0;
	virtual Mat Jres(int i)=0;
};



