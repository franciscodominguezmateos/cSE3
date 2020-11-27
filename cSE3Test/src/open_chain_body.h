/*
 * open_chain_body.h
 *
 *  Created on: 27 Nov 2020
 *      Author: Francisco Dominguez
 */
#pragma once
#include <vector>
#include "pose.h"

class OpenChainBody{
	vector<ScrewAxis> joints;
	//jacobian columns
	vector<Twist> J;
	Pose home; //"zero" position of the robot
public:
	OpenChainBody(vector<ScrewAxis> tw,Pose h):joints(tw),home(h){
		J.resize(tw.size());
	}
	Pose forwardUpTo(int N,vector<double> &states){
		Pose p;
		for(int i=0;i<N;i++){
			ScrewAxis &sa=joints[i];
			double state=states[i];
			p=p+sa*state;
		}
		return p;
	}
	Pose forward(vector<double> &states){
		Pose p=forwardUpTo(joints.size(),states);
		p=p+home;
		return p;
	}
	void buildJacobianColumns(vector<double> &states){
		for(unsigned int i=0;i<states.size();i++){
			Pose p=forwardUpTo(i,states);
			J[i]=p.Ad(joints[i]);
		}
	}
	Mat jacobian(vector<double> &states){
		buildJacobianColumns(states);
		Mat m=J[0].asMat();
		for(unsigned int i=1;i<J.size();i++){
			Mat mj=J[i].asMat();
			hconcat(m,mj,m);
		}
		return m;
	}
};
