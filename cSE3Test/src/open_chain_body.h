/*
 * open_chain_body.h
 *
 *  Created on: 27 Nov 2020
 *      Author: Francisco Dominguez
 */
#pragma once
#include <vector>
#include "pose.h"

class GeneralizedCoord : public Mat{

};
class OpenChainBody{
	vector<ScrewAxis> joints;
	//jacobian columns
	vector<Twist> J;  //Space/stationary/world/global Jacobian
	vector<Twist> Jb; //Body Jacobian
	Pose home; //"zero" position of the robot
public:
	OpenChainBody(vector<ScrewAxis> tw,Pose h):joints(tw),home(h){
		J.resize(tw.size());
	}
	/*  S P A C E  */
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
	vector<double> ikStep(vector<double> &states,Twist goal){
		Mat j=jacobian(states);
		Mat Jpinv;
		invert(j, Jpinv, DECOMP_SVD);
		Mat dState=Jpinv*goal.asMat();
		printMat(j);
		vector<double> r=states;
		for(unsigned int i=0;i<r.size();i++){
			r[i]+=dState.at<double>(i,1);
		}
		cout << dState<<endl;
		return r;
	}
	/*  B O D Y  */
	Pose backwardDownTo(int N,vector<double> &states){
		Pose p;
		for(int i=joints.size()-1;i>=N;i--){
			ScrewAxis &sa=joints[i];
			double state=states[i];
			p=p+sa*state;
		}
		return p;
	}
	Pose forwardBody(vector<double> &states){
		Pose p=forwardUpTo(joints.size(),states);
		p=home+p;
		return p;
	}
	void buildJacobianBodyColumns(vector<double> &states){
		for(unsigned int i=0;i<states.size();i++){
			Pose p=-backwardDownTo(i,states);
			Jb[i]=p.Ad(joints[i]);
		}
	}
	Mat jacobianBody(vector<double> &states){
		buildJacobianBodyColumns(states);
		Mat m=Jb[0].asMat();
		for(unsigned int i=0;i<J.size();i++){
			Mat mj=Jb[i].asMat();
			hconcat(m,mj,m);
		}
		return m;
	}
	vector<double> ikBodyStep(vector<double> &states,Twist goal){
		Mat j=jacobianBody(states);
		Mat Jpinv;
		invert(j, Jpinv, DECOMP_SVD);
		Mat dState=Jpinv*goal.asMat();
		printMat(j);
		vector<double> r=states;
		for(unsigned int i=0;i<r.size();i++){
			r[i]+=dState.at<double>(i,1);
		}
		cout << dState<<endl;
		return r;
	}
};
