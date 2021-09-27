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
		Jb.resize(tw.size());
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
		int i=joints.size()-1;
		while(i>N){
			ScrewAxis &sa=joints[i];
			double state=states[i];
			p=sa*state+p;
			i--;
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
		for(unsigned int i=1;i<J.size();i++){
			Mat mj=Jb[i].asMat();
			hconcat(m,mj,m);
		}
		return m;
	}
	vector<double> ikBodyStep(vector<double> &states,Pose goal){
		Pose&Tsd=goal;
		Pose Tsb=forwardBody(states);
		Twist tbd=Tsd-Tsb;
		cout << "tbd="<<tbd<<endl;
		Mat j=jacobianBody(states);
		Mat Jpinv;
		invert(j, Jpinv, DECOMP_SVD);
		Mat dState=Jpinv*tbd.asMat();
		cout <<"J="<<endl;
		printMat(j);
		cout << "Ji="<<endl;
		printMat(Jpinv);
		vector<double> r=states;
		for(unsigned int i=0;i<r.size();i++){
			r[i]+=dState.at<double>(i,1);
		}
		cout << "dState="<< dState<<endl;
		return r;
	}
};
class OpenChainBodyBody{
	vector<ScrewAxis> joints;
	//jacobian columns
	vector<Twist> J;  //Space/stationary/world/global Jacobian
	vector<Twist> Jb; //Body Jacobian
	Pose home; //"zero" position of the robot
public:
	OpenChainBodyBody(vector<ScrewAxis> tw,Pose h):joints(tw),home(h){
		J.resize(tw.size());
		Jb.resize(tw.size());
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
		int i=joints.size()-1;
		while(i>N){
			ScrewAxis &sa=joints[i];
			double state=states[i];
			p=p+sa*state;
			i--;
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
		for(unsigned int i=1;i<J.size();i++){
			Mat mj=Jb[i].asMat();
			hconcat(m,mj,m);
		}
		return m;
	}
	vector<double> ikBodyStep(vector<double> &states,Pose goal){
		Pose&Tsd=goal;
		Pose Tsb=forwardBody(states);
		Twist tbd=Tsd-Tsb;
		cout << "tbd="<<tbd<<endl;
		Mat j=jacobianBody(states);
		Mat Jpinv;
		invert(j, Jpinv, DECOMP_SVD);
		Mat dState=Jpinv*tbd.asMat();
		cout <<"J="<<endl;
		printMat(j);
		cout << "Ji="<<endl;
		printMat(Jpinv);
		vector<double> r=states;
		for(unsigned int i=0;i<r.size();i++){
			r[i]+=dState.at<double>(i,1);
		}
		cout << "dState="<< dState<<endl;
		return r;
	}
};
class OpenChainBodySpace{
	vector<ScrewAxis> joints;
	//jacobian columns
	vector<Pose> forwardUp;
	vector<Twist> J;  //Space/stationary/world/global Jacobian
	vector<Twist> Jb; //Body Jacobian
	Pose home; //"zero" position of the robot
public:
	OpenChainBodySpace(vector<ScrewAxis> tw,Pose h):joints(tw),home(h){
		J        .resize(tw.size());
		Jb       .resize(tw.size());
		forwardUp.resize(tw.size());
	}
	/*  S P A C E  */
	Pose forwardUpAll(vector<double> &states){
		const int& N=states.size();
		Pose p;
		for(int i=0;i<N;i++){
			ScrewAxis &sa=joints[i];
			double state=states[i];
			p=p+sa*state;
			forwardUp[i]=p;
		}
		return p;
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
		//Pose p=forwardUpTo(joints.size(),states);
		Pose p=forwardUpAll(states);
		p=p+home;
		return p;
	}
	void buildJacobianColumns(vector<double> &states){
		Pose p=forwardUpAll(states);
		for(unsigned int i=0;i<states.size();i++){
			//Pose p=forwardUpTo(i,states);
			Pose p=forwardUp[i];
			J[i]=p.Ad(joints[i]);
		}
	}
	Mat concatJacobianColumns(){
		Mat m=J[0].asMat();
		for(unsigned int i=1;i<J.size();i++){
			Mat mj=J[i].asMat();
			hconcat(m,mj,m);
		}
		return m;
	}
	Mat jacobian(vector<double> &states){
		buildJacobianColumns(states);
		Mat j=concatJacobianColumns();
		return j;
	}
	vector<double> ikStep(vector<double> states,Twist goal){
		Mat j=jacobian(states);
		Mat Jpinv;
		invert(j, Jpinv, DECOMP_SVD);
		Mat dState=Jpinv*goal.asMat();
		printMat(j);
		for(unsigned int i=0;i<states.size();i++){
			states[i]+=dState.at<double>(i,0);
		}
		cout << dState<<endl;
		return states;
	}
};
