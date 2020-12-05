/*
 * pose_problem_points.h
 *
 *  Created on: 4 Dec 2020
 *      Author: Francisco Dominguez
 */

#pragma once
#include "pose_problem.h"

class PoseProblemPoints:public PoseProblem{
public:
	//Target values
	Mat t;
	//Original values
	Mat p;
	//Transformed values
	Mat q;
	//Transformation that make t=T*p or t-T*p=0
	Pose T;
	// Residual or r=t-T*p
	Mat r;
	// Error r2=r.mul(r) function to minimize
	Mat r2;
	// Size is number of columns in all Mat the same
	int size(){return p.cols;}
	// Compute all expressions above
	void build(){
		q=T*p;
		r=t-q;
		r2=r.mul(r);
	}
	//Error less square Els=sum(r2).val(0)
	double Els(){build();return sum(r2)[0]/double(r2.cols);}
	inline Mat qCol(int i){return Mat(q,Range(0,3),Range(i,i+1));}
	Mat rCol(int i){return Mat(r,Range(0,3),Range(i,i+1));}
	Mat Jres(int i){return -T.Jaction(Vec3d(qCol(i)));}
	//Gradient of error with respecto to parameter of T
	Mat gradEls(){
		//Jacobian of Els
		Mat jEls=Mat::zeros(Size(1,6),CV_64F);
		for(int i=0;i<size();i++){
			Mat Jr=Jres(i);
			Mat gradR=Jr.t()*rCol(i);
			jEls+=gradR;
		}
		return jEls/double(size());
	}
	Pose &getPose(){return T;}
};
ostream &operator<<(ostream &os,PoseProblemPoints &ppp){
	cout << "Els="<<ppp.Els()<<endl;
	cout << "p ="<<ppp.p <<endl;
	cout << "t ="<<ppp.t <<endl;
	cout << "q ="<<ppp.q <<endl;
	cout << "r ="<<ppp.r <<endl;
	cout << "r2="<<ppp.r2<<endl;
	return os;
}

