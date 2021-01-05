/*
 * pose_solver.h
 *
 *  Created on: 4 Dec 2020
 *      Author: Francisco Dominguez
 *  TODO: with Jacobian at 0 that means alwais the same jacobian but change origin points
 */
#pragma once
#include "pose_problem.h"
class PoseSolver{
	double alpha;//learning rate
	double minE; //min error to stop=>precision
	int MAXI;    //max iterations
	int i;       //Iterations
	//Needed in Gauss-Newton
	Mat H; //J.t*J
	Mat g; //J.t*r
	double els;
	bool gaussNewton;
public:
	PoseSolver():alpha(0.01),minE(0.1),MAXI(1000),gaussNewton(false){}
	double getEls(){return els;}
	void setHalfLearningRate(){alpha/=2.0;}
	void setMethodGaussNewton()    {gaussNewton=true ;alpha=0.99;}
	void setMethodGradientDescent(){gaussNewton=false;alpha=0.001;}
	//              G A U S S - N E W T O N
	inline void Hg(PoseProblem* p){
		int k=6;
		H=Mat::zeros(k, k, CV_64F);
		g=Mat::zeros(k, 1, CV_64F);
		for(int i=0;i<p->size();i++){
			Mat J=p->Jres(i);
			Mat JT=J.t();
			Mat re=p->res(i);
			H+=JT*J;
			g+=JT*re;
		}
	}
	inline Mat gradGaussNewton(PoseProblem* p){
		Hg(p);
		Mat Hi=H.inv();
		Mat mTheta=Hi*g;
		return mTheta;
	}
	//         G R A D I E N T   D E S C E N T
	inline Mat gradDescent(PoseProblem* p){return p->gradEls();}
	//         G E N E R I C   S O L V E R
	inline Mat getGrad(PoseProblem* p){
		if(gaussNewton){
			return gradGaussNewton(p);
		}
		return gradDescent(p);
	}
	inline double step(PoseProblem* p){
		Pose &T=p->getPose();
		Mat gradientEls=getGrad(p);
		Mat dm=-gradientEls*alpha;
		Twist deltaT=dm;
		T=T+deltaT;
		double els=p->Els();
		return els;
	}
	inline bool solve(PoseProblem* p){
		i=0;
    	els =step(p);
	    while(els>minE && i<MAXI){
	    	els =step(p);
	    	i++;
	    }
	    return i==MAXI+1;
	}
};


