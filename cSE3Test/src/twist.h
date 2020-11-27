/*
 * twist.h
 *
 *  Created on: 23 Nov 2020
 *      Author: francisco
 */

#ifndef SRC_TWIST_H_
#define SRC_TWIST_H_
#include <opencv2/viz.hpp>
//#include "pose.h"

using namespace std;
using namespace cv;

class Twist{
	Vec3d w,v;
public:
	Twist():Twist(0,0,0,0,0,0){}
	Twist(double wx,double wy,double wz,
		  double vx,double vy,double vz):w(wx,wy,wz),v(vx,vy,vz){}
	Twist(Vec3d w,Vec3d v):w(w),v(v){}
	inline Vec3d getW(){return w;}
	inline Vec3d getV(){return v;}
	inline Vec3d getM(){return w;}
	inline Vec3d getF(){return v;}
	inline Mat asMat(){
		Mat mw(w);
		Mat mv(v);
		mw.push_back(mv);
		return mw;
	}
	inline Twist operator+(Twist  tw){return Twist(w+tw.w,v+tw.v);}
	inline Twist operator*(double dt){return Twist(w*dt,v*dt);}
	//inline Pose operator*(Pose p){
	//	return Pose::exp(*this)*p;
	//}
	//Pose exp(){/*return Pose::exp(*this);*/}
	Twist getScrewAxis(){
		double n=norm(w);
		if(n>0.0001){
			return Twist(w/n,v/n);
		}
		else{
			n=norm(v);
			return Twist(Vec3d(),v/n);
		}
	}
	double getTheta(){
		double n=norm(w);
		if(n<0.0001){
			n=norm(v);
		}
		return n;
	}
	friend ostream& operator<<(ostream& os,const Twist& tw);
};
ostream& operator<<(ostream& os,const Twist& tw){
	os << "w="<<tw.w<<"v="<<tw.v;
	return os;
}
typedef Twist ScrewAxis;
typedef Twist Wrench;
#endif /* SRC_TWIST_H_ */
