/*
 * twist.h
 *
 *  Created on: 23 Nov 2020
 *      Author: francisco Dominguez
 */
#pragma once
#include "poe_util.h"

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
	//Adjoint function Pose.Ad(Twist)
	inline Twist ad(Vec3d w,Vec3d v){
		Mat mw(w);
		Mat mv(v);
		Mat W=hat(w);
		Mat V=hat(v);
		Mat rw=W*mw;
		Mat rv=V*mw+W*mv;
		return Twist(rw,rv);
	}
	inline Twist ad(Twist V){return ad(V.getW(),V.getV());}
	inline Twist adT(Vec3d w,Vec3d v){
		Mat mw(w);
		Mat mv(v);
		Mat W=hat(w);
		Mat V=hat(v);
		Mat rw=W.t()*mw+V.t()*mv;
		Mat rv=W.t()*mv;
		return Twist(rw,rv);
	}
	inline Twist adT(Twist V){return adT(V.getW(),V.getV());}
	//Adjoint matrix
	inline Mat adMat(){
		Mat W=hat(w);
		Mat V=hat(v);
		Mat ad=stack4x4(W,Z,V,W);
		return ad;
	}
	inline Mat adTMat(){return adMat().t();}
	//Jacobian of pose action at 0
	inline Mat J0(Vec3d py){
		Mat yx=-hat(py);
		Mat J;
		hconcat(yx,I,J);
		return J;
	}

	friend ostream& operator<<(ostream& os,const Twist& tw);
};
ostream& operator<<(ostream& os,const Twist& tw){
	os << "w="<<tw.w<<"v="<<tw.v;
	return os;
}
typedef Twist ScrewAxis;
typedef Twist Wrench;

