/*
 * SE3Group.h
 *
 *  Created on: Oct 13, 2016
 *      Author: Francisco Dominguez
 *   This is a wraper to the lib liegroups of Ethan
 */

#ifndef SE3GROUP_H_
#define SE3GROUP_H_
#include <iostream>
//#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <liegroups/se3_io.hpp>
#include <liegroups/se3.hpp>

using namespace std;
using namespace liegroups;
using namespace cv;


template <typename S>
class se3Algebra;

//typedef se3Algebra<double> se3;

template <typename S>
class SE3Group {
    friend class se3Algebra<S>;
	SE3<S> se3;
public:
	SE3Group(){
		S t[6]={0,0,0,0,0,0};
		exp(se3,t);
	}
	SE3Group(Mat m){
		for(int i=0;i<3;i++){
			se3.R.R[3*i+0]=m.at<S>(i,0);
			se3.R.R[3*i+1]=m.at<S>(i,1);
			se3.R.R[3*i+2]=m.at<S>(i,2);
			se3.t[i]=m.at<S>(i,3);
		}
	}
	SE3Group(se3Algebra<S> &a){
		exp(se3,a.asPoiter());
	}
	SE3Group(S vx,S vy,S vz,S wx,S wy,S wz){
		S t[6]={vx,vy,vz,wx,wy,wz};
		exp(se3,t);
	}
	//virtual ~SE3Group();
	Mat asMat(){
		Mat m=Mat_<S>(4,4);
		for(int i=0;i<3;i++){
			m.at<S>(i,0)=se3.R.R[3*i+0];
			m.at<S>(i,1)=se3.R.R[3*i+1];
			m.at<S>(i,2)=se3.R.R[3*i+2];
			m.at<S>(i,3)=se3.t[i];
		}
		m.at<S>(3,0)=0;
		m.at<S>(3,1)=0;
		m.at<S>(3,2)=0;
		m.at<S>(3,3)=1;
		return m;
	}
	void invert(){
		liegroups::invert(se3);
	}
	SE3Group inverse(){
		SE3Group<S> se;
		se.se3=liegroups::inverse(se3);
		return se;
	}
	se3Algebra<S> log(){
		S t[6];
		liegroups::log(t,se3);
		return se3Algebra<S>(t[0],t[1],t[2],t[3],t[4],t[5]);
	}
	SE3Group<S> operator*(const SE3Group &s){
		//cout << "SEGroup::operator*"<<endl;
		SE3Group<S> se;
		se.se3=se3*s.se3;
		return se;
	}
	SE3Group<S> operator/(const SE3Group &s){
		SE3Group<S> se;
		multiply_a_binv(se.se3,se3,s.se3);
		return s;
	}
	SE3Group<S> operator*(const se3Algebra<S> &a){
		SE3Group<S> se;
		se.se3=se3*a.log().se3;
		return se;
	}
	friend SE3Group<S> operator*(const se3Algebra<S> &a,const SE3Group<S> s){
		SE3Group<S> se;
		se.se3=a.log().se3*s.se3;
		return se;
	}
	void operator*=(const SE3Group<S> &s){
		se3=se3*s.se3;
	}
	void operator/=(const SE3Group<S> &s){
		multiply_a_binv(se3,se3,s.se3);
	}
	void operator*=(se3Algebra<S> &aa){
		se3=se3*aa.log().se3;
	}
	Point3f operator*(Point3f &p){
		S v[3]={p.x,p.y,p.z};
		S r[3];
	    transform_point(r, se3, v);
		return Point3f(r[0],r[1],r[2]);
	}
	SE3Group<S> interpolation(SE3Group &s,S t){
		SE3<S> &a=se3;
		SE3<S> &b=s.se3;
		SE3<S> d;
		multiply_a_binv(d,b,a);
		S td[6];
		log(td,d);
		S dt[6];
		for(int i=0;i<6;i++)
			dt[i]=td[i]*t;
		SE3<S> dtg;
		exp(dtg,dt);
		SE3Group<S> ret;
		ret.se3=dtg*a;
		return ret;
	}
	se3Algebra<S> adjoint(se3Algebra<S> v){
		S r[6];
		adjoint_multiply(r,se3,v.t);
		return se3Algebra<S>(r);
	}
	se3Algebra<S> adjointT(se3Algebra<S> v){
		S r[6];
		adjoint_T_multiply(r,se3,v.t);
		return se3Algebra<S>(r);
	}
};
template <typename S>
class se3Algebra{
	S t[6];
public:
	friend class SE3Group<S>;
	se3Algebra(S s[6]){
		for(int i;i<6;i++) t[i]=s[i];
	}
	se3Algebra(){
		t[0]=t[1]=t[2]=t[3]=t[4]=t[5]=0;
	}
	se3Algebra(S vx,S vy,S vz,S wx,S wy,S wz){
		t[0]=vx;t[1]=vy;t[2]=vz;
		t[3]=wx;t[4]=wy;t[5]=wz;
	}
	se3Algebra(S &vx,S &vy,S &vz,S &wx,S &wy,S &wz){
		t[0]=vx;t[1]=vy;t[2]=vz;
		t[3]=wx;t[4]=wy;t[5]=wz;
	}
	se3Algebra(SE3Group<S> &se){
		log(t,se.se3);
	}
	se3Algebra(Mat m){
		if(m.rows==1 && m.cols==6){
			t[0]=m.at<S>(0,0);
			t[1]=m.at<S>(0,1);
			t[2]=m.at<S>(0,2);
			t[3]=m.at<S>(0,3);
			t[4]=m.at<S>(0,4);
			t[5]=m.at<S>(0,5);
		}
		else{
			t[0]=m.at<S>(0,0);
			t[1]=m.at<S>(1,0);
			t[2]=m.at<S>(2,0);
			t[3]=m.at<S>(3,0);
			t[4]=m.at<S>(4,0);
			t[5]=m.at<S>(5,0);
		}
	}
	SE3Group<S> exp(){
		return SE3Group<S>(t[0],t[1],t[2],t[3],t[4],t[5]);
	}
	se3Algebra<S> operator*(S s){
		return se3Algebra<S>(t[0]*s,t[1]*s,t[2]*s,t[3]*s,t[4]*s,t[5]*s);
	}
	friend se3Algebra<S> operator*(S s,se3Algebra<S> a){
		return a*s;
	}
	Point3f operator*(Point3f &p){
		return exp()*p;
	}
	void operator*=(S s){
		t[0]*=s;t[1]*=s;t[2]*=s;
		t[3]*=s;t[4]*=s;t[5]*=s;
	}
	se3Algebra<S> operator+(se3Algebra<S> a){
		return (exp()*a.exp()).log();
	}
	//must realize that actually the sum is this=a+this NOT this=this+a
	void operator+=(se3Algebra<S> a){
		SE3Group<S> seg=a.exp()*exp();
		SE3<S> se=seg.se3;
		log(t,se);
	}
	S* asPointer(){
		return t;
	}
	Mat asMat(){
		Mat m=Mat_<S>(6,1);
		for(int i=0;i<6;i++)
			m.at<S>(i,0)=t[i];
		return m;
	}
	//Jacobian with respect to theta at 0, of application function,  returns a 3x6 matrix
	inline Mat Jat0toP(Point3f &p){
		se3Algebra<S> &theta=*this;
		Point3f pt=theta*p;
		float &x=pt.x;
		float &y=pt.y;
		float &z=pt.z;
		Mat jT=Mat_<S>(3,6);
		jT.at<S>(0,0)=  1;jT.at<S>(0,1)=  0;jT.at<S>(0,2)=  0;jT.at<S>(0,3)=  0;jT.at<S>(0,4)=  z;jT.at<S>(0,5)= -y;
		jT.at<S>(1,0)=  0;jT.at<S>(1,1)=  1;jT.at<S>(1,2)=  0;jT.at<S>(1,3)= -z;jT.at<S>(1,4)=  0;jT.at<S>(1,5)=  x;
		jT.at<S>(2,0)=  0;jT.at<S>(2,1)=  0;jT.at<S>(2,2)=  1;jT.at<S>(2,3)=  y;jT.at<S>(2,4)= -x;jT.at<S>(2,5)=  0;
		return jT;
	}
};
#endif /* SE3GROUP_H_ */
