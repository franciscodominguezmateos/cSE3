/*
 * pose.h
 *
 *  Created on: 23 Nov 2020
 *      Author: Francisco Dominguez
 */

#ifndef POSE_H_
#define POSE_H_
#include "poe_util.h"
#include "twist.h"

class Pose{
	Mat R;// (3x3)
	Mat t;// (3x1)
public:
	Pose():R(I),t(Mat::zeros(3,1,CV_64F)){}
	Pose(double theta_deg,double x,double y){
		const double DEG_TO_RAD=3.14159265/180;
		double theta=theta_deg*DEG_TO_RAD;
		double sw=sin(theta);
		double cw=cos(theta);
		R=(Mat_<double>(3,3)<<cw,-sw,0,
				sw,cw,0,
				0,0,1);
		t=(Mat_<double>(3,1)<<x,y,0);
	}
	Pose(Mat R,Mat t):R(R),t(t){}
	Pose(Vec3d w,Vec3d v):Pose(Twist(w,v)){}
	Pose(Twist tw):Pose(exp(tw)){}
	static Pose exp(Twist tw){
		Vec3d w=tw.getW();
		Vec3d v=tw.getV();
		// w=u*theta;
		double theta=norm(w);
		//Vec3d u=w/theta;
		//Mat U=hat(u);
		//Mat U2=U*U;
		Mat W=hat(w);
		Mat W2=W*W;
		double theta2=theta*theta;
		double CA,CB,CC;
		if(theta<0.001){
			CA=1-theta2/6*(1-theta2/20*(1-theta2/42));
			CB=1/2*(1-theta2/12*(1-theta2/30*(1-theta2/56)));
			CC=1/6*(1-theta2/20*(1-theta2/42*(1-theta2/72)));
		}
		else{
			double theta3=theta2*theta;
			double sw=sin(theta);
			double cw=cos(theta);
			CA=sw/theta;
			CB=(1-cw)/theta2;
			CC=(theta-sw)/theta3;
		}
		Mat R=I+W*CA+W2*CB;
		Mat V=I+W*CB+W2*CC;
		Mat t=V*Mat(v);
		return Pose(R,t);
	}
	//Rotation
	inline Mat  getR() const{return R;}
	inline void setR(Mat R){this->R=R;}
	//Translation
	inline Mat  getT() const{return t;}
	inline void setT(Mat T){t=T;}
	inline Mat asMat()const{
		Mat T;
		hconcat(R,t,T);
		Mat m=(Mat_<double>(1,4)<<0,0,0,1);
		T.push_back(m);
		return T;
	}
	// Return transformation
	inline Mat T(){return asMat();}
	inline Twist log(){
		double tr=trace(R);
		if(tr<0.001){
			double l=norm(t);
			Mat vn=t/l;
			Vec3d v(vn.at<double>(0,0),vn.at<double>(1,0),vn.at<double>(2,0));
			return Twist(Vec3d(),v);
		}
		else{
			Vec3d w=::log(R);
			double theta=norm(w);
			Mat W=hat(w);
			Mat W2=W*W;
			double theta2=theta*theta;
			double CA,CB;
			if(theta<0.001){
				CA=1-theta2/6*(1-theta2/20*(1-theta2/42));
				CB=1/2*(1-theta2/12*(1-theta2/30*(1-theta2/56)));
			}
			else{
				double sw=sin(theta);
				double cw=cos(theta);
				CA=sw/theta;
				CB=(1-cw)/theta2;
			}
			Mat Vi=I-W*1/2+W2*(1/theta2)*(1-CA/(2*CB));
			Mat v=Vi*t;
			return Twist(w,v);
		}
	}
	// I am using + for composition and * for action.
	// Use - binary operator for right difference and - unary operator for inverse
	// then this a+-c is left difference and a-c=-c+a is right difference
	inline Pose inverse()  {return Pose(R.t(),-R.t()*t);}
	inline Pose i()        {return Pose(R.t(),-R.t()*t);}
	inline Pose operator-(){return Pose(R.t(),-R.t()*t);}
	inline Twist operator-(Pose q){
		Pose &p=*this;
		// p=q*d right operator
		Pose d=-q+p;
		Twist td=d.log();
		return td;
	}
	inline Pose operator+(Twist tw){
		Pose &self=*this;
		Pose p=self+exp(tw);
		return p;
	}
	inline Pose operator+(Pose p){
		return Pose(R*p.R,R*p.t+t);
	}
	inline Vec3d operator*(Vec3d p){
		Mat r=R*Mat(p)+t;
		return Vec3d(r);
	}
	//apply Pose action in all columns of m
	inline Mat operator*(Mat m){
		Mat r=R*m;
		//r=r+t;
		for(int j=0;j<r.cols;j++){
			r.col(j)+=t;
		}
		return r;
	}
	//Adjoint function Pose.Ad(Twist)
	inline Twist Ad(Vec3d w,Vec3d v){
		Mat mw(w);
		Mat mv(v);
		Mat rw=R*mw;
		Mat rr=hat(t)*rw+R*mv;
		return Twist(rw,rr);
	}
	inline Twist Ad(Twist V){return Ad(V.getW(),V.getV());}
	inline Twist AdT(Vec3d w,Vec3d v){
		Mat mw(w);
		Mat mv(v);
		Mat rr=hat(t)*R;
		Mat rw=R.t()*mw+rr.t()*mv;
		Mat rv=R.t()*mv;
		return Twist(rw,rv);
	}
	inline Twist AdT(Twist V){return AdT(V.getW(),V.getV());}
	//Adjoint matrix
	inline Mat AdMat(){
		Mat pr=hat(t)*R;
		Mat adj=stack4x4(R,Z,pr,R);
		return adj;
	}
	inline Mat AdTMat(){return AdMat().t();}
	//Jacobian of pose action at 0
	inline Mat J0(Vec3d py){
		Mat yx=-hat(py);
		Mat J;
		hconcat(yx,I,J);
		return J;
	}
	//Interpolation
	// this interpolated by right-> dm in on the right and update is r=T0*Tt
	// With matrices: m1=m0*dm ->dp=m0.inv()*m1;
	// With poses: p1=p0+dp  -> dp=p1-p0; //by right
	Pose interpolate(Pose &p1,double t){
		Pose &p0=*this;
		Pose dp=p1-p0;
		Twist dtw=dp.log()*t;
		Pose pr=p0+dtw;
		return pr;
	}
	//Integration
	inline static Pose integrate(Twist tw,double t){
		Pose pr=exp(tw*t);
		return pr;
	}
};
ostream& operator<<(ostream& os,const Pose& p){
	Mat m(p.asMat());
	printMat(m);
	return os;
}
inline Pose exp(Twist tw){return Pose::exp(tw);}
inline Twist log(Pose p){return p.log();}
typedef Pose Joint;
#endif /* POSE_H_ */
