/*
 * pose.h
 *
 *  Created on: 23 Nov 2020
 *      Author: Francisco Dominguez
 */

#ifndef POSE_H_
#define POSE_H_
#include <opencv2/viz.hpp>
#include "twist.h"

using namespace std;
using namespace cv;

const static Mat I=Mat::eye(3,3,CV_64F);
inline Mat hat(Vec3d w){
	double &wx=w[0];double &wy=w[1];double &wz=w[2];
	Mat W=(Mat_<double>(3,3)<<
			  0,-wz, wy,
			 wz,  0,-wx,
			-wy, wx,  0);
	return W;
}
inline Vec3d vee(Mat W){
	return Vec3d(W.at<double>(2,1),W.at<double>(0,2),W.at<double>(1,0));
}
inline double trace(Mat R){
	return R.at<double>(0,0)+R.at<double>(1,1)+R.at<double>(2,2);
}
inline Mat exp(Vec3d w){
	Mat W=hat(w);
	Mat W2=W*W;
	double theta=norm(w);
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
	return I+W*CA+W2*CB;
}
inline Vec3d log(Mat R){
	double tr=trace(R);
	if(tr>=2.999){//R==I
		return Vec3d(0,0,0);
	}
	double theta;
	//TODO
	if(tr<1+0.001 && tr>1-0.001){
		theta=3.14159265358979;
	}
	double cw=(tr-1)/2;
	theta=acos(cw);
	double sw=sin(theta);
	Mat W=(R-R.t())/(2*sw);
	Mat mw=W*theta;
	Vec3d w=vee(mw);
	return w;
}

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
	Pose(Vec3d w,Vec3d v)/*:Pose(Twist(w,v).exp())*/{
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
		R=I+W*CA+W2*CB;
		Mat V=I+W*CB+W2*CC;
		t=V*Mat(v);
	}
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
	Pose(Twist V):Pose(V.getW(),V.getV()){}
	inline Mat getR()const{return R;}
	inline Mat getT()const{return t;}
	inline Mat asMat(){
		Mat T;
		hconcat(R,t,T);
		Mat m=(Mat_<double>(1,4)<<0,0,0,1);
		T.push_back(m);
		return T;
	}
	inline void setR(Mat R){this->R=R;}
	inline void setT(Mat T){t=T;}
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
	// I am thinking about using + for composition and * for action,...
	// I could use - binary operator for right difference and - unary operator for inverse
	// then this a+-c is left difference and a-c=-c+a is right difference
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
	inline Twist adjoint(Vec3d w,Vec3d v){
		Mat mw(w);
		Mat mv(v);
		Mat rw=R*mw;
		Mat rr=hat(t)*rw+R*mv;
		return Twist(rw,rr);
	}
	inline Twist adjoint(Twist V){return adjoint(V.getW(),V.getV());}
	inline Mat adjointMat(){
		Mat adj;//=Mat::zeros(6,6,CV_64F);
		Mat pr=hat(t)*R;
		/*
		Mat m00(adj,Rect(0,0,3,3));
		Mat m10(adj,Rect(0,3,3,3));
		Mat m11(adj,Rect(3,3,3,3));
		for(int i=0;i<3;i++){
			for(int j=0;i<3;i++){
				m00.at<double>(i,j)=R .at<double>(i,j);
				m10.at<double>(i,j)=pr.at<double>(i,j);
				m11.at<double>(i,j)=R .at<double>(i,j);
			}
		}*/
		cout <<"pr="<<pr<<endl;
		Mat r1;
		hconcat(R,Mat::zeros(3,3,CV_64F),adj);
		hconcat(pr,R,r1);
		adj.push_back(r1);
		return adj;
	}
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
	os<<p.getR()<<p.getT();
	return os;
}
inline Pose exp(Twist tw){return Pose::exp(tw);}
inline Twist log(Pose p){return p.log();}
typedef Pose Joint;
#endif /* POSE_H_ */
