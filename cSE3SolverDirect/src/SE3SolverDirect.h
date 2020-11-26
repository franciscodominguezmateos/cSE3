/*
 * SE3SolverDirect.h
 *
 *  Created on: Oct 14, 2016
 *      Author: Francisco Dominguez
 */

#ifndef SE3SOLVER_DIRECT_H_
#define SE3SOLVER_DIRECT_H_
#include <iostream>
#include <vector>
//#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "depthImage.h"
#include "SE3Group.h"

using namespace std;

template <typename S>
class SE3SolverDirect {
	S alpha;//learning rate
	S minE; //min error to stop=>precision
	unsigned int MAXI;//max iterations
	Mat H; //J.t*J
	Mat g; //J.t*r
	DepthImage I2,I1;
	Mat I2p;
	vector<Point3f> P;//3Dpts in I1
	vector<Point2f> U;//2Dpts valid points in I1
	vector<Point3f> Q;//3Dpts transformed from I1
	vector<Point2f> W;//2Dpts projected from Q
public:
	typedef SE3Group<S> SE3Type;
	typedef se3Algebra<S> Tangent;
	SE3SolverDirect(){
		alpha=0.01;
		minE=0.001;
		MAXI=10000;
	}
	/*
	Point2f i2p(int i){
		int x=i % I2.getImg().cols;
		int y=i / I2.getImg().rows;
		return Point2f(x,y);
	}*/
	vector<Point2f> &warp(Tangent theta){
		transformPoints(Q,P,theta);
		W.clear();
		for(Point3f p3D:Q){
			W.push_back(I1.project(p3D));
		}
		return W;
	}
	Mat wrapI2(Tangent theta){
		I2p=Mat::zeros(I2.getImg().size(),I2.getImg().type());
		vector<Point2f> &wp=warp(theta);
		for(int i=0;i<U.size();i++){
			if(I2.is2DPointInImage(W[i])){
				Vec3b c=I2.getColor(W[i]);
				I2p.at<Vec3b>(U[i])=c;
			}//else cout <<"id2DPointInImage wrapI2"<<endl;
		}
		/*
		for(Point2f p2f:wp){
			if(I2.is2DPointInImage(p2f)){
				Vec3b c=I2.getColor(p2f);
				I2p.at<Vec3b>(p2f.y,p2f.x)=c;
			}else cout <<"id2DPointInImage wrapI2"<<endl;
		}*/
		return I2p;
	}
	inline void setAlpha(S a){alpha=a;}
	//virtual ~SE3Solver();

	//Transformation function where the parameters (theta) live
	inline Point3f T(Point3f &p,Tangent &theta){
		return theta*p;
	}
	inline void transformPoints(vector<Point3f> &Q,vector<Point3f> P,Tangent theta){
		Q.clear();
		for(unsigned int i=0;i<P.size();i++){
			Q.push_back(T(P[i],theta));
			//cout << Q[i]<<endl;
		}
	}
	//Jacobian with respect to theta at 0, of transformation function that return a 3x6 matrix
	inline Mat JT(int i){
		Point3f pt=Q[i];
		float &x=pt.x;
		float &y=pt.y;
		float &z=pt.z;
		Mat jT=Mat_<S>(3,6);
		jT.at<S>(0,0)=  1;jT.at<S>(0,1)=  0;jT.at<S>(0,2)=  0;jT.at<S>(0,3)=  0;jT.at<S>(0,4)=  z;jT.at<S>(0,5)= -y;
		jT.at<S>(1,0)=  0;jT.at<S>(1,1)=  1;jT.at<S>(1,2)=  0;jT.at<S>(1,3)= -z;jT.at<S>(1,4)=  0;jT.at<S>(1,5)=  x;
		jT.at<S>(2,0)=  0;jT.at<S>(2,1)=  0;jT.at<S>(2,2)=  1;jT.at<S>(2,3)=  y;jT.at<S>(2,4)= -x;jT.at<S>(2,5)=  0;
		return jT;
	}
	// Residual function
	inline S r(int i,Tangent &theta){
		float dg= I2.getGray(W[i])-I1.getGray(U[i]);
		return dg;
	}
	inline Mat rImg(Tangent &theta){
		Mat rm=Mat::zeros(I2.getImg().size(),CV_32F);
		for(int i=0;i<U.size();i++){
			if(I2.is2DPointInImage(W[i])){
				rm.at<float>(U[i])=r(i,theta);
			}//else cout <<"id2DPointInImage rImg"<<endl;
		}
		return rm;
	}
	inline Mat Jproj(int i){
		Point3f pt=Q[i];
		float &x=pt.x;
		float &y=pt.y;
		float &z=pt.z;
		float fx=I2.getFx();
		float fy=I2.getFy();
		Mat jp=(Mat_<S>(2,3)<< fx/z,   0,-fx*x/(z*z),
				               0   ,fy/z,-fy*y/(z*z));
		return jp;
	}
	// Jacobian of residual with respect to theta at 0
	inline Mat Jr(int i,Tangent &theta){
		Mat Ji=(Mat_<S>(1,2)<<I2.getGradXImg().at<float>(W[i]),
				              I2.getGradYImg().at<float>(W[i]));
		return Ji*Jproj(i)*JT(i);
	}
	// Residual squared is the loss function for point i
	inline S R(int i,Tangent &theta){
		S pr=r(i,theta);
		return pr*pr;
	}
	inline Mat RImg(Tangent &theta){
		Mat rm=Mat::zeros(I2.getImg().size(),CV_32F);
		for(int i=0;i<U.size();i++){
			if(I2.is2DPointInImage(W[i])){
				rm.at<float>(U[i])=sqrt(R(i,theta));
			}//else cout <<"id2DPointInImage RImg"<<endl;
		}
		return rm;
	}
	// Gradient of loss function with respect to theta in order to do gradient descent
	inline Mat GR(int i,Tangent &theta){
		Mat jrT=Jr(i,theta).t();
		Mat mr=(Mat_<S>(1,1)<<r(i,theta));
		return jrT*mr;
	}
	// Less squares error
	inline S Els(Tangent theta){
		S e=0;
		int c=0;//number of good 2D points
		for(unsigned int i=0;i<U.size();i++){
			if(I2.is2DPointInImage(W[i])){
				e+=R(i,theta);
				c++;
			}//else cout <<"id2DPointInImage Els"<<endl;
		}
		return e/c;
	}
	// Less squares gradient
	inline Mat GEls(Tangent &theta){
		Mat jels=Mat::zeros(Size(1,6),CV_64F);
		int c=0;//number of good 2D points
		for(unsigned int i=0;i<U.size();i++){
			if(I2.is2DPointInImage(W[i])){
				jels+=GR(i,theta);
				c++;
			}//else cout <<"id2DPointInImage GEls"<<endl;
		}
		return jels/c;
	}
	inline void Hg(vector<Point3f> &Q,vector<Point3f> &P,Tangent &theta){
		int k=6;
		H=Mat::zeros(k, k, CV_64F);
		g=Mat::zeros(k, 1, CV_64F);
		for(int i=0;i<Q.size();i++){
			Mat J=Jr(i,Q,P,theta);
			Mat JT=J.t();
			Mat re=Mat_<S>(3,1);
			Point3f p=r(i,Q,P,theta);
			re.at<S>(0,0)=p.x;re.at<S>(1,0)=p.y;re.at<S>(2,0)=p.z;
			H+=JT*J;
			g+=JT*re;
		}
	}/*
	def Hg(P,theta):
	    k=theta.shape[0]
	    H=np.matrix(np.eye(k,k))
	    g=np.matrix(np.zeros((k,1)))
	    T=getT(theta)
	    Q=T*P
	    n=shape(P)[1]
	    for i in range(n):
	        J=JF(i,Q,P,theta)
	        JT=J.T
	        D=SDF(P[:,i])
	        H+=JT*J
	        g+=JT*D;
	    return H,g*/
	/*
	inline void stepGradientDescent(){
        transformPoints(vp1,vp1,dThetaAlpha);
        Mat gradientEls=GEls(vp3,vp1,dThetaAlpha);
        Tangent dTheta(gradientEls);
        dThetaAlpha=dTheta*-alpha;
        theta+=dThetaAlpha;
	}*/
	//Forward Gradient Descent
	inline bool solveGradientDescent(Tangent &theta){
		S els;
		Mat gels;
		Tangent dTheta,dThetaAlpha;
		theta=Tangent();
		int i=0;
    	els =Els(dThetaAlpha);
	    while(sqrt(els)>minE && i++<MAXI){
	        gels=GEls(theta);
	        Tangent dTheta(gels);
	        dThetaAlpha=dTheta*-alpha;
	        theta+=dThetaAlpha;
	    	els =Els (dThetaAlpha);
	    }
	    cout <<"GD"<<i<<endl;
	    return i==MAXI+1;
	}
	inline Mat stepGaussNewton(vector<Point3f> &Q,vector<Point3f> &P,Tangent &dThetaAlpha){
		Hg(Q,P,dThetaAlpha);
		Mat Hi=H.inv();
		Mat mTheta=Hi*g;
		return mTheta;
	}
	inline bool solveGaussNewton(Tangent &theta){
		vector<Point3f> vp(P);
		S els;
		Mat gels;
		Tangent dTheta,dThetaAlpha;
		theta=Tangent();
		int i=0;
    	els =Els (dThetaAlpha);
	    while(sqrt(els)>minE && i++<MAXI){
	        gels=stepGaussNewton(Q,vp,dThetaAlpha);
	        Tangent dTheta(gels);
	        dThetaAlpha=dTheta*-alpha;
	        theta+=dThetaAlpha;
	        transformPoints(vp,vp,dThetaAlpha);
	    	els =Els (Q,vp,dThetaAlpha);
	    }
	    cout <<"GN"<<i<<endl;
	    return i==MAXI+1;
		//return false;
	}

	const DepthImage& getI1() const {
		return I1;
	}

	void setI1(const DepthImage& i1) {
		I1 = i1;
		U=I1.getPoints2D();
		P=I1.getPoints3D();
	}
	DepthImage& getI2(){
		return I2;
	}
	const DepthImage& getI2p() const {
		return I2p;
	}

	void setI2(const DepthImage& i2) {
		I2 = i2;
	}
};

#endif /* SE3SOLVER_H_ */
