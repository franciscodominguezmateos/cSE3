/*
 * poe_util.h
 *
 *  Created on: 28 Nov 2020
 *      Author: Francisco Dominguez
 */
#pragma once
#include <iostream>
#include <cmath>
#include <opencv2/viz.hpp>

using namespace std;
using namespace cv;

const static double PI=3.14159265358979;
const static Mat I=Mat::eye(3,3,CV_64F);
const static Mat Z=Mat::zeros(3,3,CV_64F);

void printMat(Mat m){
	int i,j;
	for(i=0;i<m.rows;i++){
		for(j=0;j<m.cols;j++){
			double d=m.at<double>(i,j);
			if(d>=0) cout <<" ";
			cout << d;
			if(j!=m.cols-1 || i!=m.rows-1) cout<<",";
		}
		cout<<endl;
	}
}

Mat stack4x4(Mat m00,Mat m01,Mat m10, Mat m11){
	Mat r0,r1;
	hconcat(m00,m01,r0);
	hconcat(m10,m11,r1);
	r0.push_back(r1);
	return r0;
}
Mat diag4x4(Mat m00, Mat m11){
	Mat Z=Mat::zeros(Size(m00.rows,m00.cols),CV_64F);
	return stack4x4(m00,Z,Z,m11);
}
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
