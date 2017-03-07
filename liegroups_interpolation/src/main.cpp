#include <iostream>
#include <liegroups/se3_io.hpp>
#include <liegroups/se3.hpp>

#include <opencv2/viz.hpp>

#include <cassert>
#include "rigid_transformation.h"

using namespace std;
using namespace liegroups;
using namespace cv;

typedef double S;
typedef SE3<S> SE3Type;
typedef S Tangent[6];

SE3Type interpolation(SE3Type &a,SE3Type &b, S t){
	SE3Type d;
	multiply_a_binv(d,b,a);
	//cout << "d" <<d <<endl;
	Tangent td;
	log(td,d);
	SE3Type de;
	exp(de,td);
	//cout << "de"<< de <<endl;
	Tangent dt;
	for(int i=0;i<6;i++) dt[i]=td[i]*t;
	SE3Type dtg;
	exp(dtg,dt);
	SE3Type ret=dtg*a;
	return ret;
}
Mat se3toMat(SE3Type se3){
	//cout << "se3"<<endl;
	//cout << se3 <<endl;
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
	//cout << "Mat"<<endl;
	//cout << m <<endl;
	return m;
}
void setPose(viz::Viz3d w,string name,SE3Type &pose){
	cv::Affine3d a(se3toMat(pose));
	w.setWidgetPose(name,a);
}
int main()
{
	viz::Viz3d myWindow("Estimation Coordinate Frame");
	//myWindow.setBackgroundColor(); // black by default
	//myWindow.registerKeyboardCallback(&keyboard_callback);
	myWindow.showWidget("Coordinate Widget", viz::WCoordinateSystem());
    viz::WCube cube_widget(Point3f(1.0,1.0,0.0), Point3f(0.0,0.0,-1.0), true, viz::Color::blue());
    cube_widget.setRenderingProperty(viz::LINE_WIDTH, 1.0);
    myWindow.showWidget("cube",cube_widget);
    SE3<double> p0,p1,pi;
    cout << "exp(log(x))==x"<<endl;
    double ap1[6]={1.3112,0.8507,1.5186,0.8851,0.2362,-0.0898};
    double ap0[6]={1.08027,2.1904,-0.270935,-2.24856,-0.600055,0.228133};
    for(int i=0;i<6;i++) cout << ap0[i]; cout << endl;
    exp(p0,ap0);
    cout << p0 << endl;
    log(ap0,p0);
    for(int i=0;i<6;i++) cout << ap0[i]; cout << endl;
    cout << "p1"<<endl;
    exp(p1,ap1);
    cout << p1 << endl;
    pi=interpolation(p0,p1,1);
    cout << "pi"<<endl;
    cout << pi << endl;
    viz::WGrid g;
    //viz::WText3D txt("hola Viz",Point3f(0,0,3));
	viz:: WCameraPosition a,b;
	viz:: WCameraPosition wcp0,wcp1,wcpi;
    /// We can get the transformation matrix from camera coordinate system to global using
    /// - makeTransformToGlobal. We need the axes of the camera
    cv::Affine3d transform = viz::makeTransformToGlobal(Vec3f(0.0f,-1.0f,0.0f), Vec3f(-1.0f,0.0f,0.0f), Vec3f(0.0f,0.0f,-1.0f),Vec3f(2,1,0));
    //cv::Affine3f pose(d.matrix3x4());
    //a.setPose(transform);
    //myWindow.showWidget("txt",txt);
    vector<Point3f> vp1,vp2;
    vp1.push_back(Point3f(0,1,0));
    vp1.push_back(Point3f(1,1,0));
    vp1.push_back(Point3f(0,1,-1));
    vp1.push_back(Point3f(0.4,1,0.2));
    vp1.push_back(Point3f(1.1,-1,0));
    vp1.push_back(Point3f(-1,1,-1));
    cout << "Before of transform"<<endl;
    cv::perspectiveTransform(vp1,vp2,transform.matrix);
    cout << "Before of reshape"<<endl;
	Mat mp1=Mat(vp1).reshape(1);
	Mat mp2=Mat(vp2).reshape(1);
    cout << mp2 << endl;
    cout << mp1 << endl;
	Mat Rll,tll;
    rigidTransformation(mp1,mp2,Rll,tll);
    cout << mp2 << endl;
    cout << mp1 << endl;
    cout << transform.matrix << endl;
    cout << Rll << endl;
    cout << tll << endl;
    myWindow.showWidget("g",g);
	myWindow.showWidget("p0",wcp0);
	myWindow.showWidget("p1",wcp1);
	myWindow.showWidget("pi",wcpi);
    setPose(myWindow,"p0",p0);
    setPose(myWindow,"p1",p1);
	//myWindow.setWidgetPose("a",transform);
	myWindow.setWidgetPose("cube",transform);

    while(!myWindow.wasStopped())
    {
    	//cout << "Izquierda"<<endl;
    	for(double t=0;t<1;t+=0.005){
    	    pi=interpolation(p0,p1,t);
    	    setPose(myWindow,"pi",pi);
            myWindow.spinOnce(10, true);
    	}
    	//cout << "Derecha"<<endl;
    	for(double t=1;t>0;t-=0.005){
    	    pi=interpolation(p0,p1,t);
    	    setPose(myWindow,"pi",pi);
            myWindow.spinOnce(10, true);
    	}
    }
    myWindow.spin();
	return 0;
}
