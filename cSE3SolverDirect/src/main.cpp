#include <iostream>
//#include <liegroups/se3_io.hpp>
//#include <liegroups/se3.hpp>

#include <sstream>
//#include <opencv2/viz.hpp>

#include <cassert>
//#include "rigid_transformation.h"
#include "SE3Group.h"
//#include "SE3Solver.h"
#include "SE3SolverDirect.h"

using namespace std;
//using namespace liegroups;
using namespace cv;

typedef double S;
typedef SE3SolverDirect<S> Solver;
typedef typename SE3SolverDirect<S>::SE3Type SE3Type;
typedef typename SE3SolverDirect<S>::Tangent Tangent;
/*
void setPose(viz::Viz3d w,string name,SE3Type pose){
	cv::Affine3d a(pose.asMat());
	w.setWidgetPose(name,a);
}
void setPose(viz::Viz3d w,string name,Tangent pose){
	setPose(w,name,pose.exp());
}
*/
string to_string(float i){
	string result;          // string which will contain the result
	ostringstream convert;   // stream used for the conversion
	convert.precision(2);
	convert << i;      // insert the textual representation of 'Number' in the characters in the stream
	result = convert.str();
	return result;
}
int main(int argc, char** argv)
{
	if ( argc != 2 )
		{
	        printf("usage: cTSDF pathToRGBDdataset \n");
	        return -1;
	}
	string basepath=argv[1];
	vector<DepthImage> vdi1(4);
	vector<DepthImage> vdi2(4);
	DepthImage di1(basepath,0);
	vdi1[0]=di1;
	vdi1[1]=vdi1[0].pyrDown(0.5);
	vdi1[2]=vdi1[1].pyrDown(0.5);
	vdi1[3]=vdi1[2].pyrDown(0.5);
	//di1=di1.pyrDown(0.5);
	di1=vdi1[1];
	DepthImage di2(basepath,2);
	vdi2[0]=di2;
	vdi2[1]=vdi2[0].pyrDown(0.5);
	vdi2[2]=vdi2[1].pyrDown(0.5);
	vdi2[3]=vdi2[2].pyrDown(0.5);
	//di2=di2.pyrDown(0.5);
	di2=vdi2[1];
	SE3SolverDirect<S> solver;
	di2.computeGrad();
	solver.setI1(di1);
	solver.setI2(di2);
	imshow("di1",solver.getI1().getImg());
	imshow("di2",solver.getI2().getImg());
	Tangent t(0.0,0,0,0,0,0);
	Mat w=solver.wrapI2(t);
	imshow("wp",w);
	Mat r=solver.RImg(t);
	double min, max;
	minMaxLoc(r, &min, &max);
	cout << "min="<<min << ",max="<<max<<endl;
	Mat rn;
	normalize(r,rn, 0, 1, CV_MINMAX);
	imshow("rn",rn);
	float els =solver.Els (t);
    Mat gels=solver.GEls(t);
    Tangent dTheta=Tangent(gels);
    Tangent dThetaAlpha=dTheta*-0.00001;
    t+=dThetaAlpha;
    cout << "Els=" << els;
    cout << " GEls="<< gels.t()<<endl;
    cout << "dThetaAlpha=" << dThetaAlpha.exp().asMat() <<endl;
    cout << "t=" << t.exp().asMat() <<endl;
	cv::waitKey(0);
	int i=0;
	Mat wpi1;
	while(sqrt(els)>0.001){
		w=solver.wrapI2(t);
		if(i % 2==0){
			cv::addWeighted(di1.getImg(),0.5,w,0.5,0,wpi1);
			imshow("wp",wpi1);
			r=solver.RImg(t);
			normalize(r,rn, 0, 1, CV_MINMAX);
			imshow("rn",rn);
			cout << "i="<<i ;
			cout << " Els=" << sqrt(els);
			cout << " GEls="<< gels.t()<<endl;
			cv::waitKey(1);
		}
		els =solver.Els (t);
		gels=solver.GEls(t);
		dTheta=Tangent(gels);
		dThetaAlpha=dTheta*-0.0002;
		t+=dThetaAlpha;
		//cout << "dThetaAlpha=" << dThetaAlpha.exp().asMat() <<endl;
		//cout << "t=" << t.exp().asMat() <<endl;
		i++;
	}
	cout << "dThetaAlpha=" << dThetaAlpha.exp().asMat() <<endl;
	cout << "t=" << t.exp().asMat() <<endl;
	cout << i << endl;
	waitKey(0);
    /*
	Solver solver;
	viz::Viz3d myWindow("Estimation Coordinate Frame");
	//myWindow.setBackgroundColor(); // black by default
	//myWindow.registerKeyboardCallback(&keyboard_callback);
	myWindow.showWidget("Coordinate Widget", viz::WCoordinateSystem());
    viz::WCube cube_widget(Point3f(1.0,1.0,0.0), Point3f(0.0,0.0,-1.0), true, viz::Color::blue());
    cube_widget.setRenderingProperty(viz::LINE_WIDTH, 1.0);
    myWindow.showWidget("cube",cube_widget);
    viz::WGrid g;
    //viz::WText3D txt("hola Viz",Point3f(0,0,3));
	viz:: WCameraPosition a,b;
	viz:: WCameraPosition wcp0,wcp1,wcpi;
    /// We can get the transformation matrix from camera coordinate system to global using
    /// - makeTransformToGlobal. We need the axes of the camera
    cv::Affine3d transform = viz::makeTransformToGlobal(Vec3f(0.0f,-1.0f,0.0f), Vec3f(-1.0f,0.0f,0.0f), Vec3f(0.0f,0.0f,-1.0f),Vec3f(2,1,0));
    SE3Type se(Mat(transform.matrix));
    Tangent te(se),theta;
    cout << "se" << se.asMat()<<endl;
    cout << "te" << te.exp().asMat()<<endl;
    //cv::Affine3f pose(d.matrix3x4());
    //a.setPose(transform);
    //myWindow.showWidget("txt",txt);
    vector<Point3f> vp1,vp2,vp3;
    vp1.push_back(Point3f(0,1,0));
    vp1.push_back(Point3f(1,1,0));
    vp1.push_back(Point3f(0,1,-1));
    vp1.push_back(Point3f(0.4,1,0.2));
    vp1.push_back(Point3f(1.1,-1,0));
    vp1.push_back(Point3f(-1,1,-1));
    vector<Point3f> P(vp1);
    vp2.resize(vp1.size());
    cout << "Before of transform"<<endl;
    cv::perspectiveTransform(vp1,vp2,transform.matrix);
    //transformPoints(vp2,vp1,theta);
    cout << "vp2"<<vp2 <<endl;
    vp3.resize(vp1.size());
    solver.transformPoints(vp3,vp1,se);
    vector<Point3f> Q(vp3);
    cout << "vp3"<<vp3 <<endl;
    cout << "Els=" << solver.Els(vp3,vp1,se)<<endl;
    cout << "Els=" << solver.Els(vp3,vp1,theta)<<endl;
    //cout << "ri=" << solver.r(3,vp3,vp1,theta)<<endl;
    //cout << "Ri=" << solver.R(3,vp3,vp1,theta)<<endl;
    //cout << "Jri=" << solver.Jr(3,vp3,vp1,theta)<<endl;
    //cout << "GRi=" << solver.GR(3,vp3,vp1,theta)<<endl;
    cout << "GEls="<< solver.GEls(vp3,vp1,theta)<<endl;
    Tangent dTheta(solver.GEls(vp3,vp1,theta));
    Tangent dThetaAlpha=dTheta*-0.05;
    theta+=dThetaAlpha;
    S els=10;
    Mat gels;
    int i=0;
    viz::WCameraPosition *wcp;
    while(sqrt(els)>0.015 && i<10000){
        solver.transformPoints(vp1,vp1,dThetaAlpha);
    	els =solver.Els (vp3,vp1,dThetaAlpha);
        gels=solver.GEls(vp3,vp1,dThetaAlpha);
        dTheta=Tangent(gels);
        dThetaAlpha=dTheta*-0.1;
        theta+=dThetaAlpha;
        if (i%50==0){
        	cout << i << endl;
        	wcp=new viz::WCameraPosition();
            myWindow.showWidget(to_string(i),*wcp);
        	setPose(myWindow,to_string(i),theta);
        	delete wcp;
        }

        cout << "Els=" << els;
        cout << " GEls="<< gels.t()<<endl;
        cout << "dThetaAlpha=" << dThetaAlpha.exp().asMat() <<endl;
        cout << "theta=" << theta.exp().asMat() <<endl;

        i++;
    }
    cout << "vp1" << vp1 << endl;
    cout << "vp3" << vp3 << endl;
    cout << "Els=" << solver.Els(vp3,vp1,dThetaAlpha)<<endl;
    cout << "Sol=" << theta.exp().asMat() <<endl;
    cout << "Soli=" << theta.exp().inverse().asMat() <<endl;

    solver.setAlpha(0.1);
    solver.solveGradientDescent(Q,P,theta);
    cout << "Sol=" << theta.exp().asMat() <<endl;
    solver.setAlpha(0.5);
    solver.solveGaussNewton(Q,P,theta);
    cout << "Sol=" << theta.exp().asMat() <<endl;

    cout << "Before of reshape"<<endl;
	Mat mp1=Mat(vp1).reshape(1);
	Mat mp2=Mat(vp2).reshape(1);
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
	//myWindow.setWidgetPose("a",transform);
	myWindow.setWidgetPose("cube",transform);

    myWindow.spin();
	return 0;
*/
}
