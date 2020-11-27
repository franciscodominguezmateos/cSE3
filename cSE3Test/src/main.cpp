#include <iostream>
#include "pose.h"
using namespace std;

class OpenChainBody{
	vector<ScrewAxis> joints;
	//jacobian columns
	vector<Twist> J;
	Pose home; //"zero" position of the robot
public:
	OpenChainBody(vector<ScrewAxis> tw,Pose h):joints(tw),home(h){
		J.resize(tw.size());
	}
	Pose forwardUpTo(int N,vector<double> &states){
		Pose p;
		for(int i=0;i<N;i++){
			ScrewAxis &sa=joints[i];
			double state=states[i];
			p=p+sa*state;
		}
		return p;
	}
	Pose forward(vector<double> &states){
		Pose p=forwardUpTo(joints.size(),states);
		p=p+home;
		return p;
	}
	void buildJacobianColumns(vector<double> &states){
		for(unsigned int i=0;i<states.size();i++){
			Pose p=forwardUpTo(i,states);
			J[i]=p.Ad(joints[i]);
		}
	}
	Mat jacobian(vector<double> &states){
		buildJacobianColumns(states);
		Mat m=J[0].asMat();
		for(unsigned int i=1;i<J.size();i++){
			Mat mj=J[i].asMat();
			hconcat(m,mj,m);
		}
		return m;
	}
};
void example_3_19(){
	Mat R,t;
	R = (Mat_<double>(3, 3) <<  0, 0,-1,
			                    0,-1, 0,
			                   -1, 0, 0);
	t = (Mat_<double>(3, 1) << 250.0,-150,200.0);
	Pose Tdb(R,t);
	R = (Mat_<double>(3, 3) <<  0, 0,-1,
			                        0,-1, 0,
			                       -1, 0, 0);
	t = (Mat_<double>(3, 1) << 300.0,100,120.0);
	Pose Tde(R,t);
	R = (Mat_<double>(3, 3) <<  0, 0,-1,
			                    0,-1, 0,
			                   -1, 0, 0);
	t = (Mat_<double>(3, 1) << 400.0, 50,300.0);
	Pose Tad(R,t);
	double q2=1/sqrt(2);
	R = (Mat_<double>(3, 3) <<  0, -q2, -q2,
			                    0,  q2, -q2,
			                    1,   0,   0);
	t = (Mat_<double>(3, 1) << 30.0,-40,25.0);
	Pose Tbc(R,t);
	Pose Tae=Tad+Tde;
	Pose Tac=Tad+Tdb+Tbc;
	Pose Taci=-Tac;
	Pose Tce=Taci+Tae;
	cout << "Tae="<< Tae << endl;
	cout << "Tac="<< Tac<<endl;
	cout << "Taci="<< Taci<<endl;
	cout << "Tce="<< Tce<<endl;
	cout << "I="<< Tac+Taci<<endl;
}
void example_3_23(){
	Twist Vs(0,0,2,-2,-4,0);
	Twist Vb(0,0,-2,2.8,4,0);
	Mat R = (Mat_<double>(3, 3) << -1, 0, 0,
			                        0, 1, 0,
			                        0, 0,-1);
	Mat t = (Mat_<double>(3, 1) << 4.0,0.4,0.0);
	Pose Tsb(R,t);
	Twist r=Tsb.Ad(Vb);
	cout << "r="<<r<<endl;
}
void example_3_26(){
	Pose Tsb=Pose(30,1,2);
	Pose Tsc(60,2,1);
	Pose Ti=Tsc+-Tsb;//Tsc=Ti+Tsb
	Twist tw=Ti.log();
	cout << "Tsb="<< Tsb<<endl;
	cout << "Tsc="<<Tsc<<endl;
	cout << "Tsc:"<<Ti+Tsb<<endl;
	cout << "Ti="<<Ti<<endl;
	cout << "Ti="<<Ti.AdMat()<<endl;
	cout << "tw="<<tw<<endl;
	cout << "sa="<<tw.getScrewAxis()<<"theta="<<tw.getTheta()<<endl;
}
void example_3_28(){
	double L1=0.1 ;//10cm
	double L2=0.15;//15cm
	Wrench Fh={0,0,0,0,-5,0};
	Wrench Fa={0,0,0,0, 0,1};
	Mat thf = (Mat_<double>(3, 1) << -L1,0.0,0.0);
	Pose Thf(I,thf);
	Mat Raf = (Mat_<double>(3, 3) <<1, 0, 0,
			                        0, 0, 1,
			                        0,-1, 0);
	Mat taf = (Mat_<double>(3, 1) << -(L1+L2),0.0,0.0);
	Pose Taf(Raf,taf);
	Wrench wFfh=Thf.AdT(Fh);
	Wrench wFfa=Taf.AdT(Fa);
	Wrench wFf=wFfh+wFfa;//Ko¿?
	Mat mFfh=Thf.AdMat().t()*Fh.asMat();
	Mat mFfa=Taf.AdMat().t()*Fa.asMat();
	Mat mFf=mFfh+mFfa;//Ok
	cout << "Thf.AdTMat="<<Thf.AdMat().t() <<endl;
	cout << "Taf.AdTMat="<<Taf.AdMat().t() <<endl;
	cout << "Fh="<<Fh<<endl;
	cout << "Fa="<<Fa<<endl;
	cout << "Thf="<<Thf<<endl;
	cout << "Taf="<<Taf<<endl;
	cout << "mFfh="<<mFfh<<endl;
	cout << "mFfa="<<mFfa<<endl;
	cout << "mFf="<<mFf<<endl;
	cout << "wFfh="<<wFfh<<endl;
	cout << "wFfa="<<wFfa<<endl;
	cout << "wFf="<<wFf<<endl;

}
int main()
{   //example_3_19();
	//example_3_23();
	//example_3_26();
	example_3_28();
	return 0;
}
