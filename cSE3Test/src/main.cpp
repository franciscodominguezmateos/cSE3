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
			J[i]=p.adjoint(joints[i]);
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
void example3_19(){
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
int main()
{
	Twist Vs(0,0,2,-2,-4,0);
	Twist Vb(0,0,-2,2.8,4,0);
	Mat R = (Mat_<double>(3, 3) << -1, 0, 0,
			                        0, 1, 0,
			                        0, 0,-1);
	Mat t = (Mat_<double>(3, 1) << 4.0,0.4,0.0);
	//Pose Tsb(R,t);
	//Twist r=Tsb.adjoint(Vb);
	Pose Tsb(30,1,2);
	Pose Tsc(60,2,1);
	Pose Ti=Tsc+-Tsb;
	Twist tw=Ti.log();
	cout << "Tsb="<< Tsb.asMat()<<endl;
	cout << "Tsc="<<Tsc.asMat()<<endl;
	cout << "Tsc:"<<Ti+Tsb<<endl;
	cout << "Ti="<<Ti<<endl;
	cout << "Ti="<<Ti.adjointMat()<<endl;
	cout << "tw="<<tw<<endl;
	cout << "sa="<<tw.getScrewAxis()<<"theta="<<tw.getTheta()<<endl;
	//example3_19();
	return 0;
}
