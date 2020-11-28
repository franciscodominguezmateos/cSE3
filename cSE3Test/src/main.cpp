#include <iostream>
#include <iomanip>
#include "open_chain_body.h"
using namespace std;

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
	cout << "Tsb="<<endl<< Tsb<<endl;
	cout << "Tsc="<<endl<<Tsc<<endl;
	cout << "Tsc:"<<endl<<Ti+Tsb<<endl;
	cout << "Ti="<<endl<<Ti<<endl;
	cout << "Ti="<<endl;
	printMat(Ti.AdMat());
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
	cout << "Thf="<<endl<<Thf<<endl;
	cout << "Taf="<<endl<<Taf<<endl;
	cout << "mFfh="<<mFfh<<endl;
	cout << "mFfa="<<mFfa<<endl;
	cout << "mFf="<<mFf<<endl;
	cout << "wFfh="<<wFfh<<endl;
	cout << "wFfa="<<wFfa<<endl;
	cout << "wFf="<<wFf<<endl;

}
void example_4_6(){
	double L1=0.425;
	double L2=0.392;
	double W1=0.109;
	double W2=0.082;
	double H1=0.089;
	double H2=0.095;
	Mat RM = (Mat_<double>(3, 3) <<-1, 0, 0,
			                        0, 0, 1,
			                        0, 1, 0);
	Mat tM = (Mat_<double>(3, 1) << L1+L2,W1+W2,H1-H2);
	Pose M(RM,tM);
	vector<ScrewAxis> axis={
			{0,0, 1,        0,    0,0    },
			{0,1, 0,      -H1,    0,0    },
			{0,1, 0,      -H1,    0,L1   },
			{0,1, 0,      -H1,    0,L1+L2},
			{0,0,-1,      -W1,L1+L2,0    },
			{0,1, 0,    H2-H1,    0,L1+L2}
	};
	OpenChainBody USR6R(axis,M);
	vector<double> state(axis.size());
	state[1]=-PI/2.0;
	state[4]= PI/2.0;
	cout << "T="<<endl<< USR6R.forward(state)<<endl;
	Pose eS2=Pose::exp(axis[1]*state[1]);
	cout << "eS2="<<endl<<eS2<<endl;
}
void example_4_7(){
	//It is end-effector frame
	double L1=0.550;
	double L2=0.300;
	double L3=0.060;
	double W1=0.045;
	Mat tM = (Mat_<double>(3, 1) << 0.0,0.0,L1+L2+L3);
	Pose M(I,tM);
	vector<ScrewAxis> axis={
			{0,0,1,        0,0,0},
			{0,1,0, L1+L2+L3,0,0},
			{0,0,1,        0,0,0},
			{0,1,0,    L2+L3,0,W1},
			{0,0,1,        0,0,0},
			{0,1,0,       L3,0,0},
			{0,0,1,        0,0,0}
	};
	OpenChainBody Wam7R(axis,M);
	vector<double> state(7);
	state[1]= PI/4.0;
	state[3]=-PI/4.0;
	state[5]=-PI/2.0;
	cout << "T="<<endl<< Wam7R.forwardBody(state)<<endl;
}
void example_6_1(){
	//It is end-effector frame
	double L1=1;
	double L2=1;
	Mat tM = (Mat_<double>(3, 1) << L1+L2,0.0,0.0);
	Pose M(I,tM);
	vector<ScrewAxis> axis={
			{0,0,1, 0,L1+L2,0},
			{0,0,1, 0,L1   ,0}
	};
	OpenChainBody RR(axis,M);
	vector<double> state(axis.size());
	state[0]= 0.0;
	state[1]= 30.0/180.0*PI;
	cout << "T="<<endl<< RR.forwardBody(state)<<endl;
	Mat Rsd = (Mat_<double>(3, 3) <<-0.5  , -0.866, 0,
			                        0.866, -0.5  , 1,
			                        0    ,  1    , 0);
	Mat tsd = (Mat_<double>(3, 1) << 0.366,1.366,0);
	Pose Tsd(Rsd,tsd);
	cout << "Tsd="<<endl<<Tsd<<endl;
	state=RR.ikStep(state,Tsd.log());
	for(double d:state) cout <<d<<",";
	cout<<endl;
	cout << "T="<<endl<<RR.forwardBody(state)<<endl;
}
int main(){
	cout.setf(ios::fixed);
	cout.setf(ios::showpoint);
	cout.precision(3);
   //example_3_19();
	//example_3_23();
	//example_3_26();
	//example_3_28();
	//example_4_6();
	//example_4_7();
	example_6_1();
	return 0;
}
