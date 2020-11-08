#include "Connection.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <kdl/jntarray.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolvervel_wdls.hpp>
#include <kdl/chainiksolverpos_nr_jl.hpp>
#include <kdl/utilities/svd_eigen_HH.hpp>
#include <eigen3/Eigen/Core>
#include <urdf_model/model.h>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <vector>

#define FILE_SIZE 1000 * 1000

using namespace std;
using namespace KDL;
using namespace boost;


void doTrajectory();
void vecNormalRandom(vector<double> &q);
void transform(JntArray jnt_arr, vector<double> &vec);
void transform(Frame frame, vector<double> &vec);
void transform(vector<double> vec, JntArray &jnt_arr);
void transform(vector<double> vec, Frame &frame);
void print(JntArray j);
void print(vector<double> vec);
void print(Frame f);

bool explanation = true;
//Initial value of ang_init is: 0, 0, 0.
vector<double> ang_init(3), ang_targ(3), pnt_targ(12);

class Robot {
	Chain chain;
	ChainFkSolverPos_recursive *fk_solver;
	ChainIkSolverVel_wdls *ik_solver;
	JntArray jnt_min, jnt_max;
	unsigned int max_iter;
	unsigned int nj;
	double eps;

public:
	Robot();
	~Robot();
	void solveFK(vector<double> &q_in, vector<double> &p_out);
	int coreFK(JntArray &q_in, Frame &p_out);
	bool solveIK(vector<double> &q_init, vector<double> &q_targ, vector<double> &p_targ);
	bool clamp(vector<double> &q_targ);
}robot;


int main (int argc, char *argv[]) {
	doTrajectory();

	return 0;
}

void doTrajectory () {
	char buffer[FILE_SIZE];
	FILE *file;
	size_t nRead;

	file = fopen("input.in", "r");
	if (file) {
   		while ((nRead = fread(buffer, 1, sizeof buffer, file)) > 0) {
   			stringstream stream(buffer);

   			int cnt = 3;
   			vector<double> coord;
   			double inp, lastTime = 0.;
    		while(stream >> inp) {
    			if (cnt--) 
    				coord.push_back(inp);
    			else {
    				sleep(inp - lastTime);
    				lastTime = inp;

    				transform(Frame(Vector(coord[0], coord[1], coord[2])), pnt_targ);
    				robot.solveIK(ang_init, ang_targ, pnt_targ);

    				//TODO: send(ang_targ); // ang_targ should be casted and encoded to desired type.

    				//TODO: receive(ang_init); // input with desired type should be casted and decoded to ang_init

    				cnt = 3;
    				coord.clear();
    			}
    		}
   		}
        	
    	fclose(file);
	}
}


Robot::Robot () {
	//Init
	max_iter = 10000;
	eps = 1e-6;

	Joint joint0(Joint::None);
	Frame frame0 = Frame(Vector(0.0, 0.0, 0.0));
	chain.addSegment(Segment(joint0, frame0));

	Joint joint1(Joint::RotZ);
	Frame frame1 = Frame::DH(10, M_PI/2., 0, 0);
	chain.addSegment(Segment(joint1, frame1));

	Joint joint2(Joint::RotZ);
	Frame frame2 = Frame::DH(5, 0, 0, 0);
	chain.addSegment(Segment(joint2, frame2));

	Joint joint3(Joint::RotZ);
	Frame frame3 = Frame::DH(5, 0, 0, 0);
	chain.addSegment(Segment(joint3, frame3));

	nj = chain.getNrOfJoints();
	jnt_min.resize(nj);
	jnt_max.resize(nj);


	jnt_min(0) = -M_PI;
	jnt_max(0) = M_PI;
	jnt_min(1) = -M_PI/2.;
	jnt_max(1) = M_PI/2.;
	jnt_min(2) = -M_PI;
	jnt_max(2) = M_PI;
	
	fk_solver = new ChainFkSolverPos_recursive(chain);
	ik_solver = new ChainIkSolverVel_wdls(chain, eps, max_iter);
}


Robot::~Robot () {
	delete fk_solver;
	delete ik_solver;
}


void Robot::solveFK (vector<double> &q_in, vector<double> &p_out) {
	JntArray q_in_jnt_array(nj);
	Frame p_out_frame;
	transform(q_in, q_in_jnt_array);
	if(coreFK(q_in_jnt_array, p_out_frame) >= 0) 
		transform(p_out_frame, p_out);
	else if (explanation)
		printf("Something went wrong while calculating FK.\n");
}

int Robot::coreFK (JntArray &q_in, Frame &p_out) {
	return fk_solver->JntToCart(q_in, p_out);
}


bool Robot::solveIK (vector<double> &q_init, vector<double> &q_targ, vector<double> &p_targ) {
	JntArray q_init_jnt_arr(nj);
	JntArray q_targ_jnt_arr(nj);
	JntArray q_delta(nj);
	Frame p_init_frame;
	Frame p_targ_frame;
	Twist p_delta;

	transform(q_init, q_init_jnt_arr);
	transform(p_targ, p_targ_frame);

	double k2 = 0.01;
	bool found = false;

	int trialCnt = 0;
	while(!found && trialCnt < 100) {
		q_targ_jnt_arr = q_init_jnt_arr;
		for(int iter = 0; iter < max_iter; ++iter) {
			//Get initial pose of EE
			coreFK(q_targ_jnt_arr, p_init_frame); 			            
			
			//Calculate the twist between init and final
			p_delta = diff(p_init_frame, p_targ_frame);				    
			
			//Calculate the joint angle change rate
			ik_solver->CartToJnt(q_init_jnt_arr, p_delta, q_delta);	    
			Multiply(q_delta, k2, q_delta);
			Add(q_targ_jnt_arr, q_delta, q_targ_jnt_arr);

			if(Equal(p_delta,Twist::Zero(),eps)) {
				transform(q_targ_jnt_arr, q_targ);
				found = clamp(q_targ);
				if(explanation && found) {
					printf("IK solution found at iteration %d.\n", iter);
					break;
				}
				else if(found)
					break;
			}
			else if(iter == max_iter-1 && explanation)
				printf("Maximum iteration count reached. ");
		}

		if(!found) {
    		vector<double> q(nj);
    		transform(q_init_jnt_arr, q);
    		vecNormalRandom(q);
    		transform(q, q_init_jnt_arr);
    		++trialCnt;
    		if(explanation) {
    			printf("Trial %d finished. Trying with new initial joint positions:\n", trialCnt);
    			print(q_init_jnt_arr);
    		}
    	}		    	
	}
	if(!found && explanation)
		printf("Could not find a valid IK solution.\n");

	return found;
}

bool Robot::clamp (vector<double> &q_targ) {
	bool found = true;
	for(int i = 0; i < q_targ.size(); ++i) {
		if(jnt_min(i) <= fmod(q_targ[i],2*M_PI) && fmod(q_targ[i],2*M_PI) <= jnt_max(i))
			q_targ[i] = fmod(q_targ[i],2*M_PI);
		else if(jnt_min(i) <= fmod(q_targ[i],2*M_PI)-2*M_PI && fmod(q_targ[i],2*M_PI)-2*M_PI <= jnt_max(i))
			q_targ[i] = fmod(q_targ[i],2*M_PI)-2*M_PI;
		else if(jnt_min(i) <= fmod(q_targ[i],2*M_PI)+2*M_PI && fmod(q_targ[i],2*M_PI)+2*M_PI <= jnt_max(i))
			q_targ[i] = fmod(q_targ[i],2*M_PI)+2*M_PI;
		else {
			if(q_targ[i] < jnt_min(i))
				q_targ[i] = jnt_min(i);
			else
				q_targ[i] = jnt_max(i);
			
			found = false;
		}
	}
	return found;		
}

void vecNormalRandom (vector<double> &q) {
	srand(time(NULL));
	for(int i = 0; i < q.size(); ++i) {
		mt19937 *rng = new mt19937();
		rng->seed(time(NULL));
		normal_distribution<> distribution(q[i], 0.001);
		variate_generator< mt19937, normal_distribution<> > dist(*rng, distribution);
		double new_angle = dist(); 
		q[i] = new_angle;
	}
}


//Transform joint array to double vector
void transform (JntArray jnt_arr, vector<double> &vec) {
	int len = jnt_arr.data.size();
	for(int i = 0; i < len; ++i)
		vec[i] = jnt_arr.data(i);
}

//Transform frame to double vector
void transform (Frame frame, vector<double> &vec) {
	for(int i = 0; i < 3; ++i)
		vec[i] = frame.p.data[i];

	for(int i = 0; i < 9; ++i)
      	vec[i+3] = frame.M.data[i];
}

//Transform double vector to joint array
void transform (vector<double> vec, JntArray &jnt_arr) {
	for(int i = 0; i < vec.size(); ++i)
		jnt_arr.data(i) = vec[i];
}

//Transform double vector to frame
void transform (vector<double> vec, Frame &frame) {
	for(int i = 0; i < 3; ++i)
   		frame.p.data[i] = vec[i];

	for(int i = 0; i < 9; ++i)
 		frame.M.data[i] = vec[i+3];
}

//Print joint array
void print (JntArray j) {
	vector<double> print_vec(j.data.size());
	transform(j, print_vec);
	print(print_vec);
}

//Print double vector
void print (vector<double> vec) {
	int len = vec.size();
	printf("[");
	for(int i = 0; i < len; ++i)
		printf(" %.5f ", vec[i]);
	printf("]\n");
}

//Print frame
void print (Frame f) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
			printf("%f ", f(i, j));
		printf("\n");
	}
	printf("\n");
}
