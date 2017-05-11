//===========================================================================
//Final_Project.cpp
//Particle Swarm Optimization
//
//Input:
//	N/A
//			
//Output:
//	Path.txt	
//
//===========================================================================

// #define NDEBUG
#define cl4ptp (double)rand()/RAND_MAX
#include "stdafx.h"
#include "stdafx.h"
#include <conio.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cassert>
#include <random>
#include <vector>
#include <iomanip>
#include <deque>
#include <algorithm>
#include <numeric>
#include <functional>
#include <Windows.h>

using namespace std;

// Classes
class bird
{
public:
	double y1;	// Self Adjustment Weight ( 0 < y1 < 1 )
	double y2;	// Social Adjustment Weight ( 0 < y2 < 1 )
	double W;	// Inertia ( 0 < W < 1 )
	double fit;
	vector<double> personal_max;
	vector<double> global_max;
	vector<double> x;
	vector<double> y;
	vector<double> v;
	

	void init(double max, double v_max);
	double fitness();
	void update(double resolution, double max);
};

// Function Prototypes
void print(vector<bird> flock, int global_bird_best, vector<double> learning_curve);

//===========================================================================
//===========================================================================					main
int main()
{
	// Initialize Random Number Generator
	srand((unsigned)time(NULL));

	// Initialize Birds
	int bird_num = 4;
	double global_m = -100000;
	int global_bird_best = 0;
	vector<double> learning_curve;
	double max = 5;		// Starting Grid is a (2*max)*(2*max) square centered at 0. Sets initialization limits.
	double v_max = 2;	// Sets maximum starting velocity
	vector<bird> flock;
	for (int i = 0; i < bird_num; i++)
	{
		bird a;
		a.init(max, v_max);
		double temp = a.fitness();
		if (temp > global_m)
		{
			global_m = temp;
			global_bird_best = i;
		}
		flock.push_back(a);
	}
	learning_curve.push_back(global_m);

	// Set Global Max for all birds
	for (int i = 0; i < bird_num; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			flock.at(i).global_max.at(j) = flock.at(global_bird_best).personal_max.at(j);
		}
	}

	// Perform Optimization
	int iteration_num = 5000;
	double resolution = .1;
	for (int jump_num = 0; jump_num < iteration_num; jump_num++)
	{
		// Update All Bird Locations, Velocities, and Fitnesses
		for (int i = 0; i < bird_num; i++)
		{
			flock.at(i).update(resolution, max);
		}

		// Find Best Bird
		for (int i = 0; i < bird_num; i++)
		{
			if (flock.at(i).personal_max.at(2) > global_m)
			{
				global_m = flock.at(i).personal_max.at(2);
				global_bird_best = i;
			}
		}
		//cout << global_m << '\t' << flock.at(global_bird_best).personal_max.at(0) << '\t' << flock.at(global_bird_best).personal_max.at(1) <<  endl;
		//Sleep(25);

		// Set Global Max for all birds
		for (int i = 0; i < bird_num; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				flock.at(i).global_max.at(j) = flock.at(global_bird_best).personal_max.at(j);
			}
		}
		learning_curve.push_back(global_m);
	}
	cout << global_m << endl;

	// Print to file
	print(flock, global_bird_best, learning_curve);

    return 0;
}

//===========================================================================					print
void print(vector<bird> flock, int global_bird_best, vector<double> learning_curve)
{
	// Create output file
	ofstream output_file;
	output_file.open("Path.txt");
	int end = flock.at(global_bird_best).x.size();
	int big_end = flock.size();

	for (int i = 0; i < end; i++)
	{
		for (int j = 0; j < big_end; j++)
		{
			output_file << flock.at(j).x.at(i) << '\t' << flock.at(j).y.at(i) << '\t';
		}
		output_file << "\tSPACE\t" << learning_curve.at(i) << endl;
	}

	//Close output file
	output_file.close();

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++					bird::init
void bird::init(double max, double v_max)
{
	x.push_back(cl4ptp * 2 * max - max);	// X location (-max < x < max)
	y.push_back(cl4ptp * 2 * max - max);	// Y location (-max < y < max)
	v.push_back(cl4ptp * 2 * v_max - v_max);	// Speed in the x
	v.push_back(cl4ptp * 2 * v_max - v_max);	// Speed in the y

	y1 = 1.5;		// Self Adjustment Weight ( 0 < y1 < 1 )
	y2 = 1.1;	// Social Adjustment Weight ( 0 < y2 < 1 )
	W  = .9;	// Inertia ( 0 < W < 1 )

	personal_max.push_back(x.back());
	personal_max.push_back(y.back());
	personal_max.push_back(fitness());

	for (int i = 0; i < 3; i++)
	{
		global_max.push_back(0);
	}

	fit = fitness();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++					bird::fitness
double bird::fitness()
{
	double a = .5; // Change number of local maxima, lower numbers increase local maxima density
	double fitness = cos(x.back() / a) * (1 - .2 * fabs(x.back()))+ cos(y.back() / a)*(1 - .2 * fabs(y.back())) + 2;
	// Global Max at (0,0)
	return fitness;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++					bird::update
void bird::update(double resolution, double max)
{
	double vx = v.at(0);
	double vy = v.at(1);
	double r1 = cl4ptp;
	double r2 = cl4ptp;

	vx = W * vx + y1 * r1 * (personal_max.at(0) - x.back()) + y2 * r2 * (global_max.at(0) - x.back());
	vy = W * vy + y1 * r1 * (personal_max.at(1) - y.back()) + y2 * r2 * (global_max.at(1) - y.back());

	x.push_back(x.back() + vx * resolution);
	if (x.back() > max)
	{
		x.pop_back();
		x.push_back(max);
	}
	if (x.back() < -1 * max)
	{
		x.pop_back();
		x.push_back(-1 * max);
	}
	y.push_back(y.back() + vy * resolution);
	if (y.back() > max)
	{
		y.pop_back();
		y.push_back(max);
	}
	if (y.back() < -1 * max)
	{
		y.pop_back();
		y.push_back(-1 * max);
	}

	v.at(0) = vx;
	v.at(1) = vy;

	double fit_temp = fitness();
	if (fit_temp > personal_max.at(2))
	{
		personal_max.at(0) = x.back();
		personal_max.at(1) = y.back();
		personal_max.at(2) = fit_temp;
	}

}