#include <ctime>
#include <sstream>
#include <iomanip>
#include <iostream>
//#include <cstdlib>

using namespace std;

float A[3][3], b[3];

float error_threshold = 0.01;

float GetRandomNumber()
{
	std::string s;
	do{
		float value = ((float)rand()/(float)(rand()));

		stringstream str;
		str << fixed << setprecision(4) << value;

		s = str.str();
	}while(stof(s) == 0.0);

	return stof(s);
}


bool IsDiagonallyDominant()
{
	bool DiagonallyDominant = false;
	//Check diagonally dominant condition
	if( fabs(A[0][0]) >= (fabs(A[0][1]) + fabs(A[0][2])))
		if( fabs(A[1][1]) >= (fabs(A[1][0]) + fabs(A[1][2])))
			if( fabs(A[2][2]) >= (fabs(A[2][0]) + fabs(A[2][1])))
				DiagonallyDominant = true;

	return DiagonallyDominant;
}



void GenerateEquations()
{
	srand (time(NULL));

	do{

		A[0][0]	= GetRandomNumber();
		A[0][1]	= GetRandomNumber();
		A[0][2]	= GetRandomNumber();

		A[1][0]	= GetRandomNumber();
		A[1][1]	= GetRandomNumber();
		A[1][2]	= GetRandomNumber();

		A[2][0]	= GetRandomNumber();
		A[2][1]	= GetRandomNumber();
		A[2][2]	= GetRandomNumber();

		b[0]	= GetRandomNumber();
		b[1]	= GetRandomNumber();
		b[2]	= GetRandomNumber();

	}while(!IsDiagonallyDominant());

	cout << "A: " << endl;
	cout << A[0][0] << "\t " << A[0][1] << "\t " << A[0][2] << endl;
	cout << A[1][0] << "\t " << A[1][1] << "\t " << A[1][2] << endl;
	cout << A[2][0] << "\t " << A[2][1] << "\t " << A[2][2] << endl;

	cout << "b:" << endl;
	cout << b[0] << "\t " << b[1] << "\t " << b[0] << endl << endl;
}



void GaussSiedelSolution()
{

	cout << "*******************Gauss Siedel Method*******************" << endl;
	float x0_old = 0;
	float x1_old = 0;
	float x2_old = 0;

	float x0_new = (b[0] - A[0][1] * x1_old - A[0][2] * x2_old) / A[0][0];
	float x1_new = (b[1] - A[1][0] * x0_new - A[1][2] * x2_old) / A[1][1];
	float x2_new = (b[2] - A[2][0] * x0_new - A[2][1] * x1_new) / A[2][2];

	x0_old = x0_new;
	x1_old = x1_new;
	x2_old = x2_new;

	float ratio0 = 0.0;
	float ratio1 = 0.0;
	float ratio2 = 0.0;

	float ratio0_old = 0.0;
	float ratio1_old = 0.0;
	float ratio2_old = 0.0;

	bool cond1 = false;
	bool cond2 = false;
	bool cond3 = false;

	int count = 0;

	int divergence_count = 0;
	int divergence_count_threshold = 10;

	//Check if there is divergence condition of new ratio increaing as compared to last ratio.
	do
	{
		count++;

		cout << "Iteration: " << count << endl;
		cout << "\t\tx = " << x0_old << endl;
		cout << "\t\ty = " << x1_old << endl;
		cout << "\t\tz = " << x2_old << endl;

		x0_new = (b[0] - A[0][1] * x1_old - A[0][2] * x2_old) / A[0][0];
		x1_new = (b[1] - A[1][0] * x0_new - A[1][2] * x2_old) / A[1][1];
		x2_new = (b[2] - A[2][0] * x0_new - A[2][1] * x1_new) / A[2][2];

		//compute ratio for error check with 4decimal precision
		ratio0 = floor((fabs(x0_new - x0_old)/fabs(x0_old)) * 10000) / 10000;
		ratio1 = floor((fabs(x1_new - x1_old)/fabs(x1_old)) * 10000) / 10000;
		ratio2 = floor((fabs(x2_new - x2_old)/fabs(x2_old)) * 10000) / 10000;
				
		x0_old = x0_new;
		x1_old = x1_new;
		x2_old = x2_new;

		if(ratio0 >= ratio0_old && ratio1 >= ratio1_old && ratio2 >= ratio2_old)
		{
			divergence_count++;

			if(divergence_count > divergence_count_threshold)
			{
				cout << "Divergence observed" << endl;
				return;
			}
		}
		else
		{
			divergence_count = 0;
		}

		//Update ratios for checking in future
		ratio0_old = ratio0;
		ratio1_old = ratio1;
		ratio2_old = ratio2;

		//Check for error threshold	
		cond1 = ratio0 > error_threshold;
		cond2 = ratio1 > error_threshold;
		cond3 = ratio2 > error_threshold;		

	}while(cond1 || cond1 || cond1);

	cout << "Iteration: " << ++count << endl;
	cout << "\t\tx = " << x0_new << endl;
	cout << "\t\ty = " << x1_new << endl;
	cout << "\t\tz = " << x2_new << endl;

	cout << "\n\nConvergence happened with following values:" << endl;
	cout << "x: " << x0_new << endl;
	cout << "y: " << x1_new << endl;
	cout << "z: " << x2_new << endl;
}



void GaussJacobiSolution()
{

	cout << "*******************Gauss Jacobi Method*******************" << endl;
	float x0_old = 0;
	float x1_old = 0;
	float x2_old = 0;

	float x0_new = (b[0] - A[0][1] * x1_old - A[0][2] * x2_old) / A[0][0];
	float x1_new = (b[1] - A[1][0] * x0_old - A[1][2] * x2_old) / A[1][1];
	float x2_new = (b[2] - A[2][0] * x0_old - A[2][1] * x1_old) / A[2][2];

	x0_old = x0_new;
	x1_old = x1_new;
	x2_old = x2_new;

	float ratio0 = 0.0;
	float ratio1 = 0.0;
	float ratio2 = 0.0;

	float ratio0_old = 0.0;
	float ratio1_old = 0.0;
	float ratio2_old = 0.0;

	bool cond1 = false;
	bool cond2 = false;
	bool cond3 = false;

	int count = 0;

	int divergence_count = 0;
	int divergence_count_threshold = 10;

	do
	{
		count++;

		cout << "Iteration: " << count << endl;
		cout << "\t\tx = " << x0_old << endl;
		cout << "\t\ty = " << x1_old << endl;
		cout << "\t\tz = " << x2_old << endl;

		x0_new = (b[0] - A[0][1] * x1_old - A[0][2] * x2_old) / A[0][0];
		x1_new = (b[1] - A[1][0] * x0_old - A[1][2] * x2_old) / A[1][1];
		x2_new = (b[2] - A[2][0] * x0_old - A[2][1] * x1_old) / A[2][2];

		//compute ratio for error check with 4decimal precision
		ratio0 = floor((fabs(x0_new - x0_old)/fabs(x0_old)) * 10000) / 10000;
		ratio1 = floor((fabs(x1_new - x1_old)/fabs(x1_old)) * 10000) / 10000;
		ratio2 = floor((fabs(x2_new - x2_old)/fabs(x2_old)) * 10000) / 10000;

		x0_old = x0_new;
		x1_old = x1_new;
		x2_old = x2_new;

		//Check if there is divergence condition of new ratio increaing as compared to last ratio.
		if(ratio0 >= ratio0_old && ratio1 >= ratio1_old && ratio2 >= ratio2_old)
		{
			divergence_count++;

			if(divergence_count > divergence_count_threshold)
			{
				cout << "Divergence observed" << endl;
				return;
			}
		}
		else
		{
			divergence_count = 0;
		}

		//Check for error threshold	
		cond1 = ratio0 > error_threshold;
		cond2 = ratio1 > error_threshold;
		cond3 = ratio2 > error_threshold;

		//Update ratios for checking in future
		ratio0_old = ratio0;
		ratio1_old = ratio1;
		ratio2_old = ratio2;

	}while(cond1 || cond1 || cond1);

	cout << "Iteration: " << ++count << endl;
	cout << "\t\tx = " << x0_new << endl;
	cout << "\t\ty = " << x1_new << endl;
	cout << "\t\tz = " << x2_new << endl;

	cout << "\n\nConvergence happened with following values:" << endl;
	cout << "x: " << x0_new << endl;
	cout << "y: " << x1_new << endl;
	cout << "z: " << x2_new << endl;

}


int main(void)
{
	//Form matrix
	GenerateEquations();

	GaussSiedelSolution();

	GaussJacobiSolution();

	getchar();

	return 0;
}

/*
Answers verified from:https://planetcalc.com/3571/
*/