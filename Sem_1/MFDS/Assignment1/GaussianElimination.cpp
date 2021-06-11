#include <time.h>
#include <sstream>
#include <iostream>

using namespace std;

double A[2][2], b[2];
double Matrix[2][2] = {};
double B[2] = {};
//int sd = 2;


double Round(double f, int sd)
{
    char buf[16];
    sprintf(buf, "%.*g", sd, f);

	std::istringstream os(buf);
	double temp;
	os>>temp;

    return temp;
}


void GenerateEquations()
{
	srand ( time(NULL) );

	int sd = 10;

	Matrix[0][0]	= Round((double)rand()/(double)(rand()), sd);
	Matrix[0][1]	= Round((double)rand()/(double)(rand()), sd);
	Matrix[1][0]	= Round((double)rand()/(double)(rand()), sd);
	Matrix[1][1]	= Round((double)rand()/(double)(rand()), sd);

	B[0]	= Round((double)rand()/(double)(rand()), sd);
	B[1]	= Round((double)rand()/(double)(rand()), sd);

	cout << "A: " << endl;
	cout << Matrix[0][0] << "\t " << Matrix[0][1] << endl;
	cout << Matrix[1][0] << "\t " << Matrix[1][1] << endl;

	cout << "B:" << endl;
	cout << B[0] << "\t " << B[1] << endl << endl;
}


void GaussianEliminationWithoutPivoting(int sd)
{
	double m = A[1][0]/A[0][0];
	A[1][0] = 0;
	A[1][1] = Round(A[1][1] - Round(m * A[0][1], sd), sd);

	b[1] = Round(b[1] - Round(m * b[0], sd), sd);
}


void ExchangeRows()
{
	double A00 = A[0][0];
	double A01 = A[0][1];
	double b0 = b[0];

	A[0][0] = A[1][0];
	A[0][1] = A[1][1];
	b[0] = b[1];

	A[1][0] = A00;
	A[1][1] = A01;
	b[1] = b0;


	cout << "After pivoting" << endl;
	cout << "A: " << endl;
	cout << A[0][0] << ", " << A[0][1] << endl;
	cout << A[1][0] << ", " << A[1][1] << endl;

	cout << "b:" << endl;
	cout << b[0] << ", " << b[1] << endl << endl;

}


void GaussianEliminationWithPivoting(int  sd)
{

	if (A[1][0] > A[0][0])
	{
		ExchangeRows();
	}

	GaussianEliminationWithoutPivoting(sd);
}


int main(void)
{
	cout << "Enter significant digits (3,4,5,6)" << endl;
	//cin >> sd;

	/*if (sd > 6 || sd < 3)
	{
		cout << "invalid significant digits" << endl;
		return 0;
	}*/

	//Form matrix
	GenerateEquations();

	for (int i = 3; i < 7; i++)
	{
		cout << "Significant digit: " << i << endl;

		cout << "******************Without Pivoting******************" << endl;

		memcpy(A, Matrix, sizeof(double) * 2 * 2);
		memcpy(b, B, sizeof(double) * 2);

		GaussianEliminationWithoutPivoting(i);

		double y = Round(b[1]/A[1][1], i);
		double x = Round(Round(b[0] - Round(y * A[0][1], i), i)/A[0][0], i);

		cout << "Solution:" << endl;
		cout << "x: " << x << endl;
		cout << "y: " << y << endl;


		cout << endl << "******************With Pivoting******************" << endl;

		memcpy(A, Matrix, sizeof(double) * 2 * 2);
		memcpy(b, B, sizeof(double) * 2);

		GaussianEliminationWithPivoting(i);

		y = Round(b[1]/A[1][1], i);
		x = Round(Round(b[0] - Round(y * A[0][1], i), i)/A[0][0], i);

		cout << "Solution:" << endl;
		cout << "x: " << x << endl;
		cout << "y: " << y << endl;
		
		cout << endl << "*************************************************" << endl;
	}

	getchar();
	return 0;
}

/*
Answers verified from:https://planetcalc.com/3571/
*/