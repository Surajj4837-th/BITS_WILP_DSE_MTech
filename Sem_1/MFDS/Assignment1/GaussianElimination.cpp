#include <time.h>
#include <sstream>
#include <iostream>

using namespace std;

double A[2][2], b[2];
int sd = 2;


double Round(double f)
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

	A[0][0]	= Round((double)rand()/(double)(rand()));
	A[0][1]	= Round((double)rand()/(double)(rand()));
	A[1][0]	= Round((double)rand()/(double)(rand()));
	A[1][1]	= Round((double)rand()/(double)(rand()));

	b[0]	= Round((double)rand()/(double)(rand()));
	b[1]	= Round((double)rand()/(double)(rand()));

	cout << "A: " << endl;
	cout << A[0][0] << "\t " << A[0][1] << endl;
	cout << A[1][0] << "\t " << A[1][1] << endl;

	cout << "b:" << endl;
	cout << b[0] << "\t " << b[1] << endl << endl;
}


void GaussianEliminationWithoutPivoting()
{
	double m = A[1][0]/A[0][0];
	A[1][0] = 0;
	A[1][1] = Round(A[1][1] - Round(m * A[0][1]));

	b[1] = Round(b[1] - Round(m * b[0]));
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


void GaussianEliminationWithPivoting()
{

	if (A[1][0] > A[0][0])
	{
		ExchangeRows();
	}

	GaussianEliminationWithoutPivoting();
}


int main(void)
{
	cout << "Enter significant digits (2,3,4,5,6)" << endl;
	cin >> sd;

	cout << "Enter 1 for Gaussian Elimination without pivoting else enter 2." << endl;
	int choice = 0;
	cin >> choice;

	if (choice < 1 || choice > 2)
	{
		cout << "invalid pivoting choice." << endl;
		return 0;
	}

	if (sd > 6 || sd < 2)
	{
		cout << "invalid significant digits" << endl;
		return 0;
	}	
	
	//Form matrix
	GenerateEquations();

	if (choice == 1)
		GaussianEliminationWithoutPivoting();
	else
		GaussianEliminationWithPivoting();

	double y = Round(b[1]/A[1][1]);
	double x = Round(Round(b[0] - Round(y * A[0][1]))/A[0][0]);

	cout << "Solution:" << endl;
	cout << "x: " << x << endl;
	cout << "y: " << y << endl;
	
	getchar();
	return 0;
}

/*
Answers verified from:https://planetcalc.com/3571/
*/