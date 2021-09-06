#include <iostream>
#include <cstring>
#include <chrono>
#include <cmath>

using namespace std;

#define N 1000

double Matrix[N][N] = {};
double B[N] = {};
double X[N] = {};


void DisplayMatrix()
{
	cout << "**************************Matrix***********************" << endl;
	for (int row = 0; row < N; row++)
	{
		for (int col = 0; col < N; col++)
		{
			cout << Matrix[row][col] << "\t";
		}
		cout << "\t|\t" << B[row] << endl;
	}
}

void GenerateEquations()
{
	int sd = 4;

	//Assign random values
	for (int counter1 = 0; counter1 < N; counter1++)
	{
		for (int counter2 = 0; counter2 < N; counter2++)
		{
			Matrix[counter1][counter2] = ((float)rand()) / pow(10, 4);
		}

		B[counter1] = ((float)rand()) / pow(10, 4);
	}
}

void ForwardElimination()
{
	double TempVec[N];

	for (int row = 0; row < N; row++)
	{
		//cout << "Processing Row: " << row << endl;
		int col = row;
		int largestRow = row;

		//Find row with largest value
		for (int row2 = row; row2 < N; row2++)
		{
			if (Matrix[largestRow][col] < Matrix[row2][col])
			{
				largestRow = row2;
			}
		}

		//Swap rows
		memcpy(TempVec, Matrix[row], sizeof(double) * N);
		memcpy(Matrix[row], Matrix[largestRow], sizeof(double) * N);
		memcpy(Matrix[largestRow], TempVec, sizeof(double) * N);

		//DisplayMatrix();

		for (int col2 = 0; col2 < N; col2++)
		{
			Matrix[row][col2] = Matrix[row][col2] / Matrix[row][col];
		}

		for (int row3 = row + 1; row3 < N; row3++)
		{
			for (int col3 = col; col3 < N; col3++)
			{
				double m = Matrix[row3][col3] / Matrix[row][col];
				Matrix[row3][col3] = m * Matrix[row3][col3] - Matrix[row][col3];

				B[row3] = m * B[row3] - B[row];
			}
			Matrix[row3][col] = 0.0;
		}

		//DisplayMatrix();
	}
}

void BackSubstitution()
{
	for (int row = N - 1; row >= 0; row--)
	{
		double temp = 0.0;
		for (int col = row + 1; col < N; col++)
		{
			temp += Matrix[row][col];
		}

		X[row] = B[row] - temp;
	}
}

int main(void)
{
	srand((unsigned)time(NULL));

	//Form matrix
	GenerateEquations();

	//DisplayMatrix();

	// Record start time
	auto start = std::chrono::high_resolution_clock::now();
	ForwardElimination();
	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time for forward elimination: " << elapsed.count() << " sec\n";

	// Record start time
	start = std::chrono::high_resolution_clock::now();
	BackSubstitution();
	// Record end time
	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start;
	std::cout << "Elapsed time for back substitution: " << elapsed.count() << " sec\n";

	cout << endl
		 << "**************************Finished*****************************" << endl;

	return 0;
}
