#include <iostream>
#include <cmath>
#include <iomanip>

#define MIN_VALUE 5
#define RANGE 5
#define NEGATIVE_ALPHAS 0.4
#define POSITIVE_ALPHAS 0.25
#define LEARNING_RATE 0.3
#define MIN_ERROR 0.02

using namespace std;

void DisplayVec(float *Vec, int n)
{
    for (int row = 0; row < n; row++)
    {
        cout << Vec[row] << "\n";
    }
}

int GetN()
{
    return MIN_VALUE + (rand() % RANGE);
}

void GetAlpha(float *Alpha, int n)
{
    int positive_values = ceil(float(n) * POSITIVE_ALPHAS);
    int negative_values = ceil(float(n) * NEGATIVE_ALPHAS);

    for (int counter = 0; counter < negative_values; counter++)
    {
        *(Alpha + counter) = -1 * ((float)rand()) / pow(10, 4);
    }

    for (int counter = negative_values; counter < positive_values + negative_values; counter++)
    {
        *(Alpha + counter) = ((float)rand()) / pow(10, 4);
    }

    for (int counter = positive_values + negative_values; counter < n; counter++)
    {
        *(Alpha + counter) = ((float)rand()) / pow(10, 4);
    }
}

void RandomStartingPoints(float *X, int n)
{
    for (int counter = 0; counter < n; counter++)
    {
        *(X + counter) = ((float)rand()) / pow(10, 4);
    }
    cout << endl;
}

float ComputeSum(float *X, float *Alpha, int n)
{
    float sum = 0.0;

    for (int counter = 0; counter < n; counter++)
    {
        sum += Alpha[counter] * X[counter];
    }

    return sum;
}

void UpdateX(float *X, int n)
{
    for (int counter = 0; counter < n; counter++)
    {
        X[counter] -= (int)(LEARNING_RATE * 2 * X[counter]* pow(10,4))/pow(10,4);
    }
}

void MinimizeFunction(float *Alpha, int n)
{
    float X[n] = {0.0};

    RandomStartingPoints(X, n);

    cout << "Starting Values: " << endl;
    DisplayVec(X, n);
    cout << "-------------------------------" << endl;

    float cost = 0.0;
    float error = 0.0;
    int iter = 1;

    float prev_cost = ComputeSum(X, Alpha, n);

    do
    {
        UpdateX(X, n);
        cost = ComputeSum(X, Alpha, n);
        error = fabs(prev_cost - cost);
        prev_cost = cost;

        cout << "Iteration " << iter++ << endl;
        cout << "Error: " << error << endl;
        cout << "Updated X: " << endl;
        DisplayVec(X, n);
        cout << "-------------------------------" << endl;

    } while (error > MIN_ERROR);
}

int main()
{
    srand((unsigned)time(NULL));
    std::cout << std::fixed;
    std::cout << std::setprecision(4);


    int n = GetN();

    cout << "n = " << n << endl;

    float Alpha[n] = {0.0};

    GetAlpha(Alpha, n);

    cout << "Alpha Values: " << endl;
    DisplayVec(Alpha, n);
    cout << "-------------------------------" << endl;

    MinimizeFunction(Alpha, n);
}
