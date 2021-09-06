#include <iostream>
#include <chrono>
#include <cmath>

using namespace std;

#define N 1000

float A[N][N] = {};
float L[N][N] = {};
float U[N][N] = {};
float B[N] = {};
float V[N] = {};
float X[N] = {};

int Lower_Triangular_Inverse(float *L, int n)
{
    int i, j, k;
    float *p_i, *p_j, *p_k;
    float sum;

    //Invert the diagonal elements of the lower triangular matrix L.

    for (k = 0, p_k = L; k < n; p_k += (n + 1), k++)
    {
        if (*p_k == 0.0)
            return -1;
        else
            *p_k = 1.0 / *p_k;
    }

    //Invert the remaining lower triangular matrix L row by row.

    for (i = 1, p_i = L + n; i < n; i++, p_i += n)
    {
        for (j = 0, p_j = L; j < i; p_j += n, j++)
        {
            sum = 0.0;
            for (k = j, p_k = p_j; k < i; k++, p_k += n)
                sum += *(p_i + k) * *(p_k + j);
            *(p_i + j) = -*(p_i + i) * sum;
        }
    }

    return 0;
}

int Upper_Triangular_Inverse(float *U, int n)
{
    int i, j, k;
    float *p_i, *p_j, *p_k;
    float sum;

    //Invert the diagonal elements of the upper triangular matrix U.

    for (k = 0, p_k = U; k < n; p_k += (n + 1), k++)
    {
        if (*p_k == 0.0)
            return -1;
        else
            *p_k = 1.0 / *p_k;
    }

    //Invert the remaining upper triangular matrix U.

    for (i = n - 2, p_i = U + n * (n - 2); i >= 0; p_i -= n, i--)
    {
        for (j = n - 1; j > i; j--)
        {
            sum = 0.0;
            for (k = i + 1, p_k = p_i + n; k <= j; p_k += n, k++)
            {
                sum += *(p_i + k) * *(p_k + j);
            }
            *(p_i + j) = -*(p_i + i) * sum;
        }
    }

    return 0;
}

void DisplayMatrix(float Matrix[N][N])
{
    for (int row = 0; row < N; row++)
    {
        for (int col = 0; col < N; col++)
        {
            cout << Matrix[row][col] << "\t";
        }
        //cout << "\t|\t" << B[row] << endl;
        cout << endl;
    }
    cout << "*************************************************" << endl;
}

void DisplayVec(float Vec[N])
{
    for (int row = 0; row < N; row++)
    {
        cout << Vec[row] << "\n";
    }
    cout << "*************************************************" << endl;
}

void CroutsMethod()
{
    int i, j, k;
    float sum = 0;
    int n = N;

    for (i = 0; i < n; i++)
    {
        U[i][i] = 1;
    }

    for (j = 0; j < n; j++)
    {
        for (i = j; i < n; i++)
        {
            sum = 0;
            for (k = 0; k < j; k++)
            {
                sum = sum + L[i][k] * U[k][j];
            }
            L[i][j] = A[i][j] - sum;
        }

        for (i = j; i < n; i++)
        {
            sum = 0;
            for (k = 0; k < j; k++)
            {
                sum = sum + L[j][k] * U[k][i];
            }
            if (L[j][j] == 0)
            {
                printf("det(L) close to 0!\n Can't divide by 0...\n");
                exit(EXIT_FAILURE);
            }
            U[j][i] = (A[j][i] - sum) / L[j][j];
        }
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
            A[counter1][counter2] = ((float)rand()) / pow(10, 4);
        }

        B[counter1] = ((float)rand()) / pow(10, 4);
    }
}

void ForwardElimination()
{
    Lower_Triangular_Inverse(&(L[0][0]), N);
    //cout << "Linv:" << endl;
    //DisplayMatrix(L);

    for (int row = 0; row < N; row++)
    {
        float temp = 0.0;
        for (int col = 0; col < N; col++)
        {
            temp += L[row][col] * B[col];
        }
        V[row] = temp;
    }

    //cout << "V:" << endl;
    //DisplayVec(V);
}

void BackSubstitution()
{
    Upper_Triangular_Inverse(&U[0][0], N);
    //cout << "Uinv:" << endl;
    //DisplayMatrix(U);

    for (int row = 0; row < N; row++)
    {
        float temp = 0.0;
        for (int col = 0; col < N; col++)
        {
            temp += U[row][col] * V[col];
        }
        X[row] = temp;
    }
}

int main(void)
{
    srand((unsigned)time(NULL));

    //Form matrix
    GenerateEquations();

    //cout << "A:" << endl;
    //DisplayMatrix(A);
    //cout << "B:" << endl;
    //DisplayVec(B);

    CroutsMethod();

    //cout << "L:" << endl;
    //DisplayMatrix(L);
    //cout << "U:" << endl;
    //DisplayMatrix(U);

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

    //cout << "X:" << endl;
    //DisplayVec(X);

    cout << endl
         << "**************************Finished*****************************" << endl;

    return 0;
}

