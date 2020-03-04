// Gauss-Seidel

#include <iostream>
#include <iomanip>
using namespace std;

/*
Converge: Sim, porque em cada linha da matriz A, o módulo do elemento
da diagonal principal é superior ao módulo da soma dos restantes elementos da linha.
*/

void gauss_seidel(matrix& m, double x1, double x2, double x3, double x4, unsigned int num_iter)
{
	cout << 0 << ":\t" << x1 << "\t" << x2 << "\t" << x3 << "\t" << x4 << endl;
	for (unsigned int i = 1; i <= num_iter; i++)
	{
		x1 = (m[0][4] - (m[0][1] * x2 + m[0][2] * x3 + m[0][3] * x4)) / m[0][0];
		x2 = (m[1][4] - (m[1][0] * x1 + m[1][2] * x3 + m[1][3] * x4)) / m[1][1];
		x3 = (m[2][4] - (m[2][0] * x1 + m[2][1] * x2 + m[2][3] * x4)) / m[2][2];
		x4 = (m[3][4] - (m[3][0] * x1 + m[3][1] * x2 + m[3][2] * x3)) / m[3][3];
		cout << i << ":\t" << x1 << "\t" << x2 << "\t" << x3 << "\t" << x4 << endl;
	}
}

int main()
{
    return 0;
}