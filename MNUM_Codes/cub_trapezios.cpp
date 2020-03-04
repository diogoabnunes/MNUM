// Regra dos Trapézios - Cubatura

#include <iostream>
#include <iomanip>
using namespace std;

double f(double x, double y) { return 2; }

double trapezio_c(double a, double A, double b, double B)
{
	double result = 0;
	for (unsigned int x = 0; x < 3; x++)
	{
		for (unsigned int y = 0; y < 3; y++)
		{
			if ((x == 0 || x == 2) && (y == 0 || y == 2)) result += f(x, y);
			else if (x == 0 || x == 2 || y == 0 || y == 2) result += f(x, y) * 2;
			else result += f(x, y) * 4;
		}
	}
	result /= 4;
	return result;
}

int main()
{
    return 0;
}