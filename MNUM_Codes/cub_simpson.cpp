// Regra de Simpson - Cubatura

#include <iostream>
#include <iomanip>
using namespace std;

double f(double x, double y) { return 2; }

double simpson_c(double a, double A, double b, double B)
{
	double hx = abs(A - a) / 2, hy = abs(B - b) / 2;
	double result = 0;
	for (unsigned int y = 0; y < 3; y++)
	{
		double temp = 0;
		for (unsigned int x = 0; x < 3; x++)
		{
			if (x == 1) temp += f(x, y) * 4;
			else temp += f(x, y);
		}
		if (y == 1) result += temp * 4;
		else result += temp;
	}
	result *= (hx * hy) / 9;
	return result;
}

int main()
{
    return 0;
}