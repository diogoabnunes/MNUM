// Regra dos Trapézios

#include <iostream>
#include <iomanip>
using namespace std;

double trapezio(double a, double b, double h)
{
	double result = 0;
	double n = abs(b - a) / h;
	for (unsigned int i = 1; i < n; i++)
		result += f(a + i * h) * 2;
	result += f(a) + f(b);
	result *= h / 2;
	return result;
}

double det_QC(double s, double sl, double sll) { return (sl - s) / (sll - sl); }

double erro(double sl, double sll) { return (sll - sl) / 15; }

int main()
{
    return 0;
}