// Regra de Simpson

#include <iostream>
#include <iomanip>
using namespace std;

double simpson(double a, double b, double h)
{
	double result = 0;
	double n = abs(b - a) / h;
	for (unsigned int i = 1; i < n; i += 2)
		result += f(a + i * h) * 4;
	for (unsigned int i = 2; i < n; i += 2)
		result += f(a + i * h) * 2;
	result += f(a) + f(b);
	result *= h / 3;
	return result;
}

double det_QC(double s, double sl, double sll) { return (sl - s) / (sll - sl); }

double erro(double sl, double sll) { return (sll - sl) / 15; }

int main()
{
    return 0;
}