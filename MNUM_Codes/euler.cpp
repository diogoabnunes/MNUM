// Euler

#include <iostream>
#include <iomanip>
using namespace std;

double df(double x, double y) { return x*y; }

void euler(double x0, double y0, double h, unsigned int num_iter)
{
	double xn = x0, yn = y0;
	cout << 0 << setw(16) << xn << setw(16) << yn << setw(16) << df(xn, yn) << endl;
	for (unsigned int i = 1; i < num_iter; i++)
	{
		yn += h * df(xn, yn);
		xn += h;
		cout << i << setw(16) << xn << setw(16) << yn << setw(16) << df(xn, yn) << endl;
	}
}

int main()
{
    return 0;
}