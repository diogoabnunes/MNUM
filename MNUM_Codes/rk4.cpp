// Range-Kutta 4

#include <iostream>
#include <iomanip>
using namespace std;

double df(double x, double y) { return x*y; }

void rk4(double x0, double y0, double h, unsigned int num_iter)
{
	double xn = x0, yn = y0;
	cout << 0 << setw(16) << xn << setw(16) << yn << setw(16) << df(xn, yn) << endl;
	for (unsigned int i = 1; i < num_iter; i++)
	{
		double dy1 = h * df(xn, yn);
		double dy2 = h * df(xn + h / 2, yn + dy1 / 2);
		double dy3 = h * df(xn + h / 2, yn + dy2 / 2);
		double dy4 = h * df(xn + h, yn + dy3);
		yn += dy1 / 6 + dy2 / 3 + dy3 / 3 + dy4 / 6;
		xn += h;
		cout << i << setw(16) << xn << setw(16) << yn << setw(16) << df(xn, yn) << endl;
	}
}

int main()
{
    return 0;
}