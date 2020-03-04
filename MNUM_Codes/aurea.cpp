// Aurea

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double f(double x)
{
	return 5 * cos(x) - sin(x);
}

void aurea()
{
	double B = (sqrt(5) - 1) / 2;
	double A = B * B;

	double x1 = 2, x2 = 4;
	double x3 = x1 + A * (x2 - x1);
	double x4 = x1 + B * (x2 - x1);

	cout << "x1: " << x1 << endl;
	cout << "x2: " << x2 << endl;
	cout << "x3: " << x3 << endl;
	cout << "x4: " << x4 << endl;
	cout << "f(x1): " << f(x1) << endl;
	cout << "f(x2): " << f(x2) << endl;
	cout << "f(x3): " << f(x3) << endl;
	cout << "f(x4): " << f(x4) << endl << endl;

	for (int i = 0; i < 2; ++i)
	{
		if (f(x3) < f(x4))
		{
			x2 = x4;
			x4 = x3;
			x3 = x1 + A * (x2 - x1);
		}
		else if (f(x4) < f(x3))
		{
			x1 = x3;
			x3 = x4;
			x4 = x1 + B * (x2 - x1);
		}

		{
			cout << "x1: " << x1 << endl;
			cout << "x2: " << x2 << endl;
			cout << "x3: " << x3 << endl;
			cout << "x4: " << x4 << endl;
			cout << "f(x1): " << f(x1) << endl;
			cout << "f(x2): " << f(x2) << endl;
			cout << "f(x3): " << f(x3) << endl;
			cout << "f(x4): " << f(x4) << endl << endl;
		}
	}
}

int main()
{
	return 0;
}
