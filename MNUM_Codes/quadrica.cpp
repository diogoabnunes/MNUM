// Quadrica

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

double f(double x, double y)
{
	return y * y - 2 * x * y - 6 * y + 2 * x * x + 12;
}

double dfx(double x, double y)
{
	return -2 * y + 4 * x;
}

double dfy(double x, double y)
{
	return 2 * y - 2 * x - 6;
}

void quadrica(double f(double, double), double x, double y, int numIt)
{
	vector<vector<double>> H = {{0.5, 0.5}, {0.5, 1}}; //matriz hessiana inversa

	double next_x = x - (H.at(0).at(0) * dfx(x, y) + H.at(0).at(1) * dfy(x, y));
	double next_y = y - (H.at(1).at(0) * dfx(x, y) + H.at(1).at(1) * dfy(x, y));

	for (int i = 0; i < numIt; ++i)
	{
		cout << "x = " << next_x << endl;
		cout << "y = " << next_y << endl;
		cout << "f(x,y) = " << f(next_x, next_y) << endl;

		if (f(next_x, next_y) < f(x, y))
		{
			x = next_x;
			y = next_y;
		}

		next_x = x - (H.at(0).at(0) * dfx(x, y) + H.at(0).at(1) * dfy(x, y));
		next_y = y - (H.at(1).at(0) * dfx(x, y) + H.at(1).at(1) * dfy(x, y));
	}
}

int main()
{
    return 0;
}