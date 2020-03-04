// Levenberg-Marquardt

#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

double f(double x, double y)
{
	return (x + 1) * (x + 1) + (y - 4) * (y - 4);
}

double dfx(double x, double y)
{
	return 2 * (x + 1);
}

double dfy(double x, double y)
{
	return 2 * (y - 4);
}

vector<double> hlm(double dfx(double, double), double dfy(double, double), double lambda, double x, double y)
{
	vector<double> hlm = {0, 0};
	vector<vector<double>> H = {{0.5, 0}, {0, 0.5}};

	hlm.at(0) = (H.at(0).at(0) * dfx(x, y) + H.at(0).at(1) * dfy(x, y)) + lambda * dfx(x, y);
	hlm.at(1) = (H.at(1).at(0) * dfx(x, y) + H.at(1).at(1) * dfy(x, y)) + lambda * dfy(x, y);

	return hlm;
}

void levenberg(double f(double, double), double dfx(double, double), double dfy(double, double), double x, double y, double lambda, int numIt)
{
	vector<double> hlm1 = hlm(dfx, dfy, lambda, x, y);
	double next_x = x - hlm1.at(0);
	double next_y = y - hlm1.at(1);

	for (int i = 0; i < numIt; ++i)
	{
		cout << "x = " << next_x << endl;
		cout << "y = " << next_y << endl;
		cout << "f(x,y) = " << f(next_x, next_y) << endl;

		if (f(next_x, next_y) < f(x, y))
		{
			lambda *= 2;
			x = next_x;
			y = next_y;
		}

		else
			lambda /= 2;

		double next_x = x - hlm1.at(0);
		double next_y = y - hlm1.at(1);
	}
}

int main()
{
	return 0;
}