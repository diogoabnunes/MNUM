// Sistema de Newton

#include <iostream>
using namespace std;

#define ERROR 0.00001

double f1(double x, double y) { return sqrt((x * (y + 5) - 1) / 2); }
double f2(double x, double y) { return sqrt(x + log10(x)); }

double df1x(double x, double y) { return (y + 5) / (pow(2, 3 / 2) * sqrt(x * (y + 5) - 1)); }
double df1y(double x, double y) { return x / (pow(2, 3 / 2) * sqrt(x * (y + 5) - 1)); }
double df2x(double x, double y) { return (1 / x + 1) / (2 * sqrt(log10(x) + x)); }
double df2y(double x, double y) { return 0; }

double J(double x, double y) { return df1x(x, y) * df2y(x, y) - df1y(x, y) * df2x(x, y); }

double h(double x, double y) { return (f1(x, y) * df2y(x, y) - f2(x, y) * df1y(x, y)) / J(x, y); }
double k(double x, double y) { return (df1x(x, y) * f2(x, y) - df2x(x, y) * f1(x, y)) / J(x, y); }

void sist_newton(double x, double y)
{
	double x0, y0;
	//do {
	for (int i = 1; i <= 2; i++) {
		x0 = x;
		y0 = y;
		x = x - h(x0, y0);
		y = y - k(x0, y0);
		//cout << "x: " << x << "\ny: " << y << "\n\n";
	} //while (abs(x0 - x) >= ERROR && abs(y0 - y) >= ERROR);
	cout << "Sistema newton:\nx: " << x << "\ny: " << y << endl;
}

int main()
{
    return 0;
}