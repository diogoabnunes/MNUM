// Teste1_MNUM.cpp | Diogo Nunes

#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;
#define ERROR 0.00001
#define MAXREP 500

void abs_error(double exato, double aproximado) { cout << (exato - aproximado) * 100 << " %" << endl; }
void rel_error(double exato, double aproximado) { cout << (exato - aproximado) / exato * 100; }

double f(double x) { return x * x * x * x + 2 * x * x * x - x - 1; }
double df(double x) { return 3 * x * x + 4 * x + 10; }

double g(double x) { return x * x - 4; }

double f1(double x, double y) { return sqrt((x * (y + 5) - 1) / 2); }
double f2(double x, double y) { return sqrt(x + log10(x)); }

double df1x(double x, double y) { return (y + 5) / (pow(2, 3 / 2) * sqrt(x * (y + 5) - 1)); }
double df1y(double x, double y) { return x / (pow(2, 3 / 2) * sqrt(x * (y + 5) - 1)); }
double df2x(double x, double y) { return (1 / x + 1) / (2 * sqrt(log10(x) + x)); }
double df2y(double x, double y) { return 0; }
double J(double x, double y) { return df1x(x, y) * df2y(x, y) - df1y(x, y) * df2x(x, y); }

double h(double x, double y) { return (f1(x, y) * df2y(x, y) - f2(x, y) * df1y(x, y)) / J(x, y); }
double k(double x, double y) { return (df1x(x, y) * f2(x, y) - df2x(x, y) * f1(x, y)) / J(x, y); }

void bisection(double a, double b)
{
	double m, count = 0;
	while (abs(a - b) >= ERROR && count < MAXREP)
	//for (int i = 1; i <=3; i++)
	{
		m = (a + b) / 2;
		//cout << a << "\t" << b << "\t" << m << "\t" << f(a) << "\t" << f(b) << "\t" << f(m) << endl;
		if (f(m) == 0.0) break;
		else if (f(a) * f(m) < 0) b = m;
		else a = m;
		count++;
	}
	cout << "Bisection: " << m << endl;
	//cout << "Abs Error: "; abs_error(1, m);
	//cout << "Rel Error: "; rel_error(1, m);
}

void falsep_corda(double a, double b)
{
	double w, count = 0;
	while (abs(a - b) >= ERROR && count < MAXREP)
	{
		w = (a * f(b) - b * f(a)) / (f(b) - f(a));
		if (f(w) == 0.0) break;
		else if (f(a) * f(w) < 0) b = w;
		else a = w;
		count++;
	}
	cout << "Falsep: " << w << endl;
	cout << "Abs Error: "; abs_error(1, w);
	cout << "Rel Error: "; rel_error(1, w);
}

void newton_tang(double guess)
{
	double x = guess, g, x0;
	do
	//for (int i = 1; i <= 3; i++)
	{
		cout << "xn: " << x << endl;
		x0 = x;
		g = f(x) / df(x);
		x = x - g;
	} while (abs(x0 - x) >= ERROR);
	// -> certo, mas exige que sucessao dos passos em cada iteracao seja monotona e decrescente
	
	cout << "Newton: " << x << endl;
	cout << "Abs Error: "; abs_error(1, x);
	cout << "Rel Error: "; rel_error(1, x);
}

void picardpeano(double guess)
{
	double x = guess, x0;
	do {
		x0 = x;
		x = g(x0);
	} while (abs(x0 - x) >= ERROR);
	cout << "Picardpeano: " << x << endl;
	cout << "Abs Error: "; abs_error(1, x);
	cout << "Rel Error: "; rel_error(1, x);
}

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

void sist_picardpeano(double x, double y)
{
	double x0, y0;
	do {
		x0 = x;
		y0 = y;
		x = f1(x0, y0);
		y = f2(x0, y0);
	} while (abs(x0 - x) >= ERROR && abs(y0 - y) >= ERROR);
	cout << "Sistema Picardpeano:\nx: " << x << "\ny: " << y << endl;
}

int main()
{
	
	return 0;
}