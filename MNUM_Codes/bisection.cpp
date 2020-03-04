// Bisection

#include <iostream>
using namespace std;

#define ERROR 0.00001

void abs_error(double exato, double aproximado) { cout << (exato - aproximado) * 100 << " %" << endl; }
void rel_error(double exato, double aproximado) { cout << (exato - aproximado) / exato * 100; }

double f(double x) { return x * x * x * x + 2 * x * x * x - x - 1; }

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

int main()
{
    return 0;
}