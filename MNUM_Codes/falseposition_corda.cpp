// False Position / Corda

#include <iostream>
using namespace std;

#define ERROR 0.00001
#define MAXREP 500

void abs_error(double exato, double aproximado) { cout << (exato - aproximado) * 100 << " %" << endl; }
void rel_error(double exato, double aproximado) { cout << (exato - aproximado) / exato * 100; }

double f(double x) { return x * x * x * x + 2 * x * x * x - x - 1; }

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

int main()
{
    return 0;
}