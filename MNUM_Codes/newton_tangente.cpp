// Newton / Tangente

#include <iostream>
using namespace std;

#define ERROR 0.00001

void abs_error(double exato, double aproximado) { cout << (exato - aproximado) * 100 << " %" << endl; }
void rel_error(double exato, double aproximado) { cout << (exato - aproximado) / exato * 100; }

double f(double x) { return x * x * x * x + 2 * x * x * x - x - 1; }
double df(double x) { return 3 * x * x + 4 * x + 10; }

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

int main()
{
    return 0;
}