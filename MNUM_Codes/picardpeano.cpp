// Picard-Peano

#include <iostream>
using namespace std;

#define ERROR 0.00001

void abs_error(double exato, double aproximado) { cout << (exato - aproximado) * 100 << " %" << endl; }
void rel_error(double exato, double aproximado) { cout << (exato - aproximado) / exato * 100; }

double g(double x) { return x * x - 4; }

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

int main()
{
    return 0;
}