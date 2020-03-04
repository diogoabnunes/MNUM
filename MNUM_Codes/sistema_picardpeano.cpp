// Sistema Picard-Peano

#include <iostream>
using namespace std;

#define ERROR 0.00001

double f1(double x, double y) { return sqrt((x * (y + 5) - 1) / 2); }
double f2(double x, double y) { return sqrt(x + log10(x)); }

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