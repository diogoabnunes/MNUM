// Euler XYZ

#include <iostream>
#include <iomanip>
using namespace std;

double dz(double t, double z)
{
	return 2 + pow(t, 2) + t * z;
}

void euler_xyz(double t0, double y0, double z0, double h0, double num_iter)
{
	double t = t0, y = y0, z = z0, h = h0, deltaz;
	for (int i = 0; i < num_iter; i++)
	{
		cout << "Iteracao " << i << endl;
		cout << "t: " << t << endl;
		cout << "y: " << y << endl;
		deltaz = dz(t, z);
		t += h;
		y += h * z;
		z += h * deltaz;
	}
}

int main()
{
    return 0;
}