// Khaletsky

#include <iostream>
#include <iomanip>
using namespace std;

double determinant(double mat[3][3])
{
	return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
		- mat[0][1] * (mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2])
		+ mat[0][2] * (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]);
}

void khaletsky(double coef[3][4])
{
	double d[3][3] = { {coef[0][0], coef[0][1], coef[0][2]},
						{coef[1][0], coef[1][1], coef[1][2]},
						{coef[2][0], coef[2][1], coef[2][2]} }; double D = determinant(d);
	double d1[3][3] = { {coef[0][3], coef[0][1], coef[0][2]},
						{coef[1][3], coef[1][1], coef[1][2]},
						{coef[2][3], coef[2][1], coef[2][2]} }; double D1 = determinant(d1);
	double d2[3][3] = { {coef[0][0], coef[0][3], coef[0][2]},
						{coef[1][0], coef[1][3], coef[1][2]},
						{coef[2][0], coef[2][3], coef[2][2]} }; double D2 = determinant(d2);
	double d3[3][3] = { {coef[0][0], coef[0][1], coef[0][3]},
						{coef[1][0], coef[1][1], coef[1][3]},
						{coef[2][0], coef[2][1], coef[2][3]} }; double D3 = determinant(d3);

	if (D != 0)
	{
		double x = D1 / D; cout << "x = " << x << endl;
		double y = D2 / D; cout << "y = " << y << endl;
		double z = D3 / D; cout << "z = " << z << endl;
	}
	else
	{
		if (D1 == 0 && D2 == 0 && D3 == 0) cout << "Infinite solutions" << endl;
		else if (D1 != 0 || D2 != 0 || D3 != 0) cout << "No solutions" << endl;
	}
}

int main()
{
    return 0;
}