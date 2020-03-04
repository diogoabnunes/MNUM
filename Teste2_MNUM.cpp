// Teste2_MNUM.cpp

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

// Gauss

typedef vector<vector<double>> matrix;

void printMatrix(matrix& m)
{
	for (unsigned int x = 0; x < m.size(); x++)
	{
		for (unsigned int y = 0; y < m[x].size(); y++)
		{
			cout << setw(10) << m[x][y];
		}
		cout << endl;
	}
	cout << endl;
}

void rowOp(matrix& m, unsigned int a, unsigned int b, double k)
{
	for (unsigned int i = 0; i < m[0].size(); i++)
	{
		m[a][i] -= k * m[b][i];
	}
}

vector<double> gauss(matrix& m)
{
	vector<double> res = { 0, 0, 0 };
	for (unsigned int i = 0; i < m.size(); i++)
	{
		rowOp(m, i, i, 1 - 1 / m[i][i]);
		for (unsigned int j = i + 1; j < m.size(); j++)
		{
			if (i != j) rowOp(m, j, i, m[j][i]);
		}
	}
	res.at(2) = m[2][3] / m[2][2];
	res.at(1) = (m[1][3] - m[1][2] * res.at(2)) / m[1][1];
	res.at(0) = (m[0][3] - m[0][2] * res.at(2) - m[0][1] * res.at(1)) / m[0][0];
	return res;
}

void printResult(vector<double> res) { cout << "(x y z) = (" << res.at(0) << " " << res.at(1) << " " << res.at(2) << ")" << endl; }

void printResidue(matrix& m, vector<double> res)
{
	cout << "1st line: " << m[0][3] - (m[0][0] * res.at(0) + m[0][1] * res.at(1) + m[0][2] * res.at(2)) << endl;
	cout << "2nd line: " << m[1][3] - (m[1][0] * res.at(0) + m[1][1] * res.at(1) + m[1][2] * res.at(2)) << endl;
	cout << "3rd line: " << m[2][3] - (m[2][0] * res.at(0) + m[2][1] * res.at(1) + m[2][2] * res.at(2)) << endl;
}

// Khaletsky

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

// Gauss Iterative

/*
Converge: Sim, porque em cada linha da matriz A, o módulo do elemento
da diagonal principal é superior ao módulo da soma dos restantes elementos da linha.
*/

void gauss_jacobi(matrix& m, double x1, double x2, double x3, double x4, unsigned int num_iter)
{
	double x1n, x2n, x3n, x4n;
	cout << 0 << ":\t" << x1 << "\t" << x2 << "\t" << x3 << "\t" << x4 << endl;
	for (unsigned int i = 1; i <= num_iter; i++)
	{
		x1n = (m[0][4] - (m[0][1] * x2 + m[0][2] * x3 + m[0][3] * x4)) / m[0][0];
		x2n = (m[1][4] - (m[1][0] * x1 + m[1][2] * x3 + m[1][3] * x4)) / m[1][1];
		x3n = (m[2][4] - (m[2][0] * x1 + m[2][1] * x2 + m[2][3] * x4)) / m[2][2];
		x4n = (m[3][4] - (m[3][0] * x1 + m[3][1] * x2 + m[3][2] * x3)) / m[3][3];
		x1 = x1n; x2 = x2n; x3 = x3n; x4 = x4n;
		cout << i << ":\t" << x1 << "\t" << x2 << "\t" << x3 << "\t" << x4 << endl;
	}
}

void gauss_seidel(matrix& m, double x1, double x2, double x3, double x4, unsigned int num_iter)
{
	cout << 0 << ":\t" << x1 << "\t" << x2 << "\t" << x3 << "\t" << x4 << endl;
	for (unsigned int i = 1; i <= num_iter; i++)
	{
		x1 = (m[0][4] - (m[0][1] * x2 + m[0][2] * x3 + m[0][3] * x4)) / m[0][0];
		x2 = (m[1][4] - (m[1][0] * x1 + m[1][2] * x3 + m[1][3] * x4)) / m[1][1];
		x3 = (m[2][4] - (m[2][0] * x1 + m[2][1] * x2 + m[2][3] * x4)) / m[2][2];
		x4 = (m[3][4] - (m[3][0] * x1 + m[3][1] * x2 + m[3][2] * x3)) / m[3][3];
		cout << i << ":\t" << x1 << "\t" << x2 << "\t" << x3 << "\t" << x4 << endl;
	}
}

// Quadratura

double f(double x) { return 2; }

double trapezio(double a, double b, double h)
{
	double result = 0;
	double n = abs(b - a) / h;
	for (unsigned int i = 1; i < n; i++)
		result += f(a + i * h) * 2;
	result += f(a) + f(b);
	result *= h / 2;
	return result;
}

double simpson(double a, double b, double h)
{
	double result = 0;
	double n = abs(b - a) / h;
	for (unsigned int i = 1; i < n; i += 2)
		result += f(a + i * h) * 4;
	for (unsigned int i = 2; i < n; i += 2)
		result += f(a + i * h) * 2;
	result += f(a) + f(b);
	result *= h / 3;
	return result;
}

double det_QC(double s, double sl, double sll) { return (sl - s) / (sll - sl); }

double erro(double sl, double sll) { return (sll - sl) / 15; }

// Cubatura

double f(double x, double y) { return 2; }

double trapezio_c(double a, double A, double b, double B)
{
	double result = 0;
	for (unsigned int x = 0; x < 3; x++)
	{
		for (unsigned int y = 0; y < 3; y++)
		{
			if ((x == 0 || x == 2) && (y == 0 || y == 2)) result += f(x, y);
			else if (x == 0 || x == 2 || y == 0 || y == 2) result += f(x, y) * 2;
			else result += f(x, y) * 4;
		}
	}
	result /= 4;
	return result;
}

double simpson_c(double a, double A, double b, double B)
{
	double hx = abs(A - a) / 2, hy = abs(B - b) / 2;
	double result = 0;
	for (unsigned int y = 0; y < 3; y++)
	{
		double temp = 0;
		for (unsigned int x = 0; x < 3; x++)
		{
			if (x == 1) temp += f(x, y) * 4;
			else temp += f(x, y);
		}
		if (y == 1) result += temp * 4;
		else result += temp;
	}
	result *= (hx * hy) / 9;
	return result;
}

// EDO

double df(double x, double y) { return x*y; }

void euler(double x0, double y0, double h, unsigned int num_iter)
{
	double xn = x0, yn = y0;
	cout << 0 << setw(16) << xn << setw(16) << yn << setw(16) << df(xn, yn) << endl;
	for (unsigned int i = 1; i < num_iter; i++)
	{
		yn += h * df(xn, yn);
		xn += h;
		cout << i << setw(16) << xn << setw(16) << yn << setw(16) << df(xn, yn) << endl;
	}
}

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

void rk4_xyz(double t0, double y0, double z0, double h0, double num_iter)
{
	double t = t0, y = y0, z = z0, h = h0;
	double z1, z2, z3, z4, dz1, dz2, dz3, dz4;
	for (int i = 0; i < num_iter; i++)
	{
		cout << "Iteracao " << i << endl;
		cout << "t: " << t << endl;
		cout << "y: " << y << endl;
		z1 = (h * z);
		dz1 = h * dz(t, z);
		z2 = (h * (z + dz1 / 2));
		dz2 = (h * (dz(t + h / 2, z + dz1 / 2)));
		z3 = (h * (z + dz2 / 2));
		dz3 = (h * (dz(t + h / 2, z + dz2 / 2)));
		z4 = (h * (z + dz3));
		dz4 = (h * (dz(t + h, z + dz3)));
		y += z1 / 6 + z2 / 3 + z3 / 3 + z4 / 6;
		z += dz1 / 6 + dz2 / 3 + dz3 / 3 + dz4 / 6;
		t += h;
	}
}

void rk2(double x0, double y0, double h, unsigned int num_iter)
{
	double xn = x0, yn = y0;
	cout << 0 << setw(16) << xn << setw(16) << yn << setw(16) << df(xn, yn) << endl;
	for (unsigned int i = 1; i < num_iter; i++)
	{
		yn += h * df(xn + h / 2, yn + h / 2 * df(xn, yn));
		xn += h;
		cout << i << setw(16) << xn << setw(16) << yn << setw(16) << df(xn, yn) << endl;
	}
}

void rk4(double x0, double y0, double h, unsigned int num_iter)
{
	double xn = x0, yn = y0;
	cout << 0 << setw(16) << xn << setw(16) << yn << setw(16) << df(xn, yn) << endl;
	for (unsigned int i = 1; i < num_iter; i++)
	{
		double dy1 = h * df(xn, yn);
		double dy2 = h * df(xn + h / 2, yn + dy1 / 2);
		double dy3 = h * df(xn + h / 2, yn + dy2 / 2);
		double dy4 = h * df(xn + h, yn + dy3);
		yn += dy1 / 6 + dy2 / 3 + dy3 / 3 + dy4 / 6;
		xn += h;
		cout << i << setw(16) << xn << setw(16) << yn << setw(16) << df(xn, yn) << endl;
	}
}

/*
Estabilidade externa gauss: Maxima
A: matrix([,,],[,,],[,,]);
db = matrix([0.1], [0.1], [0.1]);
x: matrix([],[],[]);
dA: matrix([0.1, 0.1, 0.1], [0.1, 0.1, 0.1], [0.1, 0.1, 0.1]);
dx: invert(A).(db-dA.x);
*/

int main()
{
	

	return 0;
}
