// Gauss

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

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

int main()
{
    return 0;
}