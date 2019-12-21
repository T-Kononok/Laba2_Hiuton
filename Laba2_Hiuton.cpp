#include "pch.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "Gaus.h"

using namespace std;

void Out(vector<double>);
void Out(vector<vector<double>>);
void Nev(vector<double>&, vector<double>);
void Jac(vector< vector<double>>&, vector<double>);
void Jac2(vector< vector<double>>&, vector<double>, double);
double fd1(vector<double>);
double fd2(vector<double>, vector<double>);

int main()
{
	int n = 2;
	vector<double> x;
	vector<double> x1;
	vector<double> deltx;
	vector<double> F;
	vector<int> p;
	vector<vector<double>> J;
	J.resize(n);
	for (int i = 0; i < n; i++)
	{
		F.push_back(0);
		x1.push_back(0);
		deltx.push_back(0);
		p.push_back(i);
		for (int j = 0; j < n; j++)
			J[i].push_back(0);
	}
	x.push_back(1);
	x.push_back(1);

	double e1 = 0.000000001, e2 = 0.000000001, d1 = 1, d2 = 1;
	int k = 1, NIT = 100;
	cout << "k				d1				d2" << endl;


	while ((d1 > e1 || d2 > e2) && k < NIT)
	{
		Nev(F, x);
		//Jac(J, x);
		Jac2(J, x, 0.01);

		if (!Gaus(p, deltx, J, F))
			return 0;

		for (int i = 0; i < n; i++)
			x1[i] = x[i] + deltx[i];

		d1 = fd1(x);
		d2 = fd2(x, x1);
		cout << setprecision(20) << fixed << k << "				" << d1 << "		" << d2 << endl;
		k++;
		x = x1;
		if (k >= NIT)
			cout << "ERROR: IER = 2" << endl;
	}

	Out(x);

}

void Out(vector<double> mas)
{
	for (int i = 0; i < mas.size(); i++)
		cout << mas[i] << " ";
	cout << endl << endl;
}

void Out(vector<vector<double>> mas)
{
	int n = mas[0].size();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << mas[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

void Nev(vector<double>& F, vector<double> x)
{
	F[0] = 1.5 * pow(x[0], 3) - pow(x[1], 2) - 1;
	F[1] = x[0] * pow(x[1], 3) - x[1] - 4;
}

void Jac(vector< vector<double>>& J, vector<double> x)
{
	J[0][0] = 4.5 * pow(x[0], 2);
	J[0][1] = -2 * x[1];
	J[1][0] = pow(x[1], 3);
	J[1][1] = 3 * x[0] * pow(x[1], 2) - 1;
}

void Jac2(vector< vector<double>>& J, vector<double> x, double M)
{
	J[0][0] = ((1.5 * pow(x[0] + M * x[0], 3) - pow(x[1], 2) - 1) - (1.5 * pow(x[0], 3) - pow(x[1], 2) - 1)) / (M * x[0]);
	J[0][1] = ((1.5 * pow(x[0], 3) - pow(x[1] + M * x[1], 2) - 1) - (1.5 * pow(x[0], 3) - pow(x[1], 2) - 1)) / (M * x[1]);
	J[1][0] = (((x[0] + M * x[0]) * pow(x[1], 3) - x[1] - 4) - (x[0] * pow(x[1], 3) - x[1] - 4)) / (M * x[0]);
	J[1][1] = ((x[0] * pow(x[1] + M * x[1], 3) - (x[1] + M * x[1]) - 4) - (x[0] * pow(x[1], 3) - x[1] - 4)) / (M * x[1]);
}

double fd1(vector<double> x)
{
	double x1 = 1.5 * pow(x[0], 3) - pow(x[1], 2) - 1;
	double x2 = x[0] * pow(x[1], 3) - x[1] - 4;

	if (abs(x1) > abs(x2))
		return abs(x1);
	else
		return abs(x2);
}

double fd2(vector<double> x, vector<double> x1)
{
	int n = x.size();
	double max = 0;
	for (int i = 0; i < n; i++)
	{
		if (abs(x1[i]) < 1)
		{
			if (abs(x1[i] - x[i]) > max)
				max = abs(x1[i] - x[i]);
		}
		else
		{
			if (abs((x1[i] - x[i]) / x1[i]) > max)
				max = abs((x1[i] - x[i])/ x1[i]);
		}
	}

	return max;
}

