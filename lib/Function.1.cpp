#include <vector>
#include <cmath>
#include <stdio.h>
#include "Function.1.hpp"

double func1::Function::MGF(std::vector<double> x, int l, int i)
{
	double sum = 0;
	int num = (int)(x.size() / l);
	for (int j = 0; j < num; j++) {
		sum += x[i - 1 + (j * l)];
	}
	return sum / num;
}

std::vector<std::vector<double> > func1::Function::MGF(std::vector<double> x, int n)
{
	std::vector<std::vector<double> > F;
	int m = (int)(n / 3);
	for (int l = 1; l <= m; l++)
	{
		std::vector<double> H;
		while (H.size() < n)
		{
			for (int i = 1; i <= l && H.size() < n; i++)
			{
				H.push_back(this->MGF(x, l, i));
			}
		}
		F.push_back(H);
	}
	return F;
}

std::vector<std::vector<double> > func1::Function::ArrIn(std::vector<std::vector<double> > F)
{
	int w = F.size();
	int h = F[0].size();
	std::vector<double> item(w, 0);
	std::vector<std::vector<double> > Arr(h, item);
	for (int i = 0; i < F.size(); i++)
	{
		for (int j = 0; j < F[i].size(); j++)
		{
			Arr[j][i] = F[i][j];
		}
	}
	return Arr;
}

std::vector<double> func1::Function::Differential(std::vector<double> x)
{
	std::vector<double> xd;
	for (int i = 1; i < x.size(); i++)
	{
		xd.push_back(x[i] - x[i - 1]);
	}
	return xd;
}

double func1::Function::SumAdd(std::vector<std::vector<double> > x, int l, int t)
{
	double sum = 0;
	for (int i = 0; i < t; i++)
	{
		sum += x[l][i + 1];
	}
	return sum;
}

std::vector<std::vector<double> > func1::Function::SumAdd(std::vector<std::vector<double> > x, double x1)
{
	std::vector<std::vector<double> > xa(x);
	for (int i = 0; i < xa.size(); i++)
	{
		xa[i][0] = x1;
	}
	for (int l = 0; l < xa.size(); l++)
	{
		for (int i = 1; i < xa[l].size(); i++)
		{
			xa[l][i] = x1 + this->SumAdd(x, l, i);
		}
	}
	return xa;
}

double func1::Function::CalcQ(std::vector<double> x)
{
	double xb = this->MGF(x, 1, 1);
	int n = x.size();
	double sum = 0;
	for (int t = 0; t < n; t++)
	{
		sum += (x[t] - xb) * (x[t] - xb);
	}
	return sum / n;
}

double func1::Function::CalcQ(std::vector<double> x, std::vector<double> xt)
{
	int n = xt.size();
	double sum = 0;
	for (int t = 0; t < n; t++)
	{
		sum += (x[t] - xt[t]) * (x[t] - xt[t]);
	}
	return sum / n;
}

double func1::Function::CalcS1(std::vector<double> x, std::vector<double> xk)
{
	double Qx = this->CalcQ(x);
	double Qk = this->CalcQ(x, xk);
	int n = x.size();
	double S1 = n * (1 - (Qk / Qx));
	return S1;
}

double func1::Function::CalcUV(std::vector<double> x)
{
	double sum = 0;
	double n = x.size();
	for (int i = 0; i < 0; i++)
	{
		sum += fabs(x[i]);
	}
	return sum / n;
}

int func1::Function::CalcT(std::vector<double> x, double u, int t)
{
	int T = 0;
	if (x[t] > u)
	{
		T = 0;
	}
	if (fabs(x[t]) <= u)
	{
		T = 1;
	}
	if (x[t] < -u)
	{
		T = 2;
	}
	return T;
}

std::vector<std::vector<double> > func1::Function::CalcNij(std::vector<double> x, std::vector<double> f)
{
	std::vector<double> x1 = this->Differential(x);
	std::vector<double> f1 = this->Differential(f);
	double U = this->CalcUV(x1);
	double V = this->CalcUV(f1);
	std::vector<double> Nj(3, 0);
	int i = 0, j = 0;
	std::vector<std::vector<double> > Nij(3, Nj);
	for (int t = 0; t < f1.size(); t++) {
		i = this->CalcT(x1, U, t);
		j = this->CalcT(f1, V, t);
		Nij[i][j]++;
	}
	return Nij;
}

double func1::Function::CalcS2(std::vector<double> x, std::vector<double> f) {
	std::vector<std::vector<double> > Nij = this->CalcNij(x, f);
	double R1 = 0, R2 = 0, R3 = 0, Ni = 0, Nj = 0, tem = 0;
	int n = x.size();
	for (int i = 0; i < Nij.size(); i++) {
		tem = 0;
		for (int j = 0; j < Nij[i].size(); j++) {
			if (Nij[i][j] > 0) {
				tem += Nij[i][j] * log(Nij[i][j]);
			}
		}
		R1 += tem;
	}
	for (int i = 0; i < Nij.size(); i++) {
		Ni = 0;
		for (int j = 0; j < Nij.size(); j++) {
			Ni += Nij[i][j];
		}
		if (Ni > 0) {
			R2 += Ni * log(Ni);
		}
	}
	for (int j = 0; j < Nij.size(); j++) {
		Nj = 0;
		for (int i = 0; i < Nij.size(); i++) {
			Nj += Nij[i][j];
		}
		if (Nj > 0) {
			R3 += Nj * log(Nj);
		}
	}
	return 2 * (R1 + (n - 1) * log(n - 1) - R2 - R3);
}

double func1::Function::CalcCSC(std::vector<double> x, std::vector<double> f) {
	double S1 = this->CalcS1(x, f);
	double S2 = this->CalcS2(x, f);
	return S1 + S2;
}

double func1::Function::CalcCSC(std::vector<double> x, std::vector<std::vector<double> > F) {
	std::vector<double> f = this->ComeBack(x, F);
	return this->CalcCSC(x, f);
}

std::vector<double> func1::Function::CalcCSCs(std::vector<double> x, std::vector<std::vector<double> > F) {
	std::vector<double> CSC;
	for (int i = 0; i < F.size(); i++) {
		std::vector<std::vector<double> > t;
		t.push_back(F[i]);
		CSC.push_back(this->CalcCSC(x, t));
	}
	return CSC;
}

std::vector<double> func1::Function::CalcCSCs(std::vector<double> x, std::vector<std::vector<std::vector<double> > > F) {
	std::vector<double> CSC;
	for (int i = 0; i < F.size(); i++) {
		CSC.push_back(this->CalcCSC(x, F[i]));
	}
	return CSC;
}

double func1::Function::CalcX2(std::vector<double> x, std::vector<double> f) {
	double X2 = 0;
	for (int i = 0; i < f.size(); i++) {
		if (f[i] != 0) {
			X2 += (x[i] - f[i]) * (x[i] - f[i]) / f[i];
		}
	}
	return X2;
}

std::vector<double> func1::Function::CalcX2(std::vector<double> x, std::vector<std::vector<double> > F) {
	std::vector<double> X2;
	for (int i = 0; i < F.size(); i++) {
		X2.push_back(this->CalcX2(x, F[i]));
	}
	return X2;
}

std::vector<int> func1::Function::RouSelect(std::vector<double> c, double x) {
	std::vector<int> div;
	for (int i = 0; i < c.size(); i++) {
		if (c[i] > x) {
			div.push_back(i);
		}
	}
	return div;
}

std::vector<std::vector<double> > func1::Function::RouSelect(std::vector<std::vector<std::vector<double> > > F, std::vector<std::vector<double> > CSC, double x) {
	std::vector<std::vector<double> > P;
	for (int i = 0; i < CSC.size(); i++) {
		std::vector<int> temCSC = this->RouSelect(CSC[i], x);
		for (int j = 0; j < temCSC.size(); j++) {
			P.push_back(F[i][temCSC[j]]);
		}
	}
	return P;
}

void func1::Function::TwoAddOne(std::vector<int> &two) {
	int n = two.size(), t = 0;
	two[n - 1]++;
	for (int i = n - 1; i >= 0; i--) {
		two[i] += t;
		if (two[i] > 1) {
			t = 1;
			two[i] = 0;
		}
		else {
			t = 0;
			break;
		}
	}
}

bool func1::Function::TwoIsFull(std::vector<int> two) {
	bool result = true;
	for (int i = 0; i < two.size(); i++) {
		if (two[i] == 0) {
			result = false;
			break;
		}
	}
	return result;
}

std::vector<std::vector<std::vector<double> > > func1::Function::Group(std::vector<std::vector<double> > P) {
	std::vector<std::vector<std::vector<double> > > Son;
	std::vector<int> two(P.size(), 0);
	while (!this->TwoIsFull(two)) {
		std::vector<std::vector<double> > tem;
		this->TwoAddOne(two);
		for (int i = 0; i < P.size(); i++) {
			if (two[i] == 0) {
				continue;
			}
			tem.push_back(P[i]);
		}
		Son.push_back(tem);
	}
	return Son;
}

std::vector<std::vector<double> > func1::Function::Mul(std::vector<std::vector<double> > arrA, std::vector<std::vector<double> > arrB) {
	int rowA = arrA.size();
	int colA = arrA[0].size();
	int rowB = arrB.size();
	int colB = arrB[0].size();
	std::vector<std::vector<double> > res;
	if (colA != rowB) {
		return res;
	}
	else {
		res.resize(rowA);
		for (int i = 0; i < rowA; ++i) {
			res[i].resize(colB);
		}
		for (int i = 0; i < rowA; ++i) {
			for (int j = 0; j < colB; ++j) {
				for (int k = 0; k < colA; ++k) {
					res[i][j] += arrA[i][k] * arrB[k][j];
				}
			}
		}
	}
	return res;
}

std::vector<std::vector<double> > func1::Function::Inv(std::vector<std::vector<double> > a)
{
	int n = a.size();
	std::vector<std::vector<double> > res(n);
	for (int i = 0; i < n; i++) {
		res[i].resize(n);
	}
	int *is = new int[n];
	int *js = new int[n];
	int i, j, k;
	double d, p;
	for (k = 0; k < n; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++)
			{
				p = fabs(a[i][j]);
				if (p > d)
				{
					d = p; is[k] = i; js[k] = j;
				}
			}
		if (0.0 == d)
		{
			delete[] is; delete[] js;
			return res;
		}
		if (is[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				p = a[k][j];
				a[k][j] = a[is[k]][j];
				a[is[k]][j] = p;
			}
		if (js[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				p = a[i][k];
				a[i][k] = a[i][js[k]];
				a[i][js[k]] = p;
			}
		a[k][k] = 1.0 / a[k][k];
		for (j = 0; j <= n - 1; j++)
			if (j != k)
			{
				a[k][j] *= a[k][k];
			}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
				for (j = 0; j <= n - 1; j++)
					if (j != k)
					{
						a[i][j] -= a[i][k] * a[k][j];
					}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
			{
				a[i][k] = -a[i][k] * a[k][k];
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				p = a[k][j];
				a[k][j] = a[js[k]][j];
				a[js[k]][j] = p;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				p = a[i][k];
				a[i][k] = a[i][is[k]];
				a[i][is[k]] = p;
			}
	}
	delete[] is; delete[] js;
	res = a;
	return res;
}

std::vector<std::vector<double> > Add(std::vector<std::vector<double> > arrA, std::vector<std::vector<double> > arrB)
{
	int rowA = arrA.size();
	int colA = arrA[0].size();
	int rowB = arrB.size();
	int colB = arrB[0].size();
	std::vector<std::vector<double> >  res;
	if ((colA != colB) || (rowA != rowB))
	{
		return res;
	}
	else
	{
		res.resize(rowA);
		for (int i = 0; i < rowA; ++i)
		{
			res[i].resize(colB);
		}
		for (int i = 0; i < rowA; ++i)
		{
			for (int j = 0; j < colB; ++j)
			{

				res[i][j] = arrA[i][j] + arrB[i][j];

			}
		}
	}
	return res;
}

std::vector<std::vector<double> > func1::Function::Xb(std::vector<double> x) {
	std::vector<std::vector<double> > result;
	std::vector<double> temp(x.size(), 1);
	result.push_back(temp);
	result.push_back(x);
	return ArrIn(result);
}

std::vector<std::vector<double> > func1::Function::Xb(std::vector<std::vector<double> > F) {
	std::vector<std::vector<double> > result;
	std::vector<double> temp(F[0].size(), 1);
	result.push_back(temp);
	for (int i = 0; i < F.size(); i++) {
		result.push_back(F[i]);
	}
	return ArrIn(result);
}

std::vector<std::vector<double> > func1::Function::ComeBackP(std::vector<double> y, std::vector<std::vector<double> > Xb) {
	std::vector<std::vector<double> > p;
	p.push_back(y);
	p = this->ArrIn(p);
	std::vector<std::vector<double> > XbT = this->ArrIn(Xb);
	std::vector<std::vector<double> > XbTMXb = this->Mul(XbT, Xb);
	std::vector<std::vector<double> > XbTMXbN = this->Inv(XbTMXb);
	std::vector<std::vector<double> > XbTMXbNMXbT = this->Mul(XbTMXbN, XbT);
	std::vector<std::vector<double> > res = this->Mul(XbTMXbNMXbT, p);
	return res;
}

std::vector<double> func1::Function::ComeBack(std::vector<double> x, std::vector<std::vector<double> > F) {
	std::vector<std::vector<double> > Xb = this->Xb(F);
	std::vector<std::vector<double> > a = this->ComeBackP(x, Xb);
	std::vector<double> result;
	long Frow = F.size();
	long Fcol = F[0].size();
	for (int i = 0; i < Fcol; i++) {
		double res = a[0][0];
		for (int j = 0; j < Frow; j++) {
			res += a[j + 1][0] * F[j][i];
		}
		result.push_back(res);
	}
	return result;
}

int func1::Function::FindMaxCSC(std::vector<double> CSC) {
	int res = 0;
	double maxCSC = CSC[0];
	for (int i = 1; i < CSC.size(); i++) {
		if (CSC[i] > maxCSC) {
			maxCSC = CSC[i];
			res = i;
		}
	}
	return res;
}

std::vector<double> func1::Function::Predict(std::vector<std::vector<double> > A, std::vector<std::vector<double> > Son, int q) {
	std::vector<double> result;
	long Frow = Son.size();
	long Fcol = Son[0].size();
	for (int i = 0; i < Fcol; i++) {
		double res = A[0][0];
		for (int j = 0; j < Frow; j++) {
			res += A[j + 1][0] * Son[j][i];
		}
		result.push_back(res);
	}
	return result;
}

std::vector<double> func1::Function::Predict(std::vector<double> x, int q) {
	int n = x.size();
	std::vector<std::vector<double> > F0 = this->MGF(x, n);
	std::vector<double> x1 = this->Differential(x);
	std::vector<std::vector<double> > F1 = this->MGF(x1, n);
	std::vector<double> x2 = this->Differential(x1);
	std::vector<std::vector<double> > F2 = this->MGF(x2, n);
	std::vector<std::vector<double> > F3 = this->SumAdd(F1, x[0]);
	// Predict
	std::vector<std::vector<double> > F0p = this->MGF(x, n + q);
	std::vector<double> x1p = this->Differential(x);
	std::vector<std::vector<double> > F1p = this->MGF(x1, n + q);
	std::vector<double> x2p = this->Differential(x1);
	std::vector<std::vector<double> > F2p = this->MGF(x2, n + q);
	std::vector<std::vector<double> > F3p = this->SumAdd(F1p, x[0]);
	// Predict
	std::vector<double> CSC0 = this->CalcCSCs(x, F0);
	std::vector<double> CSC1 = this->CalcCSCs(x, F1);
	std::vector<double> CSC2 = this->CalcCSCs(x, F2);
	std::vector<double> CSC3 = this->CalcCSCs(x, F3);
	std::vector<std::vector<std::vector<double> > > F;
	F.push_back(F0);
	F.push_back(F1);
	F.push_back(F2);
	F.push_back(F3);
	std::vector<std::vector<double> > CSC;
	CSC.push_back(CSC0);
	CSC.push_back(CSC1);
	CSC.push_back(CSC2);
	CSC.push_back(CSC3);
	double xx = 11.07;
	std::vector<std::vector<double> > P = this->RouSelect(F, CSC, xx);
	double xxs[27] = {40.113, 38.885, 37.652, 36.415, 35.172, 33.924, 32.761, 31.41, 30.144, 28.869, 27.587, 26.296, 24.996, 23.685, 22.262, 21.026, 19.675, 18.307, 16.919, 15.507, 14.067, 12.592, 11.07, 9.488, 7.815, 5.911, 3.841};
	for (int i = 0; i < 27; i++) {
		P = this->RouSelect(F, CSC, xxs[i]);
		if (P.size() > 6) {
			break;
		}
	}
	std::vector<std::vector<std::vector<double> > > Son = this->Group(P);
	std::vector<double> lastCSClist = this->CalcCSCs(x, Son);
	int lastCSCi = this->FindMaxCSC(lastCSClist);
	std::vector<std::vector<double> > lastSon = Son[lastCSCi];
	std::vector<std::vector<double> > lastXb = this->Xb(lastSon);
	std::vector<std::vector<double> > lastA = this->ComeBackP(x, lastXb);
	// Predict
	std::vector<std::vector<double> > newSon;
	for (int i = 0; i < lastSon.size(); i++) {
		bool isFind = false;
		for (int j = 0; j < F0.size(); j++) {
			bool hasDif = false;
			for (int k = 0; k < F0[j].size(); k++) {
				if (fabs(lastSon[i][k] - F0[j][k]) < 0.0001) {
					continue;
				}
				hasDif = true;
			}
			if (hasDif) {
				continue;
			}
			else {
				isFind = true;
				newSon.push_back(F0p[j]);
				break;
			}
		}
		if (isFind) {
			continue;
		}
		for (int j = 0; j < F1.size(); j++) {
			bool hasDif = false;
			for (int k = 0; k < F1[j].size(); k++) {
				if (fabs(lastSon[i][k] - F1[j][k]) < 0.0001) {
					continue;
				}
				hasDif = true;
			}
			if (hasDif) {
				continue;
			}
			else {
				isFind = true;
				newSon.push_back(F1p[j]);
				break;
			}
		}
		if (isFind) {
			continue;
		}
		for (int j = 0; j < F2.size(); j++) {
			bool hasDif = false;
			for (int k = 0; k < F2[j].size(); k++) {
				if (fabs(lastSon[i][k] - F2[j][k]) < 0.0001) {
					continue;
				}
				hasDif = true;
			}
			if (hasDif) {
				continue;
			}
			else {
				isFind = true;
				newSon.push_back(F2p[j]);
				break;
			}
		}
		if (isFind) {
			continue;
		}
		for (int j = 0; j < F3.size(); j++) {
			bool hasDif = false;
			for (int k = 0; k < F3[j].size(); k++) {
				if (fabs(lastSon[i][k] - F3[j][k]) < 0.0001) {
					continue;
				}
				hasDif = true;
			}
			if (hasDif) {
				continue;
			}
			else {
				isFind = true;
				newSon.push_back(F3p[j]);
				break;
			}
		}
		if (isFind) {
			continue;
		}
	}
	// Predict
	std::vector<double> result = this->Predict(lastA, newSon, q);
	return result;
}

std::vector<double> func1::Function::PredictE(std::vector<double> x, int q) {
	std::vector<double> nhz = this->Predict(x, q);
	std::vector<double> e;
	for (int i = 0; i < x.size(); i++) {
		e.push_back(x[i] - nhz[i]);
	}
	std::vector<double> ePre = this->Predict(e, q);
	std::vector<double> result;
	for (int i = x.size(); i < nhz.size(); i++) {
		result.push_back(nhz[i] + ePre[i]);
	}
	return result;
}
