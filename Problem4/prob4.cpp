#include <stdio.h>
#include "..\\cmn\\ludcom.h"
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

class MyException
{
};


double GetIntegral(CVector v, int N)
{
	double dbl = 0;
	double dblMeshSize = (double)1.0 / (double)N;
	for (int i = 0; i < v.GetSize(); i++)
	{
		dbl += v[i] * dblMeshSize;
	}

	return dbl;

}

void main()
{
	int N = 100, i = 0;
	double h = ((double)1) / ((double)N);
	double h2 = pow(h, 2);
	CVector u(N + 1), bu(N + 1) ;
	double t = 0;
	double deltaT = pow(0.1, 4);
	double dblLamda = 2;

	std::ofstream ofile;
	ofile.open("prob4.csv");

	ofile << ",";

	for (i = 0; i < N; i++)
	{
		ofile << (double)i / (double)N << ",";
	}

	ofile << "1" << std::endl;

	// Prepare the vectors of the discretized boundary value problem
	double dbl = pow(N, 2);

	for (i = 1; i < N + 1; i++)
	{
		u[i] = 1; //(double)i / (double)N;
	}

	for (t = 0; t < 1; t += deltaT)
	{
		int k = 0;
		int nDim = N - 1;
		CVector l1(nDim);
		CVector l2(nDim);
		CVector l3(nDim);
		CVector b(nDim), deltaU(nDim), g(nDim);

		printf("%.15f) \n", t);

		u[N] = 0;

		try
		{
			for (i = 0; i < nDim; i++)
			{
				double u0 = 0;
				double u1 = 0;
				double u2 = 0;
				double pu1 = 0;

				if (i == 0)
				{
					u0 = u[i];
					u1 = u[i + 1];
					u2 = u[i + 2];
				}
				else if (i + 2 <= nDim)
				{
					u0 = u[i];
					u1 = u[i + 1];
					u2 = u[i + 2];
				}
				else if (i + 1 == nDim)
				{
					u0 = u[nDim - 1];
					u1 = u[nDim];
					u2 = u[N];
				}

				bu[i+1] = u1 + deltaT *(u0 - 2 * u1 + u2) / h2 + deltaT*dblLamda*u1*(1-u1);
				if (_finite(bu[i]) == 0)
				{
					throw MyException();
				}

			}
		}
		catch (MyException & e)
		{
			cout << "*****************************************************************" << endl;
			break;
		}

		u = bu;

		printf("\tT = %.15f \t deltaT = %.15f\n", t, deltaT);
		ofile << std::setprecision(17) << t << ",";
		for (i = 0; i < N; i++)
		{
			ofile << std::setw(15) << u[i] << ",";
		}
		ofile << "0," << GetIntegral(u, N) << std::endl;
	}

	ofile.close();

	printf("\n");
}