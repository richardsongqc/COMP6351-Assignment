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
	int N = 40, i = 0;
	CVector u(N + 1), bu(N+1);
	double dblLamda = 0;

	std::ofstream ofile;
	ofile.open("prob3.csv");

	ofile << ",";

	for (i = 0; i < N; i++)
	{
		ofile << (double)i / (double)N << ",";
	}

	ofile << "1" << std::endl;

	// Prepare the vectors of the discretized boundary value problem
	double dbl = pow(N, 2);

	for (i = 0; i < N+1; i++)
	{
		u[i] = 0;// (double)i / (double)N;
	}

	double deltaLamda = 1;
	for (dblLamda = 0; dblLamda < 5; dblLamda += deltaLamda)
	{
		int k = 0;
		int nDim = N - 1;
		CVector l1(nDim);
		CVector l2(nDim);
		CVector l3(nDim);
		CVector b(nDim), deltaU(nDim), g(nDim);

		printf("%.15f) \n", dblLamda);

		try
		{
			bu = u;
			do
			{
				u[N] = 0;


				for (i = 0; i < nDim; i++)
				{
					double u0 = 0;
					double u1 = 0;
					double u2 = 0;

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

					b[i] = -dbl* (u0 - 2 * u1 + u2) - dblLamda*(exp(u0) + 10 * exp(u1) + exp(u2)) / ((double)12);

					double d0 = dbl + dblLamda*(exp(u0) / ((double)12));
					double d1 = -2 * dbl + 10 * dblLamda* (exp(u1) / ((double)12));
					double d2 = dbl + dblLamda* (exp(u2) / ((double)12));

					l2[i] = d1;

					if (i > 0)
					{
						l3[i - 1] = d2;
					}

					if (i + 1 < nDim)
					{
						l1[i + 1] = d0;
					}
				}

				// Apply LU-Decomp to solve it
				LUDecompTridiagnoal(l1, l2, l3, b, deltaU, g);

				//printf("\nk = %d\n", k);
				for (i = 0; i < nDim; i++)
				{
					u[i + 1] = deltaU[i] + u[i + 1];
					if (_finite(u[i + 1]) == 0)
					{
						throw MyException();
					}
					//printf("deltaU[%d]=%15.11f\tu[%d]=%12.11f\n", i, deltaU[i], i, u[i+1] );
					// Output the results to one excel file, later we can analyze this file.
					//ofile << std::setw(25) << u[i + 1] << ",";
				}

				if (deltaU < 0.00000000001)
				{
					printf("\nk = %d\n", k);
					break;
				}

				//printf("\n");

				k++;
			} while (k < 1000);
		}
		catch (MyException & e)
		{
			cout << "Caught an exception of type: " << typeid(e).name() << endl;
			cout << "Lamda = " << dblLamda << endl;
			dblLamda -= deltaLamda;
			deltaLamda /= (double)10;
			u = bu;
			if (deltaLamda < pow(0.1,15))
			{
				break;
			}
			continue;
		}

		ofile << std::setprecision(17) << dblLamda << ",";
		for (i = 0; i < N; i++)
		{
			ofile << std::setw(15) << u[i] << ",";
		}
		ofile << "0," << GetIntegral(u,N) << std::endl;


	}

	ofile.close();

	printf("\n");
}