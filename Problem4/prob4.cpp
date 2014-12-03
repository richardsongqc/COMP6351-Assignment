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
	CVector u(N + 1), bu(N + 1), pu(N + 1);
	double t = 0;
	double deltaT = pow(0.1, 2);
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

	pu = u;
	int nIdx = 0;

	for (t = 0; t < 5; t += deltaT, nIdx++)
	{
		
		int k = 0;
		int nDim = N - 1;
		CVector l1(nDim);
		CVector l2(nDim);
		CVector l3(nDim);
		CVector b(nDim), deltaU(nDim), g(nDim);

		printf("%.15f) \n", t);

		try
		{

			do
			{
				u[N] = 0;

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
						pu1 = pu[i + 1];
					}
					else if (i + 2 <= nDim)
					{
						u0 = u[i];
						u1 = u[i + 1];
						u2 = u[i + 2];
						pu1 = pu[i + 1];
					}
					else if (i + 1 == nDim)
					{
						u0 = u[nDim - 1];
						u1 = u[nDim];
						u2 = u[N];
						pu1 = pu[nDim];
					}

					b[i] = deltaT*u0 / h2 - (1 + 2 * deltaT / h2)*u1 + deltaT*u2 / h2 + pu1 + deltaT*dblLamda*u1*(1-u1);

					double d0 = -deltaT / h2;
					double d1 = 1 + 2 * deltaT / h2 - deltaT*dblLamda*(1-2*u1);

					l2[i] = d1;

					if (i > 0)
					{
						l3[i - 1] = d0;
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
			cout << "*****************************************************************" << endl;
			cout << "\t!!!Lamda = " << dblLamda << "!!!" << endl;
			//dblLamda -= deltaLamda;
			//deltaLamda /= (double)10;
			//u = bu;
			//if (deltaLamda < pow(0.1,15))
			//{
			//	printf("\n============================================================\nLamda = %.15f\n", dblLamda);
			//	break;
			//}
			continue;
		}

		pu = u;

		printf("\tT = %.15f \t deltaT = %.15f\n", t, deltaT);

		if (nIdx % 2 == 0)
		{
			ofile << std::setprecision(17) << t << ",";
			for (i = 0; i < N; i++)
			{
				ofile << std::setw(15) << u[i] << ",";
			}
			ofile << "0," << GetIntegral(u, N) << std::endl;
		}

		bu = u;
	}

	ofile.close();

	printf("\n");
}