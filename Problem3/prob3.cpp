#include <stdio.h>
#include "..\\cmn\\ludcom.h"
#include <fstream>
#include <iomanip>

void main()
{
	int N = 100, i = 0;
	CVector u(N + 1);
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

	for (dblLamda = 0; dblLamda < 40; dblLamda += 1)
	{
		int k = 0;
		int nDim = N - 1;
		CVector l1(nDim);
		CVector l2(nDim);
		CVector l3(nDim);
		CVector b(nDim), deltaU(nDim), g(nDim);

		printf("%.3f) \n", dblLamda);

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
					u1 = u[i+ 1];
					u2 = u[i + 2];
				}
				else if (i + 1 == nDim)
				{
					u0 = u[nDim-1];  
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
				u[i+1] = deltaU[i] + u[i+1];
				if (_finite(u[i + 1]) == 0)
				{
					throw;
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
		} 
		while (k < 1000  );

		ofile << std::setw(25) << dblLamda << ",";
		for (i = 0; i < N; i++ )
		{
			ofile << std::setw(25) << u[i] << ",";
		} 
		ofile << "0" << std::endl;
	}

	ofile.close();

	printf("\n");
}