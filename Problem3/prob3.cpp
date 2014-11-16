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

	u[0]  = 0.0000000000000;
	u[1]  = 0.2823891844000;
	u[2]  = 0.5616385780200;
	u[3]  = 0.8367352184100;
	u[4]  = 1.1063627597000;
	u[5]  = 1.3047245315000;
	u[6]  = 1.4981548973000;
	u[7]  = 1.6856021186000;
	u[8]  = 1.8658348095000;
	u[9]  = 2.0001139772000;
	u[10] = 2.1284244625000;
	u[11] = 2.2499505967000;
	u[12] = 2.3638180656000;
	u[13] = 2.4820436474000;
	u[14] = 2.5879537854000;
	u[15] = 2.6801791150000;
	u[16] = 2.7574051180000;
	u[17] = 2.8110651920000;
	u[18] = 2.8519926112000;
	u[19] = 2.8796584800000;
	u[20] = 2.8936924388000;
	u[21] = 2.8938958517000;
	u[22] = 2.8802500390000;
	u[23] = 2.8529420916000;
	u[24] = 2.8123377648000;
	u[25] = 2.7585628902000;
	u[26] = 2.6925236742000;
	u[27] = 2.6150004774000;
	u[28] = 2.5268465560000;
	u[29] = 2.4103991867000;
	u[30] = 2.2818895009000;
	u[31] = 2.1427664169000;
	u[32] = 1.9944047168000;
	u[33] = 1.7777560044000;
	u[34] = 1.5489690335000;
	u[35] = 1.3105196672000;
	u[36] = 1.0644562450000;
	u[37] = 0.8045977215000;
	u[38] = 0.5398528351600;
	u[39] = 0.2713590271900;
	u[40] = 0.0000000000000;

	double deltaLamda = 0.0001;
	for (dblLamda = 2; dblLamda < 5; dblLamda += deltaLamda)
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
			} 
			while (k < 1000);
		}
		catch (MyException & e)
		{
			//cout << "*****************************************************************" << endl;
			//cout << "\t!!!Lamda = " << dblLamda << "!!!" << endl;
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

		printf("\nLamda = %.15f \t delta = %.15f\n", dblLamda, deltaLamda);
		ofile << std::setprecision(17) << dblLamda << ",";
		for (i = 0; i < N; i++)
		{
			ofile << std::setw(15) << u[i] << ",";
		}
		ofile << "0," << GetIntegral(u,N) << std::endl;

		bu = u;
	}

	ofile.close();

	printf("\n");
}