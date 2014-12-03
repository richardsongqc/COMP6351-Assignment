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
	CVector u(N + 1), bu(N+1), pu(N+1);
	double t = 0;
	double deltaT = pow(0.1, 4);
	double dblLamda = 4;

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

	//for (i = 1; i < N+1; i++)
	//{
	//	u[i] = 0; //(double)i / (double)N;
	//}

	u[0] = 0;
	u[1] = 0.089475435;
	u[2] = 0.178607629;
	u[3] = 0.267450426;
	u[4] = 0.356054866;
	u[5] = 0.444469263;
	u[6] = 0.532739251;
	u[7] = 0.620907795;
	u[8] = 0.709015164;
	u[9] = 0.79709887;
	u[10] = 0.885193558;
	u[11] = 0.973330858;
	u[12] = 1.061539183;
	u[13] = 1.149843472;
	u[14] = 1.238264878;
	u[15] = 1.326820382;
	u[16] = 1.415522335;
	u[17] = 1.504377924;
	u[18] = 1.59338854;
	u[19] = 1.682549048;
	u[20] = 1.771846948;
	u[21] = 1.861261419;
	u[22] = 1.950762229;
	u[23] = 2.040308515;
	u[24] = 2.129847411;
	u[25] = 2.219312539;
	u[26] = 2.308622345;
	u[27] = 2.3976783;
	u[28] = 2.486362975;
	u[29] = 2.574538015;
	u[30] = 2.662042065;
	u[31] = 2.748688676;
	u[32] = 2.834264312;
	u[33] = 2.918526516;
	u[34] = 3.001202403;
	u[35] = 3.08198761;
	u[36] = 3.160545907;
	u[37] = 3.236509674;
	u[38] = 3.309481461;
	u[39] = 3.379036865;
	u[40] = 3.444728924;
	u[41] = 3.506094187;
	u[42] = 3.562660541;
	u[43] = 3.613956754;
	u[44] = 3.659523561;
	u[45] = 3.698925955;
	u[46] = 3.731766159;
	u[47] = 3.757696608;
	u[48] = 3.776432173;
	u[49] = 3.787760793;
	u[50] = 3.791551726;
	u[51] = 3.787760793;
	u[52] = 3.776432173;
	u[53] = 3.757696608;
	u[54] = 3.731766159;
	u[55] = 3.698925955;
	u[56] = 3.659523561;
	u[57] = 3.613956754;
	u[58] = 3.562660541;
	u[59] = 3.506094187;
	u[60] = 3.444728924;
	u[61] = 3.379036865;
	u[62] = 3.309481461;
	u[63] = 3.236509674;
	u[64] = 3.160545907;
	u[65] = 3.08198761;
	u[66] = 3.001202403;
	u[67] = 2.918526516;
	u[68] = 2.834264312;
	u[69] = 2.748688676;
	u[70] = 2.662042065;
	u[71] = 2.574538015;
	u[72] = 2.486362975;
	u[73] = 2.3976783;
	u[74] = 2.308622345;
	u[75] = 2.219312539;
	u[76] = 2.129847411;
	u[77] = 2.040308515;
	u[78] = 1.950762229;
	u[79] = 1.861261419;
	u[80] = 1.771846948;
	u[81] = 1.682549048;
	u[82] = 1.59338854;
	u[83] = 1.504377924;
	u[84] = 1.415522335;
	u[85] = 1.326820382;
	u[86] = 1.238264878;
	u[87] = 1.149843472;
	u[88] = 1.061539183;
	u[89] = 0.973330858;
	u[90] = 0.885193558;
	u[91] = 0.79709887;
	u[92] = 0.709015164;
	u[93] = 0.620907795;
	u[94] = 0.532739251;
	u[95] = 0.444469263;
	u[96] = 0.356054866;
	u[97] = 0.267450426;
	u[98] = 0.178607629;
	u[99] = 0.089475435;
	u[100] = 0;





	pu = u;

	for (t = 0.865; t < 0.9; t += deltaT)
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

					b[i] = deltaT*u0 / h2 - (1 + 2 * deltaT / h2)*u1 + deltaT*u2 / h2 + pu1 + deltaT*dblLamda*exp(u1);

					double d0 = -deltaT / h2;
					double d1 = 1 + 2 * deltaT / h2 - deltaT*dblLamda*exp(u1);

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
			} 
			while (k < 1000);
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
		ofile << std::setprecision(17) << t << ",";
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