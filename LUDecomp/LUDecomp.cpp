#include <stdio.h>
#include "Matrix.h"
#include "Vector.h"

void LUDecomp(CMatrix a, CVector b, CMatrix& l, CMatrix& u, CVector& x)
{
	CVector  z(b.GetSize());
	int i = 0, j = 0, k = 0, n = b.GetSize();

	// computations of L and U matrices
	for (i = 0; i < n; i++)
	{
		l[i][0] = a[i][0];
	}

	for (j = 1; j < n; j++)
	{
		u[0][j] = a[0][j] / l[0][0];
	}

	for (i = 0; i < n; i++)
	{
		u[i][i] = 1;
	}

	for (i = 1; i < n; i++)
	{
		for (j = 1; j < n; j++)
		{
			if (i >= j)
			{
				l[i][j] = a[i][j];
				for (k = 0; k <= j - 1; k++)
				{
					l[i][j] -= l[i][k] * u[k][j];
				}
			}
			else
			{
				u[i][j] = a[i][j];
				for (k = 0; k <= i - 1; k++)
				{
					u[i][j] -= l[i][k] * u[k][j];
				}
				u[i][j] /= l[i][i];
			}
		}
	}

	// solve lz=b by forward substitution
	z[0] = b[0] / l[0][0];
	for (i = 1; i < n; i++)
	{
		z[i] = b[i];
		for (j = 0; j <= i - 1; j++)
		{
			z[i] -= l[i][j] * z[j];
		}
		z[i] /= l[i][i];
	}

	//solve Ux=z by backward substitution 
	//x[n] = z[n];
	x = z;
	for (i = n - 1; i >= 0; i--)
	{
		x[i] = z[i];
		for (j = i + 1; j < n; j++)
		{
			x[i] -= u[i][j] * x[j];
		}
	}
}

// Tridiagnonal Systems
void LUDecompTridiagnoal(CMatrix a, CVector b, CMatrix& l, CMatrix& u, CVector& x)
{
	int i = 0, j = 0, k = 0, n = b.GetSize();
	CVector l1(b.GetSize());		//	vector a
	CVector l2(b.GetSize());		//	vector b
	CVector l3(b.GetSize());		//	vector c

	// Get the elements and put them into the vector a, b, c
	for (i = 0; i < n; i++)
	{
		if (i > 0)
		{
			l3[i-1] = a[i][i-1];
		}

		if (i + 1 < n)
		{
			l1[i + 1] = a[i + 1][i];
		}

		l2[i] = a[i][i];
		
	}

	CVector alpha(b.GetSize());
	CVector beta(b.GetSize());
	CVector gama(b.GetSize());
	CVector g(b.GetSize());

	// Apply the algorithm of Gauss elimination ( Page 37 )
	beta[0] = l2[0];
	g[0]    = b[0];

	for (i = 1; i < n; i++)
	{
		gama[i] = l1[i] / beta[i-1];
		beta[i] = l2[i] - gama[i] * l3[i-1];
		g[i] = b[i] - gama[i]*g[i-1];
	}

	// computation of L matrix and U matrix
	for (i = 0; i < n; i++)
	{
		l[i][i] = 1;
		u[i][i] = beta[i];

		if (i + 1 < n)
		{
			l[i + 1][i] = gama[i + 1];
			u[i][i + 1] = l3[i];
		}
	}

	// backsubstitution algorithm ( Page 38 )
	x[n - 1] = g[n - 1] / beta[n - 1];
	for (i = n - 2; i >= 0; i-- )
	{
		x[i] = (g[i] - l3[i] * x[i + 1]) / beta[i];
	}

	return;
}


void main()
{
	int n = 4;
	CMatrix a(n, n);
	//a[0][0] = 4;     a[0][1] = 2;    a[0][2] = 1;
	//a[1][0] = 2;     a[1][1] = 5;    a[1][2] = -2;
	//a[2][0] = 1;     a[2][1] = -2;   a[2][2] = 7;
	a[0][0] = 3;     a[0][1] = 1;    a[0][2] = 0;	a[0][3] = 0;
	a[1][0] = 1;     a[1][1] = 3;    a[1][2] = 1;	a[1][3] = 0;
	a[2][0] = 0;     a[2][1] = 1;    a[2][2] = 3;	a[2][3] = 1;
	a[3][0] = 0;     a[3][1] = 0;    a[3][2] = 1;	a[3][3] = 3;

	CMatrix l(n, n);
	CMatrix u(n, n);

	CVector b(n), x(b.GetSize());
	b[0] = 4;        b[1] = 5;       b[2] = 5;		b[3] = 4;
	//b[0] = 3;        b[1] = 4;       b[2] = 5;

	int i = 0, j = 0, k = 0;
	printf("\nEnter the size of the coeficient matrix : %d\n", n);

	printf("Enter the elements rowwise \n");

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%6.2f\t", a[i][j]);
		}
		printf("\n");
	}

	printf("\nEnter the right hand vector\n");

	for (i = 0; i < n; i++)
	{
		printf("%6.2f\t", b[i]);
	}
	printf("\n");

	LUDecompTridiagnoal(a, b, l, u, x);

	printf("\nThe lower triangular matrix L\n");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j <= i; j++)
		{
			printf("%6.2f\t", l[i][j]);
		}
		printf("\n");
	}

	printf("\nThe upper triangular matrix U\n");

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < i; j++)
		{
			printf("      \t" );
		}

		for (j = i; j < n; j++)
		{
			printf("%6.2f\t", u[i][j]);
		}
		
		printf("\n");
	}

	printf("\n********************************\nThe solution is \n");
	for (i = 0; i < n; i++)
	{
		printf("%6.2f ", x[i]);
	}

	printf("\n");
}