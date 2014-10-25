#pragma once

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <vector>


class CVector 
{
private:
	// size of vector
	int m_maxSize;

	// data of vector
	std::vector<double> m_vec;		

public:
	// Constructor, zero-length  vector
	CVector();

	// Constructor, vector of length MaxSize
	CVector(int MaxSize);

	// Copy constructor
	CVector(const CVector& vec);							
	~CVector();												

	// Assignment operator 
	CVector& operator =(const CVector& vec);				
	                   
	bool operator==(double dbl) const;
	bool operator<(double dbl) const;

	// access to the elements of the CMatrix  
	double& operator [](int nPos);							
	const double& operator [](int nPos) const;

	// validate whether the size of vector is legal
	int CheckSize(int nSize);

	// get Rows and Columns of CMatrix
	int GetSize() const;
};



	
	
#endif //__VECTOR_H__