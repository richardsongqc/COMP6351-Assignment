#pragma once

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>

typedef std::vector<std::vector<double> >    MATRIXELEMENTTYPE;

class CMatrix
{
private:
	int m_row;												// rows of CMatrix
	int m_col;												// columns of CMatrix
	MATRIXELEMENTTYPE   m_vec;        						// data of matrix

public:
	// Constructor, zero-length  vector
	CMatrix();				

	// Constructor, vector of length MaxSize
	CMatrix(int row, int col);		

	// Copy constructor	
	CMatrix(const CMatrix& mat);                   
	~CMatrix();												// Destructor        
	CMatrix& operator =(const CMatrix& mat);				// Assignment operator  
	CMatrix& operator +=(const CMatrix& mat);

	// access to the elements of the CMatrix                     
	std::vector<double>& operator [](int row);
	const std::vector<double>& operator [](int row) const;

	// validate whether the size of vector is legal
	int CheckRowCol(int row, int col);
	
	// get Rows and Columns of CMatrix
	int GetRows() const;           
	int GetCols() const;           
};

#endif //__MATRIX_H__