#include "Vector.h"


CVector::CVector() : 
	m_maxSize(1)
{
	m_vec.resize(m_maxSize);
}

CVector::CVector(int MaxSize) :
	m_maxSize(MaxSize)
{
	m_vec.resize(m_maxSize);
}

CVector::CVector(const CVector& vec)
{
	m_maxSize = vec.m_maxSize;
	m_vec     = vec.m_vec;
}

CVector::~CVector()
{
	m_vec.clear();
}

CVector& CVector::operator =(const CVector& vec) 				// Assignment operator 
{
	m_maxSize = vec.m_maxSize;
	m_vec = vec.m_vec;

	return *this;
}

bool CVector::operator==(double dbl) const
{
	for (std::vector<double>::const_iterator iter = m_vec.begin();
		iter != m_vec.end();
		iter++)
	{
		if (*iter != dbl)
		{
			return false;
		}
	}

	return true;
}

bool CVector::operator<(double dbl) const
{
	for (std::vector<double>::const_iterator iter = m_vec.begin();
		iter != m_vec.end();
		iter++)
	{
		if ( abs(*iter) >= dbl)
		{
			return false;
		}
	}

	return true;
}

// access to the elements of the CMatrix                     
double& CVector::operator [](int nPos)					
{
	return m_vec[nPos];
}

const double& CVector::operator [](int nPos) const
{
	return m_vec[nPos];
}

int CVector::CheckSize(int nSize)
{
	return ( (nSize >= 0) && (nSize < m_maxSize) ) ? 1 : 0; 
}

// get Rows and Columns of CMatrix
int CVector::GetSize() const
{
	return m_maxSize;
}