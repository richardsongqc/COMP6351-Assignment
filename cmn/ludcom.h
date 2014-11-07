#pragma once

#ifndef __LUDCOMP_H__
#define __LUDCOMP_H__

#include "..\\cmn\\Vector.h"
#include "..\\cmn\\Matrix.h"


extern void LUDecompTridiagnoal(
	CMatrix a, 
	CVector b, 
	CVector& x, CVector& g);

	
// Tridiagnonal Systems (Vector Form)
extern void LUDecompTridiagnoal(
	CVector l1, CVector l2, CVector l3,							// vector a, b, c in Tridiagnonal matrix
	CVector b,													
	CVector& x, CVector& g)	;
	
// Tridiagnonal Systems (Vector Form)
extern void LUDecompTridiagnoal(
	CVector l1, CVector l2, CVector l3,							// vector a, b, c in Tridiagnonal matrix
	CVector b,
	CVector& x );

#endif //__LUDCOMP_H__