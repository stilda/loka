#ifndef __LOKA3_MATH_HEADER__
#define __LOKA3_MATH_HEADER__

#include "Loka3.h"
//#include "Vector3.h"
#include <vector>

namespace Loka3
{
	const double PI = 3.1415926535897932384626433832795;

	const L3 normalize( const L3 & n );
//	Vector3::Vect3d L3_to_vect( const L3 & n );
	// Convert dekart orthogonal coordinates of a point to L3 representation.
//	L3 vect_to_L3( const Vector3::Vect3d & v );
	double exponent_period();
	bool near_equal( const L3 & n1, const L3 & n2, const double tol = 1e-12 );
	bool near_zero( const L3 & n1, const double tol = 1e-12 );
	L3 exponent( const L3 & n, const int maxiters = 500, const double tol = 1e-12 );
    L3 logarithm( const L3 & n, const int maxiters = 500, const double tol = 1e-12 );
	L3 lin_comb( const L3 & n1, const L3 & n2, const L3 & k1, const L3 & k2 );
	L3 rotate( const L3 & pt, const L3 & center, const L3 & r );
	// calculate a table of values for func Exp[omega*t+fi]
	std::vector<L3> exponent_table( const L3 & omega, const L3 & fi, const double t0, const double t1, const int n );
    L3 pow_e(double n); // n^E
    L3 pow_a(double n); // n^A
    L3 pow_b(double n); // n^B
    L3 ratio(double a, double b, double e); // e * a^A * b^B
}

#endif
