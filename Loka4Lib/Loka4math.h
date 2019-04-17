#include "Loka4.h"
//#include "Vector3.h"
#include <vector>

namespace Loka4
{	
	const double PI = 3.1415926535897932384626433832795;

	const L4 normalize( const L4 & n );
//	Vector3::Vect3d L4_to_vect( const L4 & n );
	// Convert dekart orthogonal coordinates of a point to L4 representation.
//	L4 vect_to_L4( const Vector3::Vect3d & v );
	double exponent_period();
	bool near_equal( const L4 & n1, const L4 & n2, const double tol = 1e-12 );
	bool near_zero( const L4 & n1, const double tol = 1e-12 );
	L4 exponent( const L4 & n, const int maxiters = 500, const double tol = 1e-12 );
	//L4 lin_comb( const L4 & n1, const L4 & n2, const L4 & k1, const L4 & k2 );
	//L4 rotate( const L4 & pt, const L4 & center, const L4 & r );
	// calculate a table of values for func Exp[omega*t+fi]
	std::vector<L4> exponent_table( const L4 & omega, const L4 & fi, const double t0, const double t1, const int n );
}
