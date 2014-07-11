#include "Loka3math.h"
#include "Vector3.h"

using namespace Vector3;

namespace Loka3
{
	static const Vect3d ortA( 1, 0, 0 );
	static const Vect3d ortB( cos(2*PI/3), sin(2*PI/3), 0 );
	static const Vect3d ortE( cos(4*PI/3), sin(4*PI/3), 0 );
	
	const L3 normalize( const L3 & n )
	{
		return n / norm(n);
	}

	Vect3d L3_to_vect( const L3 & n )
	{
		return ortA * n.a() + ortB * n.b() + ortE * n.e();
	}

	L3 vect_to_L3( const Vect3d & v )
	{
		const double e = 1e-9;
		const double x = v.x;
		const double y = v.y;

		// Point is in origin
		if( fabs(x) < e && fabs(y) < e )
			return L3(0,0,0);
		
		// Point is on ortA
		if( x >= 0 && fabs(y) < e )
			return A(x);

		// Point is on ortB
		const double ortB_deviation = ortB.y * x - ortB.x * y;
		if( fabs(ortB_deviation) < e )
			return B( sqrt(x*x+y*y) );

		// Point is on ortE
		const double ortE_deviation = - ortE.y * x + ortE.x * y;
		if( fabs(ortE_deviation) < e )
			return E( sqrt(x*x+y*y) );

		// Point between A and B orts
		if( y >= 0 && ortB_deviation >= 0 )
		{
			return L3( x - y * ortB.x / ortB.y,  y / ortB.y, 0 );
		}

		// Point between E and A orts
		if( y <= 0 && ortE_deviation >= 0 )
		{
			return L3( x - y * ortE.x / ortE.y, 0, y / ortE.y );
		}

		// Point between B and E orts
		const double k = ortB_deviation / ( ortE.x * ortB.y - ortE.y * ortB.x );
		return L3( 0, (x - k * ortE.x) / ortB.x, k );
	}

	double exponent_period()
	{
		static const double p = 4 * PI / sqrt( 3 );
		return p;
	}

	bool near_equal( const L3 & n1, const L3 & n2, const double tol /*= 1e-12*/ )
	{
		return (fabs(n1.a()-n2.a()) < tol) && (fabs(n1.b()-n2.b()) < tol) && (fabs(n1.e()-n2.e()) < tol);
	}

	bool near_zero( const L3 & n1, const double tol /* = 1e-12 */ )
	{
		return (fabs(n1.a()) < tol) && (fabs(n1.b()) < tol) && (fabs(n1.e()) < tol);
	}

	L3 lin_comb( const L3 & n1, const L3 & n2, const L3 & k1, const L3 & k2 )
	{
		return n1 * k1 + n2 * k2;
	}

	L3 rotate( const L3 & pt, const L3 & center, const L3 & r )
	{
		return lin_comb( pt, center, r, E(1) + opp(r) );
	}

	L3 exponent( const L3 & n, const int maxiters /*= 500*/, const double tol /*= 1e-12*/ )
	{
		L3 summand(0,0,1);
		L3 sum(0,0,1);
		int iter = 1;
		while( (iter < maxiters) && !near_zero(summand,tol) )
		{
			summand *= n / double(iter);
			sum += summand;
			iter++;
		}
		return sum;
	}

	// calculate a table of values for func Exp[omega*t+fi]
	std::vector<L3> exponent_table( const L3 & omega, const L3 & fi, const double t0, const double t1, const int n )
	{
		using namespace std;
		vector<L3> tbl;
		tbl.reserve( n );
		const double dt = ( t1 - t0 ) / ( n - 1 );
		double t = t0;
		tbl.push_back( exponent(omega*t0+fi) );
		for( int i = 1; i < n-1; ++i )
		{
			t += dt;
			tbl.push_back( exponent(omega*t+fi) );
		}
		tbl.push_back( exponent(omega*t1+fi) );
		return tbl;
	}
}