#include "Loka4math.h"
#include "Vector3.h"

using namespace Vector3;

namespace Loka4
{
	const static Vect3d ortA = normalize( Vect3d(1,-1,-1) );
	const static Vect3d ortB = normalize( Vect3d(-1,1,-1) );
	const static Vect3d ortC = normalize( Vect3d(-1,-1,1) );
	const static Vect3d ortE = normalize( Vect3d(1,1,1) );

	const L4 normalize( const L4 & n )
	{
		return n / norm(n);
	}

	Vect3d L4_to_vect( const L4 & n )
	{
		return ortA*n.a() + ortB*n.b() + ortC*n.c() + ortE*n.e();
	}

	L4 vect_to_L4( const Vect3d & v )
	{
		throw "vect_to_L4 not implemented";
		return L4();
	}

	double exponent_period()
	{
		throw "Loka4::exponent_period() not implemented";
		return 1;
	}

	bool near_equal( const L4 & n1, const L4 & n2, const double tol /*= 1e-12*/ )
	{
		return (fabs(n1.a()-n2.a()) < tol) && (fabs(n1.b()-n2.b()) < tol) && (fabs(n1.c()-n2.c()) < tol) && (fabs(n1.e()-n2.e()) < tol);
	}

	bool near_zero( const L4 & n1, const double tol /*= 1e-12*/ )
	{
		return (fabs(n1.a()) < tol) && (fabs(n1.b()) < tol) && (fabs(n1.c()) < tol) && (fabs(n1.e()) < tol);
	}

	L4 exponent( const L4 & n, const int maxiters /*= 500*/, const double tol /*= 1e-12*/ )
	{
		L4 summand(0,0,0,1);
		L4 sum(0,0,0,1);
		int iter = 1;
		while( (iter < maxiters) && !near_zero(summand,tol) )
		{
			summand *= n / double(iter);
			sum += summand;
			++iter;
		}
		return sum;
	}

	// calculate a table of values for func Exp[omega*t+fi]
	std::vector<L4> exponent_table( const L4 & omega, const L4 & fi, const double t0, const double t1, const int n )
	{
		using namespace std;
		vector<L4> tbl;
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