#pragma once

#include <Windows.h>
// math
#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include "Vector3.h"

namespace Loka3
{
	using namespace Vector3;

	Vect3 ortX(1,0,0), ortY(cos(2*M_PI/3),sin(2*M_PI/3),0), ortZ(cos(4*M_PI/3),sin(4*M_PI/3),0);

	bool NearZero( const double d )
	{
		return (d > 0) ? (d < 1e-9) : (-d < 1e-9);
	}

	class L3
	{
		double _a, _b, _e;

		void compensate()
		{
			const double mx = min( min( _a, _b ), _e );
			_a -= mx;
			_b -= mx;
			_e -= mx;
		}

	public:
		L3() : _a(0), _b(0), _e(0) {}
		L3( const double a, const double b, const double e ) : _a(a), _b(b), _e(e) 
		{
			compensate();
		}
		L3( const L3 & n ) : _a(n.a()), _b(n.b()), _e(n.e())
		{
			compensate();
		}

		~L3() {}

		double a() const { return _a; }
		double b() const { return _b; }
		double e() const { return _e; }

		std::string show() const
		{
			std::ostringstream o;
			if( !NearZero(_e) )
			{
				o << _e;
				if( !NearZero(_a) || !NearZero(_b) )
					o << " ";
			}
			if( !NearZero(_a) )
			{
				o << _a << "A";
				if( !NearZero(_b) )
					o << " ";
			}
			if( !NearZero(_b) )
				o << _b << "B";
			if( NearZero(_a) && NearZero(_b) && NearZero(_e) )
				o << "0";
			return o.str();
		}

		std::string showb() const
		{
			return "(" + show() + ')';
		}
	};

	L3 A( const double a ) 
	{ 
		return L3(a,0,0); 
	}
	L3 B( const double b ) 
	{ 
		return L3(0,b,0); 
	}
	L3 E( const double e ) 
	{ 
		return L3(0,0,e); 
	}

	const L3 operator + (const L3 & n, const L3 & m )
	{
		return L3( n.a() + m.a(), n.b() + m.b(), n.e() + m.e() );
	}

	const L3 operator * (const L3 & n, const L3 & m )
	{
		return L3( 
			n.a()*m.e() + m.a()*n.e() + n.b()*m.b(), 
			n.b()*m.e() + m.b()*n.e() + n.a()*m.a(), 
			n.e()*m.e() + n.a()*m.b() + n.b()*m.a() );
	}

	const L3 operator * (const L3 & n, const double d )
	{
		return L3( n.a()*d, n.b()*d, n.e()*d );
	}

	const L3 operator * (const double d,const L3 & n )
	{
		return L3( n.a()*d, n.b()*d, n.e()*d );
	}

	const L3 conjugate1( const L3 & n )
	{
		return L3( n.b(), n.a(), n.e() );
	}

	const L3 conjugate2( const L3 & n )
	{
		return L3( 0, 0, n.a() + n.b() + n.e() );
	}

	const double abs3( const L3 & n )
	{
		return n.a()*n.a()*n.a() + n.b()*n.b()*n.b() + n.e()*n.e()*n.e();
	}

	const double abs3c( const L3 & n )
	{
		return (n*conjugate1(n)*conjugate2(n)).e();
	}

	const double abs( const L3 & n )
	{
		return pow( abs3(n), 1.0/3.0 );
	}

	const L3 normalize( const L3 & n )
	{
		return n * (1.0/abs(n));
	}

	Vect3 tovect( const L3 & n )
	{
		return ortX*n.a() + ortY*n.b() + ortZ*n.e();
	}

	bool NearEqual( const L3 & n1, const L3 & n2, const double tol = 1e-6 )
	{
		return (fabs(n1.a()-n2.a()) < tol) && (fabs(n1.b()-n2.b()) < tol) && (fabs(n1.e()-n2.e()) < tol);
	}

	bool NearZero( const L3 & n1, const double tol = 1e-6 )
	{
		return (fabs(n1.a()) < tol) && (fabs(n1.b()) < tol) && (fabs(n1.e()) < tol);
	}

	L3 exponent( const L3 & n, const int maxiters = 100, const double tol = 1e-9 )
	{
		L3 summand(0,0,1);
		L3 sum(0,0,1);
		int iter = 1;
		while( (iter < maxiters) && (!NearZero(summand,tol)) )
		{
			summand = summand * n * (1.0/iter);
			sum = sum + summand;
			iter++;
		}
		return sum;
	}
}