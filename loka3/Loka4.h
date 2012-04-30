#pragma once

#include <Windows.h>
// math
#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include "Vector3.h"


namespace Loka4
{
	using namespace Vector3;

	bool NearZero( const double d )
	{
		return (d > 0) ? (d < 1e-9) : (-d < 1e-9);
	}

	class L4
	{
		double _a, _b, _c, _e;

		void compensate()
		{
			const double mx = min( min( _a, _b ), min( _c, _e ) );
			_a -= mx;
			_b -= mx;
			_c -= mx;
			_e -= mx;
		}

	public:
		L4() : _a(0), _b(0), _c(0), _e(0) {}
		L4( const double a, const double b, const double c, const double e ) : _a(a), _b(b), _c(c), _e(e) 
		{
			compensate();
		}
		L4( const L4 & n ) : _a(n.a()), _b(n.b()), _c(n.c()), _e(n.e())
		{
			compensate();
		}

		~L4() {}

		double a() const { return _a; }
		double b() const { return _b; }
		double c() const { return _c; }
		double e() const { return _e; }

		std::string show() const
		{
			std::ostringstream o;
			if( !NearZero(_e) )
			{
				o << _e;
				if( !NearZero(_a) || !NearZero(_b) || !NearZero(_c) )
					o << " ";
			}
			if( !NearZero(_a) )
			{
				o << _a << "A";
				if( !NearZero(_b) || !NearZero(_c) )
					o << " ";
			}
			if( !NearZero(_b) )
			{
				o << _b << "B";
				if( !NearZero(_c) )
					o << " ";
			}
			if( !NearZero(_c) )
			{
				o << _c << "C";
			}

			if( NearZero(_a) && NearZero(_b) && NearZero(_c) && NearZero(_e) )
				o << "0";
			return o.str();
		}

		std::string showb() const
		{
			return "(" + show() + ')';
		}
	};

	L4 A( const double a ) 
	{ 
		return L4(a,0,0,0); 
	}
	L4 B( const double b ) 
	{ 
		return L4(0,b,0,0); 
	}
	L4 C( const double c ) 
	{ 
		return L4(0,0,c,0); 
	}
	L4 E( const double e ) 
	{ 
		return L4(0,0,0,e); 
	}

	const L4 operator + (const L4 & n, const L4 & m )
	{
		return L4( n.a() + m.a(), n.b() + m.b(), n.c() + m.c(), n.e() + m.e() );
	}

	const L4 operator * (const L4 & n, const L4 & m )
	{
		/*	(x,y,z,p)*(u,v,w,k) = 
			x*k + y*w + z*v + p*u, 
			x*u + y*k + z*w + p*v, 
			x*v + y*u + z*k + p*w,
			x*w + y*v + z*u + p*k */
		return L4( 
			n.a()*m.e() + n.b()*m.c() + n.c()*m.b() + n.e()*m.a(), 
			n.a()*m.a() + n.b()*m.e() + n.c()*m.c() + n.e()*m.b(), 
			n.a()*m.b() + n.b()*m.a() + n.c()*m.e() + n.e()*m.c(),
			n.a()*m.c() + n.b()*m.b() + n.c()*m.a() + n.e()*m.e() );
	}

	const L4 operator * (const L4 & n, const double d )
	{
		return L4( n.a()*d, n.b()*d, n.c()*d, n.e()*d );
	}

	const L4 operator * (const double d,const L4 & n )
	{
		return n*d;
	}

	const L4 conjugate1( const L4 & n )
	{
		return L4( n.b(), n.c(), n.e(), n.a() );
	}

	const L4 conjugate2( const L4 & n )
	{
		return L4( n.c(), n.e(), n.a(), n.b() );
	}

	const L4 conjugate3( const L4 & n )
	{
		return L4( n.e(), n.a(), n.b(), n.c() );
	}

	const double abs4( const L4 & n )
	{
		return n.a()*n.a()*n.a()*n.a() + n.b()*n.b()*n.b()*n.b() + n.c()*n.c()*n.c()*n.c() + n.e()*n.e()*n.e()*n.e();
	}

	const double abs4c( const L4 & n )
	{
		return (n*conjugate1(n)*conjugate2(n)*conjugate3(n)).e();
	}

	const double abs( const L4 & n )
	{
		return pow( abs4(n), 0.25 );
	}

	const L4 normalize( const L4 & n )
	{
		return n * (1.0/abs(n));
	}

	Vect3 tovect( const L4 & n )
	{
		const static Vect3 ortA = normalize(Vect3(1,-1,-1));
		const static Vect3 ortB = normalize(Vect3(-1,1,-1));
		const static Vect3 ortC = normalize(Vect3(-1,-1,1));
		const static Vect3 ortE = normalize(Vect3(1,1,1));

		return ortA*n.a() + ortB*n.b() + ortC*n.c() + ortE*n.e();
	}

	bool NearEqual( const L4 & n1, const L4 & n2, const double tol = 1e-6 )
	{
		return (fabs(n1.a()-n2.a()) < tol) && (fabs(n1.b()-n2.b()) < tol) && (fabs(n1.c()-n2.c()) < tol) && (fabs(n1.e()-n2.e()) < tol);
	}

	bool NearZero( const L4 & n1, const double tol = 1e-6 )
	{
		return (fabs(n1.a()) < tol) && (fabs(n1.b()) < tol) && (fabs(n1.c()) < tol) && (fabs(n1.e()) < tol);
	}

	L4 exponent( const L4 & n, const int maxiters = 100, const double tol = 1e-9 )
	{
		L4 summand(0,0,0,1);
		L4 sum(0,0,0,1);
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