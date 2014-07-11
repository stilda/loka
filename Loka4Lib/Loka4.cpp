#pragma once

#include "Loka4.h"
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>

namespace Loka4
{
	using namespace std;

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

	void L4::compensate()
	{
		const double mx = min( min( _a, _b ), min( _c, _e ) );
		_a -= mx;
		_b -= mx;
		_c -= mx;
		_e -= mx;
	}

	void L4::check_positive()
	{
		const double eps = 1e-12;
		if( _a < -eps || _b < -eps || _c < -eps || _e < -eps )
			throw "L4 Number Has To Have All Coeficients Positive";
	}

	double L4::operator[] ( int i ) const
	{
		assert( i >= 0 && i < 4 );
		switch( i )
		{
		case 0: return _a;
		case 1: return _b;
		case 2: return _c;
		case 3: return _e;
		}
		return _e; // never get here
	}

	L4 & L4::operator += ( const L4 & n )
	{
		_a += n.a();
		_b += n.b();
		_c += n.c();
		_e += n.e();
		compensate();
		return *this;
	}

	L4 & L4::operator *= ( const L4 & n )
	{
		*this = *this * n;
		return *this;
	}

	L4 & L4::operator /= ( const L4 & n )
	{
		*this = *this / n;
		return *this;
	}

	L4 & L4::operator *= ( const double d )
	{
		assert( d >= -1e-12 );
		_a *= d;
		_b *= d;
		_c *= d;
		_e *= d;
		return *this;
	}

	L4 & L4::operator /= ( const double d )
	{
		assert( d >= 1e-12 );
		*this *= 1 / d;
		return *this;
	}

	string L4::show() const
	{
		auto NearZero = []( const double d ) { return (d > 0) ? (d < 1e-9) : (-d < 1e-9); };

		ostringstream o;
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

	string L4::showb() const
	{
		return "(" + show() + ')';
	}

	L4 operator + ( const L4 & n, const L4 & m )
	{
		return L4( n.a() + m.a(), n.b() + m.b(), n.c() + m.c(), n.e() + m.e() );
	}

	L4 operator * ( const L4 & n, const L4 & m )
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

	L4 operator * ( const L4 & n, const double d )
	{
		return L4( n.a()*d, n.b()*d, n.c()*d, n.e()*d );
	}

	L4 operator * ( const double d,const L4 & n )
	{
		return n * d;
	}

	L4 operator / ( const L4 & n, const double d )
	{
		return n * (1/d);
	}

	L4 operator / ( const double d, const L4 & n )
	{
		return d * inv(n);
	}

	L4 operator / ( const L4 & n, const L4 & m )
	{
		return n * inv(m);
	}

	L4 conjugate1( const L4 & n )
	{
		return L4( n.b(), n.c(), n.a(), n.e() );
	}

	L4 conjugate2( const L4 & n )
	{
		return L4( n.c(), n.a(), n.b(), n.e() );
	}

	L4 conjugate3( const L4 & n )
	{
		return L4( 0, 0, 0, n.e() + n.a() + n.b() + n.c() );
	}

	double norm4( const L4 & n )
	{
		return n.a()*n.a()*n.a()*n.a() + n.b()*n.b()*n.b()*n.b() + n.c()*n.c()*n.c()*n.c() + n.e()*n.e()*n.e()*n.e();
	}

	double norm4c( const L4 & n )
	{
		return (n*conjugate1(n)*conjugate2(n)*conjugate3(n)).e();
	}

	double norm( const L4 & n )
	{
		return pow( norm4(n), 0.25 );
	}

	L4 inv( const L4 & n )
	{
		throw "inv(L4) not implemented";
	}

	L4 opp( const L4 & n )
	{
		const double mx = max( max( max( n.a(), n.b() ), n.c() ), n.e() );
		return L4( mx - n.a(), mx - n.b(), mx - n.c(), mx - n.e() );
	}

}