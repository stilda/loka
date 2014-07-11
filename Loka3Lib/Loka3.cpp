#pragma once

#include "Loka3.h"
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>

namespace Loka3
{
	using namespace std;

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

	void L3::compensate()
	{
		const double mx = min( min( _a, _b ), _e );
		_a -= mx;
		_b -= mx;
		_e -= mx;
	}

	void L3::check_positive()
	{
		const double eps = 1e-12;
		if( _a < -eps || _b < -eps || _e < -eps )
			throw "L3 Number Has To Have All Coeficients Positive";
	}

	double L3::operator[] (int i) const
	{
		assert( i >= 0 && i < 3 );
		switch( i )
		{
		case 0: return _a;
		case 1: return _b;
		case 2: return _e;
		}
		return _e; // never get here
	}


	L3 & L3::operator += ( const L3 & n )
	{
		_a += n.a();
		_b += n.b();
		_e += n.e();
		compensate();
		return *this;
	}

	L3 & L3::operator *= ( const L3 & n )
	{
		*this = *this * n;
		return *this;
	}

	L3 & L3::operator /= ( const L3 & n )
	{
		*this = *this / n;
		return *this;
	}

	L3 & L3::operator *= ( const double d )
	{
		assert( d >= -1e-12 );
		_a *= d;
		_b *= d;
		_e *= d;
		return *this;
	}

	L3 & L3::operator /= ( const double d )
	{
		assert( d >= 1e-12 );
		_a /= d;
		_b /= d;
		_e /= d;
		return *this;
	}

	string L3::show() const
	{
		auto NearZero = []( const double d ) { return (d > 0) ? (d < 1e-9) : (-d < 1e-9); };

		stringstream o;
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

	string L3::showb() const
	{
		return "(" + show() + ')';
	}

	L3 operator + ( const L3 & n, const L3 & m )
	{
		return L3( n.a() + m.a(), n.b() + m.b(), n.e() + m.e() );
	}

	L3 operator * ( const L3 & n, const L3 & m )
	{
		return L3( 
			n.a()*m.e() + m.a()*n.e() + n.b()*m.b(), 
			n.b()*m.e() + m.b()*n.e() + n.a()*m.a(), 
			n.e()*m.e() + n.a()*m.b() + n.b()*m.a() );
	}

	L3 operator * ( const L3 & n, const double d )
	{
		return L3( n.a()*d, n.b()*d, n.e()*d );
	}

	L3 operator * ( const double d,const L3 & n )
	{
		return L3( n.a()*d, n.b()*d, n.e()*d );
	}

	L3 operator / ( const L3 & n, const double d )
	{
		return n * (1/d);
	}

	L3 operator / ( const double d, const L3 & n )
	{
		return d * inv(n);
	}

	L3 operator / ( const L3 & n, const L3 & m )
	{
		return n * inv(m);
	}

	L3 conjugate1( const L3 & n )
	{
		return L3( n.b(), n.a(), n.e() );
	}

	L3 conjugate2( const L3 & n )
	{
		return L3( 0, 0, n.a() + n.b() + n.e() );
	}

	double norm3( const L3 & n )
	{
		return n.a()*n.a()*n.a() + n.b()*n.b()*n.b() + n.e()*n.e()*n.e();
	}

	double norm3c( const L3 & n )
	{
		return (n*conjugate1(n)*conjugate2(n)).e();
	}

	double norm( const L3 & n )
	{
		return pow( norm3(n), 1.0/3.0 );
	}

	L3 inv( const L3 & n )
	{
		throw "inv(L3) not implemented";
	}

	L3 opp( const L3 & n )
	{
		const double mx = max( max( n.a(), n.b() ), n.e() );
		return L3( mx - n.a(), mx - n.b(), mx - n.e() );
	}

}