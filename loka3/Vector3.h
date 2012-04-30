#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

namespace Vector3
{
	class Vect3
	{
	public:
		double x, y, z;
		Vect3() : x(0), y(0), z(0) {}
		Vect3( const double a, const double b, const double c ) : x(a), y(b), z(c) {}
		~Vect3() {}
	};

	Vect3 operator + ( const Vect3 & v, const Vect3 & u )
	{
		return Vect3( v.x+u.x, v.y+u.y, v.z+u.z );
	}

	Vect3 operator - ( const Vect3 & v, const Vect3 & u )
	{
		return Vect3( v.x-u.x, v.y-u.y, v.z-u.z );
	}

	Vect3 operator * ( const Vect3 & v, const double d )
	{
		return Vect3( v.x*d, v.y*d, v.z*d );
	}

	double abs( const Vect3 & v )
	{
		return sqrt( v.x*v.x + v.y*v.y + v.z*v.z );
	}

	Vect3 normalize( const Vect3 & v )
	{
		return v * (1.0/abs(v));
	}
}

