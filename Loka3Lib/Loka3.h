#ifndef __LOKA3_HEADER__
#define __LOKA3_HEADER__

#include <string>

namespace Loka3
{
	// There is no negative values here. There is no subtraction and no unary mines.
	// All scalar values (doubles) have to be positive. Class forces this.
	class L3
	{
		double _a, _b, _e;

		void compensate();
		void check_positive();
	public:
		L3() : _a(0), _b(0), _e(0) 
		{
		}

		L3( const double a, const double b, const double e ) : _a(a), _b(b), _e(e) 
		{
			check_positive();
			compensate();
		}

		L3( const L3 & n ) : _a(n.a()), _b(n.b()), _e(n.e())
		{
			check_positive();
		}

		~L3() {}

		double a() const { return _a; }
		double b() const { return _b; }
		double e() const { return _e; }

		double operator[] ( int i ) const;
		L3 & operator += ( const L3 & n );
		L3 & operator *= ( const L3 & n );
		L3 & operator /= ( const L3 & n );
		L3 & operator *= ( const double d );
		L3 & operator /= ( const double d );

		std::string show() const;
		std::string showb() const;
	};

	L3 A( const double a );
	L3 B( const double b );
	L3 E( const double e );

	L3 operator + ( const L3 & n, const L3 & m );
	L3 operator * ( const L3 & n, const L3 & m );
	L3 operator * ( const L3 & n, const double d );
	L3 operator * ( const double d,const L3 & n );
	L3 operator / ( const L3 & n, const double d );
	L3 operator / ( const double d, const L3 & n );
	L3 operator / ( const L3 & n, const L3 & m );

	L3 conjugate1( const L3 & n );
	L3 conjugate2( const L3 & n );
	double norm3( const L3 & n );
	double norm3c( const L3 & n );
	double norm( const L3 & n );
	L3 inv( const L3 & n ); // inverse: n*inv(n)=1
	L3 opp( const L3 & n ); // opposite: n+opp(n)=0
}

#endif
