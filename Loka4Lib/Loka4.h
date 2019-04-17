#include <string>

namespace Loka4
{

	class L4
	{
		double _a, _b, _c, _e;

		void compensate();
		void check_positive();
	public:

		L4() : _a(0), _b(0), _c(0), _e(0) 
		{
		}

		L4( const double a, const double b, const double c, const double e ) : _a(a), _b(b), _c(c), _e(e) 
		{
			compensate();
			check_positive();
		}

		L4( const L4 & n ) : _a(n.a()), _b(n.b()), _c(n.c()), _e(n.e())
		{
			check_positive();
		}

		~L4() {}

		double a() const { return _a; }
		double b() const { return _b; }
		double c() const { return _c; }
		double e() const { return _e; }

		double operator[] ( int i ) const;
		L4 & operator += ( const L4 & n );
		L4 & operator *= ( const L4 & n );
		L4 & operator /= ( const L4 & n );
		L4 & operator *= ( const double d );
		L4 & operator /= ( const double d );

		std::string show() const;
		std::string showb() const;
	};

	L4 A( const double a );
	L4 B( const double b );
	L4 C( const double c );
	L4 E( const double e );

	L4 operator + ( const L4 & n, const L4 & m );
	L4 operator * ( const L4 & n, const L4 & m );
	L4 operator * ( const L4 & n, const double d );
	L4 operator * ( const double d,const L4 & n );
	L4 operator / ( const L4 & n, const double d );
	L4 operator / ( const double d, const L4 & n );
	L4 operator / ( const L4 & n, const L4 & m );

	L4 conjugate1( const L4 & n );
	L4 conjugate2( const L4 & n );
	L4 conjugate3( const L4 & n );
	double norm4( const L4 & n );
	double norm4c( const L4 & n );
	double norm( const L4 & n );
	L4 inv( const L4 & n ); // inverse: n*inv(n)=1
	L4 opp( const L4 & n ); // opposite: n+opp(n)=0

}
