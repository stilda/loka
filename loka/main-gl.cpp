// loka3.cpp : Defines the entry point for the console application.
//

#include <Windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "glaux.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "../Loka3Lib/Loka3.h"
#include "../Loka3Lib/Loka3math.h"
#include "../Loka4Lib/Loka4.h"
#include "../Loka4Lib/Loka4math.h"
#include "Vector3.h"
#include "LinearAlgebra.h"

std::vector<Loka3::L3> Add( const std::vector<Loka3::L3> & v1, const std::vector<Loka3::L3> & v2 )
{
	std::vector<Loka3::L3> r(v1);
	for( size_t i = 0; i < v2.size(); ++i )
		r[i] += v2[i];
	return r;
}

void DrawCircle2D( const double x, const double y, const double r )
{
	glPushMatrix();
	glTranslated( x, y, 0 );
	glBegin(GL_POLYGON);
	for( int i = 0; i < 100; i++ ) 
	{ 
		double angle = i * 2.0 * M_PI * 0.01; 
		glVertex3d( cos(angle) * r, sin(angle) * r, 0 ); 
	} 
	glEnd();
	glPopMatrix();
}

void DrawCurve( const std::vector<Loka3::L3> & curve, const bool points )
{
	using namespace Vector3;

	// Curve
	glLineWidth( 2 );
	if( !points )
		glBegin( GL_LINE_STRIP );
	else
		glBegin( GL_POINTS );
	std::for_each( begin(curve), end(curve), [](const Loka3::L3 & l3){ 
		const Vect3d v = Loka3::L3_to_vect(l3);
		glVertex3d( v.x, v.y, v.z ); 
	} );
	glEnd();

	// Curve start
	glLineWidth(1);
	glColor3d(1,1,0);
	const Vect3d v0 = Loka3::L3_to_vect(curve[0]);
	DrawCircle2D( v0.x, v0.y, 0.01 );

	// Curve direction
	Vect3d dir = Loka3::L3_to_vect(curve[1]) - v0;
	dir = dir * (0.2 / abs(dir));
	const Vect3d ptEnd = v0 + dir;
	glBegin( GL_LINES );
	glVertex3d( v0.x, v0.y, v0.z );
	glVertex3d( ptEnd.x, ptEnd.y, ptEnd.z );
	glEnd();
}

void DrawCurve( const std::vector<Loka4::L4> & curve, const bool points )
{
	using namespace Vector3;

	// Curve
	glLineWidth( 2 );
	if( !points )
		glBegin( GL_LINE_STRIP );
	else
		glBegin( GL_POINTS );
	std::for_each( begin(curve), end(curve), [](const Loka4::L4 & l4){ 
		const Vect3d v = Loka4::L4_to_vect(l4);
		glVertex3d( v.x, v.y, v.z ); 
	} );
	glEnd();

	// Curve start
	glLineWidth(1);
	glColor3d(1,1,0);
	const Vect3d v0 = Loka4::L4_to_vect(curve[0]);
	//DrawCircle2D( v0.x, v0.y, v0.z );

	// Curve direction
	Vect3d dir = Loka4::L4_to_vect(curve[1]) - v0;
	dir = dir * (0.2 / abs(dir));
	const Vect3d ptEnd = v0 + dir;
	glBegin( GL_LINES );
	glVertex3d( v0.x, v0.y, v0.z );
	glVertex3d( ptEnd.x, ptEnd.y, ptEnd.z );
	glEnd();
}

void DrawDekartProjectionCurve( const std::vector<Loka3::L3> & curve )
{
	using namespace Vector3;

	// Curve
	glLineWidth( 2 );
	glBegin( GL_LINE_STRIP );
	//glBegin( GL_POINTS );
	double t = 0;
	const double dt = 1.0 / curve.size();
	std::for_each( begin(curve), end(curve), [&t,dt](const Loka3::L3 & l3){ 
		const Vect3d v = Loka3::L3_to_vect(l3);
		t += dt;
		glVertex3d( t, v.x, 0 ); 
	} );
	glEnd();
}

void DrawPolarityProjectionCurve( const std::vector<Loka3::L3> & curve, const int coordinate_index )
{
	using namespace Vector3;

	// Curve
	glLineWidth( 2 );
	glBegin( GL_LINE_STRIP );
	//glBegin( GL_POINTS );
	double t = 0;
	const double dt = 1.0 / curve.size();
	std::for_each( begin(curve), end(curve), [&t,dt,coordinate_index](const Loka3::L3 & l3){ 
		t += dt;
		glVertex3d( t, l3[coordinate_index] / 3, 0 ); 
	} );
	glEnd();
}

void DrawCurve( const std::vector<Vector3::Vect3d> & curve_xyz )
{
	using namespace Vector3;

	// Curve
	glLineWidth( 2 );
	glBegin( GL_LINE_STRIP );
	//glBegin( GL_POINTS );
	std::for_each( begin(curve_xyz), end(curve_xyz), 
		[](const Vect3d & v){ glVertex3d( v.x, v.y, v.z ); } );
	glEnd();

	// Curve start
	glLineWidth(1);
	glColor3d(1,1,0);
	DrawCircle2D( curve_xyz[0].x, curve_xyz[0].y, 0.01 );

	// Curve direction
	Vect3d dir = curve_xyz[1] - curve_xyz[0];
	dir = dir * (0.2 / abs(dir));
	Vect3d ptEnd = curve_xyz[0] + dir;
	glBegin( GL_LINES );
	glVertex3d( curve_xyz[0].x, curve_xyz[0].y, curve_xyz[0].z );
	glVertex3d( ptEnd.x, ptEnd.y, ptEnd.z );
	glEnd();
}

void DrawRotationTrace( const Loka3::L3 & n )
{
	using namespace Loka3;
	using namespace Vector3;

	// Generate curve
	std::vector<L3> curve;
	L3 acc(0,0,1);
	for( int t = 0; t < 100; ++t )
	{
		curve.push_back( acc );		
		acc = acc * n;
		//std::cout << norm3(acc) << std::endl;
	}

	// Convert it to vectors
	std::vector<Vect3d> curve_xyz;
	for( auto it = curve.begin(); it != curve.end(); ++it )
	{
		curve_xyz.push_back( L3_to_vect(*it) );
	}

	DrawCurve( curve_xyz );
}

void CALLBACK resize(int width,int height)
{
	glViewport(0,0,width,height);
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glOrtho(-1,2, -1,1, 2,12);
	gluLookAt( 0,0,5, 0,0,0, 0,1,0 );
	glMatrixMode( GL_MODELVIEW );
}

std::vector<Loka3::L3> ExponentTrace( const Loka3::L3 & omega, const Loka3::L3 & fi, const int count )
{
	std::vector<Loka3::L3> curve(count);
	double t = -0.1;
	std::generate_n( begin(curve), count, 
		[&t,omega,fi]()->Loka3::L3 { t+=0.1; return exponent(t*omega+fi); } );
	return curve;
}

std::vector<Loka4::L4> ExponentTrace( const Loka4::L4 & omega, const Loka4::L4 & fi, const int count )
{
	std::vector<Loka4::L4> curve(count);
	double t = -0.1;
	std::generate_n( begin(curve), count, 
		[&t,omega,fi]()->Loka4::L4 { t+=0.1; return exponent(t*omega+fi); } );
	return curve;
}

void DrawExponentTrace( const Loka3::L3 & omega, const Loka3::L3 & fi )
{
	using namespace Loka3;
	using namespace Vector3;

	// Generate curve
	std::vector<L3> curve = ExponentTrace(omega,fi,300);

	// Convert it to vectors
	std::vector<Vect3d> curve_xyz( curve.size() );
	std::transform( begin(curve), end(curve), begin(curve_xyz), 
		[](const L3 & n) { return L3_to_vect(n); } );
	DrawCurve( curve_xyz );

	// Convert it to vectors
	std::vector<Vect3d> curve_abs(curve.size());
	{
		double dt = 1.0 / curve.size();
		double t = -dt;
		std::transform( begin(curve), end(curve), begin(curve_abs), 
			[&t,dt](const L3 & n)->Vect3d { t+=dt; return Vector3::Vect3d(t, norm(n), 0); } );
	}
	DrawCurve( curve_abs );
}

void DrawExponentTrace( const Loka4::L4 & omega, const Loka4::L4 & fi )
{
	using namespace Loka4;
	using namespace Vector3;

	// Generate curve
	std::vector<L4> curve = ExponentTrace(omega,fi,600);

	// Convert it to vectors
	std::vector<Vect3d> curve_xyz( curve.size() );
	std::transform( begin(curve), end(curve), begin(curve_xyz), 
		[](const L4 & n) { return L4_to_vect(n); } );
	DrawCurve( curve_xyz );
}

void DrawMovingConvolution()
{
	using namespace std;
	using namespace Loka3;
	using namespace Vector3;

	// # signal
	vector<double> signal;
	signal.resize(300);
	for( unsigned i = 0; i < 300; ++i )
		signal[i] = 0.3 * (2 + sin(2*2*M_PI/300*i) + 1*sin(20*2*M_PI/300*i));
	vector<Vect3d> signal_xyz(signal.size());
	{
		double dt = 1.0 / signal.size();
		double t = -dt;
		transform( begin(signal), end(signal), begin(signal_xyz), 
			[&t,dt](const double & d)->Vect3d { t+=dt; return Vector3::Vect3d(t, d, 0); } );
	}
	glColor3f(1,0,0);
	DrawCurve( signal_xyz );

	// # exp curv
	vector<L3> expcurve = ExponentTrace( L3(1,0,0), L3(0.7,0,0), 300 );
	// Convert it to vector
	std::vector<Vect3d> expcurve_xyz( expcurve.size() );
	std::transform( begin(expcurve), end(expcurve), begin(expcurve_xyz), 
		[](const L3 & n) { return L3_to_vect(n); } );
	glColor3f(0,1,0);
	DrawCurve( expcurve_xyz );

	// # convolution
	vector<L3> conv;
	for(unsigned i = 0; i < signal.size(); ++i )
	{
		conv.push_back( convolution( signal, expcurve, L3() ) );
		rotate( begin(signal), ++begin(signal), end(signal) );
	}
	// Convert it to vector
	vector<Vect3d> conv_xyz( conv.size() );
	transform( begin(conv), end(conv), begin(conv_xyz), 
		[](const L3 & n) { return L3_to_vect(n)*0.05; } );
	glColor3f(0,0,1);
	DrawCurve( conv_xyz );
}

void displayLoka3()
{
	using namespace std;
	using namespace Loka3;
	using namespace Vector3;

	Vect3d ortX = L3_to_vect( A(1) );
	Vect3d ortY = L3_to_vect( B(1) );
	Vect3d ortZ = L3_to_vect( E(1) );

	// Axes
	glLineWidth(1);
	glBegin( GL_LINES );
	glColor3d(1,0,0);
	glVertex3d(0,0,0);
	glVertex3d(ortX.x,ortX.y,ortX.z);
	glColor3d(0,1,0);
	glVertex3d(0,0,0);
	glVertex3d(ortY.x,ortY.y,ortY.z);
	glColor3d(0,0,1);
	glVertex3d(0,0,0);
	glVertex3d(ortZ.x,ortZ.y,ortZ.z);
	glEnd();

	// Curve
	const L3 omega = L3(1,0,0.0)*4;
	const vector<L3> tbl1 = exponent_table( omega/norm(omega), E(0), 0, exponent_period(), 500 );
	const vector<L3> tbl2 = exponent_table( omega, E(0), 0, exponent_period(), 500 );
	glColor3d(1,1,1);
	DrawCurve( tbl1, true );
	glColor3d(0,1,0);
	DrawCurve( tbl2, true );

	DrawDekartProjectionCurve( tbl1 );
	glColor3d(0,1,1);
	DrawPolarityProjectionCurve( tbl1, 0 );
	glColor3d(1,0,1);
	DrawPolarityProjectionCurve( tbl1, 1 );
	glColor3d(1,1,0);
	DrawPolarityProjectionCurve( tbl1, 2 );

	// Rotate object
	{
		const L3 n_c = vect_to_L3( Vect3d(0,0.3,0) );
		vector<L3> rot_trace( 100 );
		const L3 rot = exponent( L3(1,0,0.5)/40 );
		rot_trace[0] = vect_to_L3( Vect3d(0,0.8,0) );
		for( size_t i = 1; i < rot_trace.size(); ++i )
			rot_trace[i] = rotate( rot_trace[i-1], n_c, rot );
		glColor3d(1,0.5,0.2);
		DrawCurve( rot_trace, true );
	}

	//DrawExponentTrace( normalize(L3(2,0,1)), L3(1,1,0) );
	//glColor3d(1,0,1);
	//DrawExponentTrace( normalize(L3(2,0,1)), L3(0.0,0.1,0.1) );
	//glColor3d(0,1,1);
	//DrawExponentTrace( normalize(L3(0,6,3.1)), );

	//glColor3d(0,1,1);
	//DrawRotationTrace( exponent(L3(0.2,0,0.1)) );

	//DrawMovingConvolution();

}

void displayLoka4()
{
	//glRotated( 2, 1, 0, 0 );
	glRotated( 3, 0, 1, 0.5 );
	//glScaled(0.99,0.99,0.99);

	using namespace std;
	using namespace Loka4;
	using namespace Vector3;

	const Vect3d vA = L4_to_vect( A(1) );
	const Vect3d vB = L4_to_vect( B(1) );
	const Vect3d vC = L4_to_vect( C(1) );
	const Vect3d vE = L4_to_vect( E(1) );

	// Axes
	glLineWidth(1);
	glBegin( GL_LINES );
	glColor3d(1,0,0);
	glVertex3d(0,0,0);
	glVertex3d(vA.x,vA.y,vA.z);
	glColor3d(0,1,0);
	glVertex3d(0,0,0);
	glVertex3d(vB.x,vB.y,vB.z);
	glColor3d(0,0,1);
	glVertex3d(0,0,0);
	glVertex3d(vC.x,vC.y,vC.z);
	glColor3d(1,1,1);
	glVertex3d(0,0,0);
	glVertex3d(vE.x,vE.y,vE.z);
	glEnd();

	/*{
		const L4 omega = L4(0.31,0.01,0.21,0);
		const L4 fi = L4(0,0.5,1,0);
		const vector<L4> tbl = exponent_table( omega/norm(omega), fi, 0, 10, 500 );
		glColor3d(1,0.5,0.5);
		DrawCurve( tbl, false );
	}*/
	//{
	//	// circle
	//	const L4 omega = L4(0.0,0.1,0.0,0.0);
	//	const L4 fi = L4(0,0.5,1,0);
	//	const vector<L4> tbl = exponent_table( omega/norm(omega), fi, 0, 10, 500 );
	//	glColor3d(0.5,1,0.5);
	//	DrawCurve( tbl, false );
	//}
	{
		// circle
		const L4 omega = L4( 0.1, 0.0, 1.0, 0.0 );
		const L4 fi = L4( 0, 0, 1, 0 );
		const vector<L4> tbl = exponent_table( omega/norm(omega), fi, 0, 15, 500 );
		glColor3d(1,0.2,0.2);
		DrawCurve( tbl, false );
	}
	{
		// circle
		const L4 omega = L4( 0.0, 0.1, 0.0, 0.0 );
		const L4 fi = L4( 0.5, 0.5, 0.5, 0 );
		const vector<L4> tbl = exponent_table( omega/norm(omega), fi, 0, 15, 500 );
		glColor3d(0.2,1,0.2);
		DrawCurve( tbl, false );
	}
	{
		// circle
		const L4 omega = L4( 0.0, 0.0, 0.1, 0.0 );
		const L4 fi = L4( 0.5, 0.5, 0.5, 0 );
		const vector<L4> tbl = exponent_table( omega/norm(omega), fi, 0, 15, 500 );
		glColor3d(0.2,0.2,1);
		DrawCurve( tbl, false );
	}
	{
		// circle
		const L4 omega = L4( 0.0, 0.0, 0.0, 0.1 );
		const L4 fi = L4( 0.5, 0.5, 0.5, 0 );
		const vector<L4> tbl = exponent_table( omega/norm(omega), fi, 0, 15, 500 );
		glColor3d(0.2,0.2,0.2);
		DrawCurve( tbl, false );
	}

	//{
	//	// circle
	//	const L4 omega = L4( 0.1, 0.0, 0.0, 0.0 );
	//	const L4 fi = L4( 0.0, 0.0, 0.0, 0.5 );
	//	const vector<L4> tbl = exponent_table( omega/norm(omega), fi, 0, 15, 500 );
	//	glColor3d(1,0.5,0.0);
	//	DrawCurve( tbl, false );
	//}

	/*{
		const L4 omega = L4(0.1,0.0,0.15,0.1);
		const L4 fi = L4(0.2,0.0,1,0);
		const vector<L4> tbl = exponent_table( omega/norm(omega), fi, 0, 10, 500 );
		glColor3d(0.5,0.5,1);
		DrawCurve( tbl, false );
	}*/
}

void CALLBACK display(void)
{
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	//displayLoka3();
	displayLoka4();

	auxSwapBuffers();
}

void ShowOpenGLWindow()
{
	auxInitPosition( 10, 10, 1200, 800);
	auxInitDisplayMode( AUX_RGB | AUX_DEPTH | AUX_DOUBLE );
	auxInitWindow( L"Glaux Template" );
	auxIdleFunc(display);
	auxReshapeFunc(resize);
	/*glEnable(GL_ALPHA_TEST);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	float pos[4] = {3,3,3,1};
	float dir[3] = {-1,-1,-1};    
	glLightfv(GL_LIGHT0, GL_POSITION, pos);
	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir);*/
	auxMainLoop(display);
}

void TestLoka4()
{
	using namespace std;
	using namespace Loka4;

	cout << " ============== Loka4 ================= " << endl;
	L4 n = A(2) + B(3) + C(4) + E(5);
	cout << "n = " << n.showb() << endl;

	L4 m = A(3) + E(2), p = B(3) + E(2), q = C(3) + E(2), r = E(3) + E(2);
    L4 mi = A(11./7.) + B(5./7.) + C(2) + E(1);
	L4 m1 = conjugate1(m);
	L4 m2 = conjugate2(m);
	L4 m3 = conjugate3(m);
	cout << "m*m1*m2*m3 = " << (m*m1*m2*m3).show() << endl;
	cout << "m*p [2a 2b 5c] = " << (m*p).show() << endl;
	cout << "p*q [5a 2b 2c] = " << (p*q).show() << endl;
	cout << "m*q [13 6a 6c] = " << (m*q).show() << endl;
	cout << "m*p*q [4a 13b 4c] = " << (m*p*q).show() << endl;
	cout << "m*p*q*r [...] = " << (m*p*q*r).show() << endl;
	cout << "m*mi [...] = " << (m*mi).show() << endl;

}

int main()
{
	using namespace Loka3;

	std::cout << "Hello world. Loka3 is here." << std::endl;

	std::vector<std::pair<L3,L3>> operands;
	std::vector<L3> results;

	// Multiplication
	operands.push_back( std::pair<L3,L3>( L3(2,0,0), L3(3,0,0) ) );
	operands.push_back( std::pair<L3,L3>( L3(0,2,0), L3(0,3,0) ) );
	operands.push_back( std::pair<L3,L3>( L3(0,0,2), L3(0,0,3) ) );
	operands.push_back( std::pair<L3,L3>( L3(4,0,0), L3(0,4,0) ) );
	operands.push_back( std::pair<L3,L3>( L3(0,4,0), L3(4,0,0) ) );
	operands.push_back( std::pair<L3,L3>( L3(5,0,0), L3(0,0,5) ) );
	operands.push_back( std::pair<L3,L3>( L3(0,0,5), L3(5,0,0) ) );
	operands.push_back( std::pair<L3,L3>( L3(0,6,0), L3(0,0,6) ) );
	operands.push_back( std::pair<L3,L3>( L3(0,0,6), L3(0,6,0) ) );
	operands.push_back( std::pair<L3,L3>( L3(2,3,4), L3(5,6,7) ) );
	for( auto it = operands.begin(); it != operands.end(); ++it )
	{
		results.push_back( it->first * it->second );
	}
	for( unsigned int i = 0; i < operands.size(); ++i )
	{
		std::cout << operands[i].first.showb() <<" * " << operands[i].second.showb() << " = " << results[i].show() << std::endl;
	}

	// Addition
	results.clear();
	for( auto it = operands.begin(); it != operands.end(); ++it )
	{
		results.push_back( it->first + it->second );
	}
	for( unsigned int i = 0; i < operands.size(); ++i )
	{
		std::cout << operands[i].first.showb() <<" + " << operands[i].second.showb() << " = " << results[i].show() << std::endl;
	}

	// Compensation
	results.clear();
	for( auto it = operands.begin(); it != operands.end(); ++it )
	{
		results.push_back( it->first );
	}
	for( unsigned int i = 0; i < operands.size(); ++i )
	{
		std::cout << operands[i].first.show() << " = " << results[i].show() << std::endl;
	}

	// Abs
	L3 n1(1,2,3);
	std::cout << "n*conjugate1*conjugate2 = " << (n1*conjugate1(n1)*conjugate2(n1)).show() << std::endl;
	std::cout << "norm3 = " << norm3(n1) << std::endl;

	// Exp
	{
		L3 n2 = L3(2,0,1);
		L3 n1 = L3(0,2,1);
		L3 n0 = L3(0,0,3);
		for( double t = 0; t < 1; t += 0.01 )
		{
			L3 exp2 = exponent( t * n2 );
			std::cout << "exponent = " << exp2.show() << " | " << norm3(exp2) << std::endl;
		}

		std::cout << "exponent(n0) = " << exponent(n0).show() << std::endl;
		std::cout << "exponent(n1) = " << exponent(n1).show() << std::endl;
		std::cout << "exponent(n2) = " << exponent(n2).show() << std::endl;
	}

	//
	{
		//L3 a = A(2), b = B(3), e = E(4);
		L3 a = A(2)+B(1), b = B(3)+E(5), e = E(4)+A(2);
		std::cout << "Abs(a*b*e) = " << norm3(a*b*e) << std::endl;
		std::cout << "Abs(a)*Abs(b)*Abs(e) = " << norm3(a)*norm3(b)*norm3(e) << std::endl;
		std::cout << "Absc(a*b*e) = " << norm3c(a*b*e) << std::endl;
		std::cout << "Absc(a)*Absc(b)*Absc(e) = " << norm3c(a)*norm3c(b)*norm3c(e) << std::endl;
	}

	{
		L3 a = A(2)+B(1), b = B(3)+E(5), e = E(4)+A(2);
		std::cout << "Exp(a+b) = Exp(a)*Exp(b) " << near_equal(exponent(a+b),exponent(a)*exponent(b)) << std::endl;
		std::cout << "Exp(a+e) = Exp(a)*Exp(e) " << near_equal(exponent(a+e),exponent(a)*exponent(e)) << std::endl;
		std::cout << "Exp(b+e) = Exp(b)*Exp(e) " << near_equal(exponent(b+e),exponent(b)*exponent(e)) << std::endl;
	}
	//std::cout << "Abs(normalize(q1)) = " << norm3(normalize(q1)*normalize(q1)) << std::endl;

	{
		L3 a = A(2)+B(1), b = B(3)+E(5), e = E(4)+A(2);
		std::cout << "Opp(2a+b) = " << opp(a).show() << std::endl;
		std::cout << "Opp(3b+5e) = " << opp(b).show() << std::endl;
		std::cout << "Opp(4e+2a) = " << opp(e).show() << std::endl;
	}

	{ // Convertion Vector3d <-> L3 number
		std::cout << "vect_to_L3( L3_to_vect(A(0)) ) = " << vect_to_L3( L3_to_vect(A(0)) ).show() << std::endl;
		std::cout << "vect_to_L3( L3_to_vect(A(2)) ) = " << vect_to_L3( L3_to_vect(A(2)) ).show() << std::endl;
		std::cout << "vect_to_L3( L3_to_vect(B(2)) ) = " << vect_to_L3( L3_to_vect(B(2)) ).show() << std::endl;
		std::cout << "vect_to_L3( L3_to_vect(E(2)) ) = " << vect_to_L3( L3_to_vect(E(2)) ).show() << std::endl;
		std::cout << "vect_to_L3( L3_to_vect(A(2)+B(3)) ) = " << vect_to_L3( L3_to_vect(A(2)+B(3)) ).show() << std::endl;
		std::cout << "vect_to_L3( L3_to_vect(B(2)+E(4)) ) = " << vect_to_L3( L3_to_vect(B(2)+E(4)) ).show() << std::endl;
		std::cout << "vect_to_L3( L3_to_vect(E(2)+A(1)) ) = " << vect_to_L3( L3_to_vect(E(2)+A(1)) ).show() << std::endl;

	}


	TestLoka4();

	ShowOpenGLWindow();
	//	int i;
	//	std::cin >> i;
	return 0;
}


