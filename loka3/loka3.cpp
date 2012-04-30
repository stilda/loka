// loka3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <Windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "glaux.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "Loka3.h"
#include "Loka4.h"
#include "Vector3.h"
#include "LinearAlgebra.h"

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

void DrawCurve( const std::vector<Vector3::Vect3> & curve_xyz )
{
	using namespace Vector3;

	// Curve
	glLineWidth( 2 );
	glBegin( GL_LINE_STRIP );
	//glBegin( GL_POINTS );
	std::for_each( begin(curve_xyz), end(curve_xyz), 
		[](const Vect3 & v){ glVertex3d( v.x, v.y, v.z ); } );
	glEnd();

	// Curve start
	glLineWidth(1);
	glColor3d(1,1,0);
	DrawCircle2D( curve_xyz[0].x, curve_xyz[0].y, 0.01 );
	
	// Curve direction
	Vect3 dir = curve_xyz[1] - curve_xyz[0];
	dir = dir * (0.2 / abs(dir));
	Vect3 ptEnd = curve_xyz[0] + dir;
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
		//std::cout << abs3(acc) << std::endl;
	}

	// Convert it to vectors
	std::vector<Vect3> curve_xyz;
	for( auto it = curve.begin(); it != curve.end(); ++it )
	{
		curve_xyz.push_back( tovect(*it) );
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
	std::vector<Vect3> curve_xyz( curve.size() );
	std::transform( begin(curve), end(curve), begin(curve_xyz), 
		[](const L3 & n) { return tovect(n); } );
	DrawCurve( curve_xyz );

	// Convert it to vectors
	std::vector<Vect3> curve_abs(curve.size());
	{
		double dt = 1.0 / curve.size();
		double t = -dt;
		std::transform( begin(curve), end(curve), begin(curve_abs), 
			[&t,dt](const L3 & n)->Vect3 { t+=dt; return Vector3::Vect3(t, abs(n), 0); } );
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
	std::vector<Vect3> curve_xyz( curve.size() );
	std::transform( begin(curve), end(curve), begin(curve_xyz), 
		[](const L4 & n) { return tovect(n); } );
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
		signal[i] =  sin(5*2*M_PI/300*i) + 0.2*sin(20*2*M_PI/300*i);
	vector<Vect3> signal_xyz(signal.size());
	{
		double dt = 1.0 / signal.size();
		double t = -dt;
		transform( begin(signal), end(signal), begin(signal_xyz), 
			[&t,dt](const double & d)->Vect3 { t+=dt; return Vector3::Vect3(t, d, 0); } );
	}
	glColor3f(1,0,0);
	DrawCurve( signal_xyz );

	// # exp curv
	vector<L3> expcurve = ExponentTrace( L3(1,0,0), L3(0.7,0,0), 300 );
	// Convert it to vector
	std::vector<Vect3> expcurve_xyz( expcurve.size() );
	std::transform( begin(expcurve), end(expcurve), begin(expcurve_xyz), 
		[](const L3 & n) { return tovect(n); } );
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
	vector<Vect3> conv_xyz( conv.size() );
	transform( begin(conv), end(conv), begin(conv_xyz), 
		[](const L3 & n) { return tovect(n)*0.05; } );
	glColor3f(0,0,1);
	DrawCurve( conv_xyz );
}

void displayLoka3()
{
	using namespace Loka3;
	using namespace Vector3;
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
	//glColor3d(1,1,1);
	//DrawExponentTrace( normalize(L3(2,0,1)), L3(0,0,0) );
	//glColor3d(1,0,1);
	//DrawExponentTrace( normalize(L3(0.2,0.1,0.1)), L3(0.0,0.1,0.1) );
	//glColor3d(0,1,1);
	//DrawExponentTrace( normalize(L3(0,6,3.1)), );
	
	//glColor3d(0,1,1);
	//DrawRotationTrace( exponent(L3(0.2,0,0.1)) );

	DrawMovingConvolution();

}

void displayLoka4()
{
	glRotated( 2, 1, 0, 0 );
	glRotated( 3, 0, 1, 0.5 );
	glScaled(0.99,0.99,0.99);
	

	using namespace Loka4;
	using namespace Vector3;

	const Vect3 vA = tovect( A(1) );
	const Vect3 vB = tovect( B(1) );
	const Vect3 vC = tovect( C(1) );
	const Vect3 vE = tovect( E(1) );

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

	// Curve
	glColor3d(1,0.5,0.5);
	DrawExponentTrace( L4(0.31,0.01,0.21,0), L4(0,0.5,1,0) );
	glColor3d(0.5,1,0.5);
	DrawExponentTrace( L4(0.0,0.1,0.2,0.1), L4(0,0.5,1,0) );
	glColor3d(0.5,0.5,1);
	DrawExponentTrace( L4(0.0,0.1,0.2,0.1), L4(0,0.5,1,0) );
	//glColor3d(1,0,1);
	//DrawExponentTrace( normalize(L3(0.2,0.1,0.1)), L3(0.0,0.1,0.1) );
	//glColor3d(0,1,1);
	//DrawExponentTrace( normalize(L3(0,6,3.1)), );
	
	//glColor3d(0,1,1);
	//DrawRotationTrace( exponent(L3(0.2,0,0.1)) );

}

void CALLBACK display(void)
{
	displayLoka3();

	// Draw
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
	L4 m1 = conjugate1(m);
	L4 m2 = conjugate2(m);
	L4 m3 = conjugate3(m);
	cout << "m*m1*m2*m3 = " << (m*m1*m2*m3).show() << endl;
	cout << "m*p [2a 2b 5c] = " << (m*p).show() << endl;
	cout << "p*q [5a 2b 2c] = " << (p*q).show() << endl;
	cout << "m*q [13 6a 6c] = " << (m*q).show() << endl;
	cout << "m*p*q [4a 13b 4c] = " << (m*p*q).show() << endl;
	cout << "m*p*q*r [...] = " << (m*p*q*r).show() << endl;

}

int _tmain(int argc, _TCHAR* argv[])
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
	std::cout << "abs3 = " << abs3(n1) << std::endl;

	// Exp
	{
		L3 n2 = L3(2,0,1);
		L3 n1 = L3(0,2,1);
		L3 n0 = L3(0,0,3);
		for( double t = 0; t < 1; t += 0.01 )
		{
			L3 exp2 = exponent( t * n2 );
			std::cout << "exponent = " << exp2.show() << " | " << abs3(exp2) << std::endl;
		}
		
		std::cout << "exponent(n0) = " << exponent(n0).show() << std::endl;
		std::cout << "exponent(n1) = " << exponent(n1).show() << std::endl;
		std::cout << "exponent(n2) = " << exponent(n2).show() << std::endl;
	}

	//
	{
		//L3 a = A(2), b = B(3), e = E(4);
		L3 a = A(2)+B(1), b = B(3)+E(5), e = E(4)+A(2);
		std::cout << "Abs(a*b*e) = " << abs3(a*b*e) << std::endl;
		std::cout << "Abs(a)*Abs(b)*Abs(e) = " << abs3(a)*abs3(b)*abs3(e) << std::endl;
		std::cout << "Absc(a*b*e) = " << abs3c(a*b*e) << std::endl;
		std::cout << "Absc(a)*Absc(b)*Absc(e) = " << abs3c(a)*abs3c(b)*abs3c(e) << std::endl;
	}
	
	{
		L3 a = A(2)+B(1), b = B(3)+E(5), e = E(4)+A(2);
		std::cout << "Exp(a+b) = Exp(a)*Exp(b) " << NearEqual(exponent(a+b),exponent(a)*exponent(b)) << std::endl;
		std::cout << "Exp(a+e) = Exp(a)*Exp(e) " << NearEqual(exponent(a+e),exponent(a)*exponent(e)) << std::endl;
		std::cout << "Exp(b+e) = Exp(b)*Exp(e) " << NearEqual(exponent(b+e),exponent(b)*exponent(e)) << std::endl;
	}
	//std::cout << "Abs(normalize(q1)) = " << abs3(normalize(q1)*normalize(q1)) << std::endl;

	{
		using namespace std;
		using namespace Loka3;

		vector<double> signal;
		signal.resize(300);
		for( unsigned i = 0; i < 300; ++i )
			signal[i] =  sin(2*M_PI/50*i);

		vector<L3> expcurve = ExponentTrace( normalize(L3(0.1,0,0)), L3(0,0,0), 300 );

		L3 conv = convolution( signal, expcurve, L3() );

		std::cout << "conv = " << conv.show() << std::endl;
	}

	TestLoka4();

	ShowOpenGLWindow();
//	int i;
//	std::cin >> i;
	return 0;
}


