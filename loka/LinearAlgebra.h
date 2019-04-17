#pragma once

#include <numeric>

template< typename TL, typename TR, typename T0 >
T0 convolution( const TL & al, const TR & ar, const T0 & s0 )
{
	return inner_product( begin(al), end(al), begin(ar), s0 );
}