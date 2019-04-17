#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>

#include "../Loka3Lib/Loka3.h" 
#include "../Loka3Lib/Loka3math.h"

int main()
{
    using namespace Loka3;
    
    L3 n1(1,2,0);
    L3 n2(1,0,2);
    std::cout << n1.showb() << " * " << n2.showb() << " = " << (n1 * n2).showb() << std::endl;
    
    std::cout << "pow_a(2) = " << pow_a(2).showb() << std::endl;
    std::cout << "pow_b(2) = " << pow_b(2).showb() << std::endl;
    
    std::cout << "pow_a(2)*pow_b(2) = " << (pow_a(2)*pow_b(2)).showb() << std::endl;
    std::cout << "ratio(2,2,2) = " << ratio(2,2,2).showb() << std::endl;
    std::cout << "ratio(2,3,4)*ratio(5,7,3) = " << (ratio(2,3,4)*ratio(5,7,3)).showb() << std::endl;
    std::cout << "ratio(10,21,12) = " << ratio(10,21,12).showb() << std::endl;
    
    
    std::vector<L3> trajectory;
    for( int i = 1; i < 50; ++i )
        for( int j = 1; j < 50; ++j )
            for( int k = 1; k < 50; ++k )
                trajectory.push_back(ratio(i,j,k));
    
    {
        std::ofstream f("traj.txt");
        for( const auto & v : trajectory )
            f << v.e() << ", " << v.a() << ", " << v.b() << "\n"; 
    }
    
    return 0;
}
