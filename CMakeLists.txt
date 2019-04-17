cmake_minimum_required(VERSION 2.8)

project(loka)
add_executable(${PROJECT_NAME} 
	"loka/main.cpp" 
	"Loka3Lib/Loka3.h" 
	"Loka3Lib/Loka3.cpp" 
	"Loka3Lib/Loka3math.h" 
	"Loka3Lib/Loka3math.cpp"
	"Loka4Lib/Loka4.h" 
	"Loka4Lib/Loka4.cpp" 
	"Loka4Lib/Loka4math.h" 
	"Loka4Lib/Loka4math.cpp")

