cmake_minimum_required(VERSION 3.0)

project(RGB_EPP)
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

set(CMAKE_CXX_STANDARD 17)

add_executable(splitfasta splitfasta.cpp)
add_executable(sortdiamond sortdiamond.cpp)
	
target_compile_options(splitfasta PRIVATE -Wall -Wextra -pedantic -O3)
target_compile_options(sortdiamond PRIVATE -Wall -Wextra -pedantic -O3)

install(FILES RGBEPP.sh DESTINATION /usr/bin/RGBEPP.sh)
install(FILES splitfasta DESTINATION /usr/bin/splitfasta)
install(FILES sortdiamond DESTINATION /usr/bin/sortdiamond)
