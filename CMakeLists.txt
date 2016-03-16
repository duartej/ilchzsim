cmake_minimum_required(VERSION 2.8)
project(ilchzsim)
set(CMAKE_INSTALL_PREFIX $ENV{HOME}/.local)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC")

add_definitions(-std=c++11)

add_executable(ilchz ilchz.cc)

# ROOT dependencies
set(ROOTLIBSTOLINK Tree RIO)
find_package(ROOT REQUIRED ${ROOTLIBSTOLINK})
include_directories(${ROOT_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ROOT_LIBRARIES})

find_package(Pythia8 REQUIRED)
include_directories(${PYTHIA8_INCLUDE_DIRS})
set(LIBS ${LIBS} ${PYTHIA8_LIBRARIES})

target_link_libraries(ilchz ${LIBS})

install(TARGETS ilchz 
    RUNTIME DESTINATION bin)
install(PROGRAMS utils/processhzroot.py 
    DESTINATION bin
    RENAME processhzroot)
install(PROGRAMS utils/hzplots.py 
    DESTINATION bin
    RENAME hzplots)
install(PROGRAMS utils/kinplots.py 
    DESTINATION bin
    RENAME kinplots)

