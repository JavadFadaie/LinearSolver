cmake_minimum_required(VERSION 3.12)

project(runTest)

#Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Link runTests with what we want to test and the GTest and pthread library
add_executable(${PROJECT_NAME} VectorTest.cpp CSRMatricTest.cpp)


target_link_libraries(${PROJECT_NAME} ${GTEST_LIBRARIES} pthread)


