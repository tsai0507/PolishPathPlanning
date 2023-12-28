cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(PathPlanning)

find_package(PCL 1.10 REQUIRED )
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)

set(CMAKE_C_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
# add_compile_options(-Wall -lX11)
add_compile_options(-lX11)

include_directories(include ${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})


# add_executable(slicing_method src/slicing_method.cpp)
# target_link_libraries(slicing_method  ${PCL_LIBRARIES} GSL::gsl GSL::gslcblas)

add_executable(main src/main.cpp src/Path_Generation.cpp)
target_link_libraries(main ${PCL_LIBRARIES} GSL::gsl GSL::gslcblas)

add_executable(connect src/connect.cpp 
    src/path_slicing_alg.cpp
    src/path_translation_alg.cpp
    src/path_dynamic_alg.cpp 
)
target_link_libraries(connect 
    ${PCL_LIBRARIES} 
    GSL::gsl GSL::gslcblas 
    ${X11_LIBRARIES}
)

add_executable(ft_test_path src/ft_test_path.cpp)
target_link_libraries(ft_test_path Eigen3::Eigen)

add_executable(contour src/contour.cpp src/contour_alg.cpp)
target_link_libraries(contour ${PCL_LIBRARIES} GSL::gsl GSL::gslcblas)
