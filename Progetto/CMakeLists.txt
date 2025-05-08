cmake_minimum_required(VERSION 3.9)

project(ProgettoPCS2025 LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

# Add dependencies
################################################################################

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND ProgettoPCS2025_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads
find_package(Threads REQUIRED)
list(APPEND ProgettoPCS2025_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest
find_package(GTest REQUIRED)

list(APPEND ProgettoPCS2025_LINKED_LIBRARIES PRIVATE GTest::gtest)
list(APPEND ProgettoPCS2025_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
list(APPEND ProgettoPCS2025_LINKED_LIBRARIES PRIVATE GTest::gmock)
list(APPEND ProgettoPCS2025_LINKED_LIBRARIES PRIVATE GTest::gmock_main)

# SubDirectory
################################################################################

add_subdirectory(src)

list(APPEND ProgettoPCS2025_sources ${polyhedralMesh_sources})
list(APPEND ProgettoPCS2025_headers ${polyhedralMesh_headers})
list(APPEND ProgettoPCS2025_includes ${polyhedralMesh_includes})

add_subdirectory(ExportParaview)

list(APPEND ProgettoPCS2025_sources ${src_paraview_sources})
list(APPEND ProgettoPCS2025_headers ${src_paraview_headers})
list(APPEND ProgettoPCS2025_includes ${src_paraview_includes})

# Create executable
################################################################################

add_executable(${PROJECT_NAME} 
	main.cpp 
	${${PROJECT_NAME}_sources}
	${${PROJECT_NAME}_headers})


#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedralMesh/Cell0Ds.csv  ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedralMesh/Cell1Ds.csv  ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedralMesh/Cell2Ds.csv
      #DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(${PROJECT_NAME} ${ProgettoPCS2025_LINKED_LIBRARIES})
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${${PROJECT_NAME}_includes})