cmake_minimum_required(VERSION 3.9.5)
 
project(spine_ext)
 
FIND_PACKAGE(ITK)

IF(ITK_FOUND)
    
  INCLUDE(${ITK_USE_FILE})

ELSE(ITK_FOUND)
    
  MESSAGE(FATAL_ERROR
    "ITK not found. Please set ITK_DIR.")

ENDIF(ITK_FOUND)
 
add_executable(spine_ext spine_ext.cxx)
target_link_libraries(spine_ext ${ITK_LIBRARIES})