#g++ -I include -I /mnt/crick/pelina/bin/boost_1_54_0 
#     -L /home/rickards/testprojects/rnaseqpipeline/programs/bamtools/lib 
#     -I /home/rickards/testprojects/rnaseqpipeline/programs/bamtools/include 
#     -std=c++0x 
#     -fopenmp 
#     ProcessBAMFiles.cpp alglibinternal.cpp alglibmisc.cpp ap.cpp dataanalysis.cpp diffequations.cpp fasttransforms.cpp integration.cpp interpolation.cpp linalg.cpp optimization.cpp solvers.cpp specialfunctions.cpp statistics.cpp 
#     -o HiCapAnalysis_perChr 
#     -lz 
#     -lbamtools

#First steps
#mkdir /tmp/build
#cd /tmp/build
#cmake '-DCMAKE_PREFIX_PATH=/mnt/crick/pelina/bin/HiCapAnalysis-readfrombam/include;/mnt/crick/pelina/bin/boost_1_54_0;/home/rickards/testprojects/rnaseqpipeline/programs/bamtools;/home/rickards/testprojects/rnaseqpipeline/programs/bamtools/include/' '-DCMAKE_INSTALL_PATH=/mnt/kauffman/pelina/bin/HiCapAnalysis-readfrombam/src'
#make
#make install

cmake_minimum_required (VERSION 3.1)

project(AnalyseHiCap_perChr CXX)

find_package(OpenMP REQUIRED)
include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

  # set our library and executable destination dirs
  set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )
  
   add_executable(AnalyseHiCap_perChr ProcessBAMFiles.cpp alglibinternal.cpp alglibmisc.cpp ap.cpp dataanalysis.cpp diffequations.cpp fasttransforms.cpp integration.cpp interpolation.cpp linalg.cpp optimization.cpp solvers.cpp specialfunctions.cpp statistics.cpp)
   
  set_property(TARGET AnalyseHiCap_perChr PROPERTY CXX_STANDARD 11)
  set_property(TARGET AnalyseHiCap_perChr PROPERTY CXX_STANDARD_REQUIRED ON)

  find_package(ZLIB REQUIRED)
  add_library(libbamtools STATIC IMPORTED)
  set_target_properties(libbamtools PROPERTIES IMPORTED_LOCATION /sw/apps/bioinfo/bamtools/2.3.0/milou/lib/libbamtools.so.2.3.0)	
  include_directories(/sw/apps/bioinfo/bamtools/2.3.0/milou/include)
  include_directories(/pica/h1/pelin/External/boost_1_58_0)
  link_directories(/sw/apps/bioinfo/bamtools/2.3.0/milou/lib)
  target_link_libraries( AnalyseHiCap_perChr libbamtools ${ZLIB_LIBRARIES})

