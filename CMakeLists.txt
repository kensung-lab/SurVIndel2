cmake_minimum_required(VERSION 3.13)
project(SurVIndel2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

include_directories(htslib-1.17)
include_directories(htslib-1.17/include)

set(READS_CATEGORIZER_SRC reads_categorizer.cpp)
add_executable(reads_categorizer ${READS_CATEGORIZER_SRC})

set(CLIP_CONSENSUS_BUILDER_SRC libs/ssw.c libs/ssw_cpp.cpp clip_consensus_builder.cpp)
add_executable(clip_consensus_builder ${CLIP_CONSENSUS_BUILDER_SRC})

set(MERGE_IDENTICAL_CALLS merge_identical_calls.cpp)
add_executable(merge_identical_calls ${MERGE_IDENTICAL_CALLS})

set(DP_CLUSTERER_SRC libs/ssw.c libs/ssw_cpp.cpp dp_clusterer.cpp)
add_executable(dp_clusterer ${DP_CLUSTERER_SRC})

set(NORMALISE_INDELS_SRC normalise.cpp)
add_executable(normalise ${NORMALISE_INDELS_SRC})

find_library(HTS_LIB hts PATHS htslib-1.17/lib NO_DEFAULT_PATH)
target_link_libraries(reads_categorizer "${HTS_LIB}")
target_link_libraries(clip_consensus_builder "${HTS_LIB}")
target_link_libraries(merge_identical_calls "${HTS_LIB}")
target_link_libraries(dp_clusterer "${HTS_LIB}")
target_link_libraries(normalise "${HTS_LIB}")
