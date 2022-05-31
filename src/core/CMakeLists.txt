project(core)

include_directories(${HOME_DIR}/include)

file(GLOB SRC_LIST src/*.cpp)

set(LIBRARY_OUTPUT_PATH ${HOME_DIR}/lib)

link_directories(${HOME_DIR}/lib /usr/lib64)

link_libraries(gadget time)

add_library(core SHARED ${SRC_LIST})