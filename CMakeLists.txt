cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
get_filename_component(name ${CMAKE_SOURCE_DIR} NAME_WE)
project(${name})

add_compile_options(-fPIC)

find_package(Geant4 REQUIRED)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# QT: 5
add_compile_options(-march=native -Ofast)
#add_compile_options(-O0)
add_compile_options(-ggdb3)

#set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fsanitize=address" )
#set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address" )

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


find_package(Qt5Core REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5Gui REQUIRED)
# TODO apply to all headers: 
qt5_wrap_cpp(WSRCS ${headers})
include_directories(${Qt5Widgets_INCLUDE_DIRS})
include_directories(${Qt5Gui_INCLUDE_DIRS})
include_directories(${Qt5Core_INCLUDE_DIRS})





add_library(_base OBJECT ${sources} ${headers} ${WSRCS})


# Build and install all targets
file(GLOB mainfiles ${PROJECT_SOURCE_DIR}/*.cc)
foreach(mainf ${mainfiles})
    get_filename_component(main ${mainf} NAME_WE)
    add_executable(${main} ${mainf} $<TARGET_OBJECTS:_base>)
    target_link_libraries(${main} ${Geant4_LIBRARIES})
    install(TARGETS ${main} DESTINATION bin)
endforeach()

