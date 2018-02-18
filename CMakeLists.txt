cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
get_filename_component(name ${CMAKE_SOURCE_DIR} NAME_WE)
project(${name})

add_compile_options(-fPIC)

# Geant4 overwrites flags :-(
set(M_CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}")
set(M_CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}")
set(M_CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL}")
set(M_CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}")
set(M_CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(M_CMAKE_EXE_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS}")
find_package(Geant4 REQUIRED)

find_package(Qt5Core REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5Gui REQUIRED)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# add_compile_options(-march=native -Ofast)
add_compile_options(-ggdb3)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

qt5_wrap_cpp(WSRCS ${headers})
include_directories(${Qt5Widgets_INCLUDE_DIRS})
include_directories(${Qt5Gui_INCLUDE_DIRS})
include_directories(${Qt5Core_INCLUDE_DIRS})

add_library(_base OBJECT ${sources} ${headers} ${WSRCS})

file(GLOB mainfiles ${PROJECT_SOURCE_DIR}/*.cc)
foreach(mainf ${mainfiles})
    get_filename_component(main ${mainf} NAME_WE)
    add_executable(${main} ${mainf} $<TARGET_OBJECTS:_base>)
    target_link_libraries(${main} ${Geant4_LIBRARIES} Qt5::Widgets Qt5::Core Qt5::Gui)
    install(TARGETS ${main} DESTINATION bin)
endforeach()

# Geant4 overwrites flags, so we add back our own
set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} ${M_CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} ${M_CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} ${M_CMAKE_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} ${M_CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${M_CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_EXE_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS} ${M_CMAKE_EXE_LINKER_FLAGS}")

