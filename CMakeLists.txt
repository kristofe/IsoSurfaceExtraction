cmake_minimum_required  (VERSION 2.8)

set(app_name quadratic_iso)
project(${app_name})
include(GenerateExportHeader)

FILE(GLOB MAIN_SOURCE Src/*.*)

set(SOURCES
   ${MAIN_SOURCE}
)

set(HEADERS
)

set(RESOURCES
)

set(ASSETS
   )

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
option(CMAKE_BUILD_TYPE  "DEBUG or RELEASE" RELEASE)
string( TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_LOWER )
if(BUILD_TYPE_LOWER STREQUAL "debug")
   set(DEBUG 1)
endif()

find_package(Threads)
list(APPEND EXTRA_LIBS ${CMAKE_THREAD_LIBS_INIT} )

###############################################################################
###############################################################################
# -g adds debugging symbols
# -Wall turns on all warnings
# -Wextra turns on a lot of warnings (but not too pedantic)
if(APPLE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wextra -Wno-ignored-qualifiers -Wno-unused-parameter -std=c++11")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -g -Wextra -Wno-ignored-qualifiers -Wno-unused-parameter")
  add_definitions(-DTARGETMAC)
ENDIF(APPLE)

if(WIN32)
  #LIST(APPEND SOURCES source/PlatformClasses/Windows/Platform_Win32.cpp)
  #LIST(APPEND HEADERS source/PlatformClasses/Windows/Platform_Win32.h)
  #include_directories(source/PlatformClasses/Windows)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  add_definitions(-D_CRT_SECURE_NO_WARNINGS) #disable sprintf deprecation
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS 1)
ENDIF(WIN32)


if(BUILD_TYPE_LOWER STREQUAL "debug")
    message("adding debug flags")
    add_definitions(-DDEBUG)
    add_definitions(-D_DEBUG)
    add_definitions(-DBREAK_ON_EXCEPTION)
endif(BUILD_TYPE_LOWER STREQUAL "debug")

include_directories(source)
include_directories(source/PlatformClasses)

if(APPLE)
   FIND_LIBRARY(IOKIT_LIBRARY IOKit)
   FIND_LIBRARY(COCOA_LIBRARY Cocoa)
   MARK_AS_ADVANCED(
                     IOKIT_LIBRARY
                     COCOA_LIBRARY
                   )
  LIST(
      APPEND
      EXTRA_LIBS
      ${IOKIT_LIBRARY}
      ${COCOA_LIBRARY}
      )
ENDIF(APPLE)


if(UNIX AND NOT APPLE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wextra -Wno-ignored-qualifiers -Wno-unused-parameter -std=c++11")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -g -Wextra -Wno-ignored-qualifiers -Wno-unused-parameter")
  add_definitions(-DTARGEUNIX)
  #LIST(APPEND SOURCES source/PlatformClasses/Unix/Platform_Unix.cpp)
  #LIST(APPEND HEADERS source/PlatformClasses/Unix/Platform_Unix.h)
  #include_directories(source/PlatformClasses/Unix)
ENDIF(UNIX AND NOT APPLE)

if(BUILD_TYPE_LOWER STREQUAL "debug")
  message("Creating Debug Build")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -O0")
else()
  message("Creating Release Build")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -O3")
endif()


#message("EXTRA LIBS " ${EXTRA_LIBS})

# EXAMPLE CALL (Xcode)
#cmake .. -G"Xcode"
#Then build the quadratic_iso target

# EXAMPLE CALL (Unix Makefiles)
#cmake .. -G"Unix Makefiles"
#Then call: make   -j8
##############################################################################
##############################################################################
#  WARNING!!!!!!!!!!!!!!!!!!!!
#  WARNING!!!!!!!!!!!!!!!!!!!!
#  WARNING!!!!!!!!!!!!!!!!!!!!
#  On unix you have to create a build/source directory or it won't compile
##############################################################################
##############################################################################

# TODO:  Remove the "d" from _isosurface.so (only happens in Xcode)

add_library(${app_name} SHARED ${SOURCES})

target_link_libraries(${app_name} ${EXTRA_LIBS})

