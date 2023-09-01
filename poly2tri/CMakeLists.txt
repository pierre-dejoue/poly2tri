#
# Poly2tri library
#
file(GLOB SOURCES src/common/*.cc src/sweep/*.cc)
file(GLOB PUBLIC_HEADERS include/*.h include/common/*.h include/sweep/*.h)
file(GLOB PRIVATE_HEADERS src/common/*.h src/sweep/*.h)

add_library(poly2tri ${SOURCES} ${PUBLIC_HEADERS} ${PRIVATE_HEADERS})

target_include_directories(poly2tri
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

get_target_property(poly2tri_target_type poly2tri TYPE)
if(poly2tri_target_type STREQUAL SHARED_LIBRARY)
  target_compile_definitions(poly2tri PRIVATE P2T_SHARED_EXPORTS)
else()
  target_compile_definitions(poly2tri PUBLIC P2T_STATIC_EXPORTS)
endif()
