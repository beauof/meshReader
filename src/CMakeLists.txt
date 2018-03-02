add_library(meshReader meshReader.F90)

target_include_directories(meshReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS meshReader DESTINATION ${PROJECT_SOURCE_DIR}/lib)


