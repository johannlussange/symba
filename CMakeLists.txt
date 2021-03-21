cmake_minimum_required(VERSION 3.17.0)
project(symba VERSION 0.1.0)

option(SYMBA_LINK_STATIC "Link Symba libraries statically" ON)

find_package(GSL REQUIRED)
find_package(CLI11 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

add_executable(symba)
file(GLOB_RECURSE symbaSources 
    LIST_DIRECTORIES true 
    CONFIGURE_DEPENDS
    "src/*.cpp"
)
target_sources(symba PRIVATE ${symbaSources})
target_include_directories(symba PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_features(symba PRIVATE cxx_std_20)
set_target_properties(symba PROPERTIES
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>"
)
target_link_libraries(symba PRIVATE -static
    GSL::gsl
    GSL::gslcblas
    CLI11::CLI11
    nlohmann_json::nlohmann_json
)

install(TARGETS symba)