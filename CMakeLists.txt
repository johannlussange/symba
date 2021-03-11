cmake_minimum_required(VERSION 3.17.0)
project(symba VERSION 0.1.0)

find_package(GSL REQUIRED)
find_package(CLI11 CONFIG REQUIRED)

add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE
	"src/symba.cpp"
	"src/order_book.cpp"
	"src/agent.cpp"
	"src/utility.cpp"
	"src/random_walk.cpp"
	"src/share.cpp"
)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)
set_target_properties(${PROJECT_NAME} PROPERTIES
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>"
)
target_link_libraries(${PROJECT_NAME} PRIVATE GSL::gsl GSL::gslcblas CLI11::CLI11)

install(TARGETS ${PROJECT_NAME})