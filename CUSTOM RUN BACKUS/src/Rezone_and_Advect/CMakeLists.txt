include_directories($ENV{JSONFORTRAN})
link_directories($ENV{JSONFORTRAN}/lib)

add_library(Rezone_and_Advect advect.f90 rezone.f90)
target_link_libraries(Rezone_and_Advect jsonfortran)
