include_directories($ENV{JSONFORTRAN})
link_directories($ENV{JSONFORTRAN}/lib)
add_library(Mesh  mesh_base.f90 mesh.f90 mesh_3d.f90)
target_link_libraries(Mesh jsonfortran)

