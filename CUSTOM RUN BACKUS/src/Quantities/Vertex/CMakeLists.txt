add_library(Q_Vertex vertex_quantity.f90 acceleration.f90 velocity.f90 momentum.f90 vertex_mass.f90 coordinates.f90 )
target_link_libraries(Q_Vertex Quantities)
