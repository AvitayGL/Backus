
module slip_cell_3d_module
    use cell_boundary_condition_module, only : cell_boundary_condition_t
    use data_module    , only : data_t
    use quantity_module, only : quantity_t

    implicit none
    private

    type, extends(cell_boundary_condition_t), public :: slip_cell_3d_t
        private
    contains

        procedure, public :: Calculate => Slip_cell_3d_calculate
    end type slip_cell_3d_t

contains
    subroutine Slip_cell_3d_calculate (this, c_quantity, edge_num, is_offload)
        class (slip_cell_3d_t) , intent (in out)     :: this      
        class(quantity_t)   , intent (in out)     :: c_quantity  
        integer             , intent (in)         :: edge_num  
        logical, optional :: is_offload

        real(8), dimension(:,:,:), pointer :: values
        real(8), dimension(:,:,:,:), pointer :: values_4d
        integer :: i, j, k, nxp, nyp, nzp,m, nmats

        logical :: is_offload_local

        if (.not. present(is_offload)) then
        
            !write(*,*) "Slip_cell_3d_calculate is_offload_local NOT PRESENT!"
            is_offload_local = .False.  
        
        else
            is_offload_local = is_offload  
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! yoni for now to make sure:
        !is_offload_local = .True.
        !write(*,*) "Slip_cell_3d_calculate called, is_offload_local=", is_offload_local
        nxp = c_quantity%d1 + 1
        nyp = c_quantity%d2 + 1
        nzp = c_quantity%d3 + 1
        if (associated(c_quantity%data_4d)) then
            call c_quantity%Point_to_data(values_4d)
            nmats = c_quantity%data_4d%nmats
            select case(edge_num)
                case(1)
                    i = 1
                    !$omp target teams distribute parallel do collapse(2) firstprivate(i) if(is_offload_local)
                    do k = 0, nzp
                        do j = 0, nyp
                            do m = 1, nmats
                                values_4d(m, i-1, j, k) = values_4d(m, i, j, k)
                            end do
                        end do
                    end do
                    !$omp end target teams distribute parallel do


                case(2)
                    i = nxp
                    !$omp target teams distribute parallel do collapse(2) firstprivate(i) if(is_offload_local) 
                    do k = 0, nzp
                        do j = 0, nyp
                            do m = 1, nmats
                                values_4d(m, i, j, k) = values_4d(m, i-1, j, k)
                            end do
                        end do
                    end do
                    !$omp end target teams distribute parallel do

                case(3)
                    j = 1
                    !$omp target teams distribute parallel do collapse(2) firstprivate(j) if(is_offload_local)
                    do k = 0, nzp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j-1, k) = values_4d(m, i, j, k)
                            end do
                        end do
                    end do
                    !$omp end target teams distribute parallel do

                case(4)               
                    j = nyp
                    !$omp target teams distribute parallel do collapse(2) firstprivate(j) if(is_offload_local)
                    do k = 0, nzp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j, k) = values_4d(m, i, j-1, k)
                            end do
                        end do
                    end do
                    !$omp end target teams distribute parallel do

                case(5)
                    k = 1
                    !$omp target teams distribute parallel do collapse(2) firstprivate(k) if(is_offload_local)
                    do j = 0, nyp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j, k-1) = values_4d(m, i, j, k)
                            end do
                        end do
                    end do
                    !$omp end target teams distribute parallel do

                case(6)
                    k = nzp
                    !$omp target teams distribute parallel do collapse(2) firstprivate(k) if(is_offload_local)
                    do j = 0, nyp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j, k) = values_4d(m, i, j, k-1)
                            end do
                        end do
                    end do
                    !$omp end target teams distribute parallel do
            end select
        else
            call c_quantity%Point_to_data(values)

            select case(edge_num)
                case(1)
                    i = 1
                    !$omp target teams distribute parallel do collapse(2) firstprivate(i) if(is_offload_local)
                    do k = 0, nzp
                        do j = 0, nyp
                            values(i-1, j, k) = values(i, j, k)
                        end do
                    end do
                    !$omp end target teams distribute parallel do

                case(2)
                    i = nxp
                    !$omp target teams distribute parallel do collapse(2) firstprivate(i) if(is_offload_local)
                    do k = 0, nzp
                        do j = 0, nyp
                            values(i, j, k) = values(i-1, j, k)
                        end do
                    end do
                    !$omp end target teams distribute parallel do

                case(3)
                    j = 1
                    !$omp target teams distribute parallel do collapse(2) firstprivate(j) if(is_offload_local)
                    do k = 0, nzp
                        do i = 0, nxp
                            values(i, j-1, k) = values(i, j, k)
                        end do
                    end do
                    !$omp end target teams distribute parallel do

                case(4)
                    j = nyp
                    !$omp target teams distribute parallel do collapse(2) firstprivate(j) if(is_offload_local)
                    do k = 0, nzp
                        do i = 0, nxp
                            values(i, j, k) = values(i, j-1, k)
                        end do
                    end do
                    !$omp end target teams distribute parallel do

                case(5)
                    k = 1
                    !$omp target teams distribute parallel do collapse(2) firstprivate(k) if(is_offload_local)
                    do j = 0, nyp
                        do i = 0, nxp
                            values(i, j, k-1) = values(i, j, k)
                        end do
                    end do
                    !$omp end target teams distribute parallel do

                case(6)
                    k = nzp
                    !$omp target teams distribute parallel do collapse(2) firstprivate(k) if(is_offload_local)
                    do j = 0, nyp
                        do i = 0, nxp
                            values(i, j, k) = values(i, j, k-1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
            end select
        end if
        return
    end subroutine

end module slip_cell_3d_module
