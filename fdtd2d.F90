program fdtd2d
    use constants
    use HDF5
    use MPI
    use output
    use hdfio
    use pml2d
    use pwave
    use tfsf_gausian
    use tfsf_dgausian
    use tfsf_sin
    use tfsf_RCP
    use tfsf_LCP
    use omp_lib
    use farfield
    use directivity
    
    implicit none
    real(kind=8) :: time0,time1
    integer :: nthread
    comm = mpi_comm_world
    info = mpi_info_null
    call mpi_init(mpierr)
    call mpi_comm_rank(comm,myrank,mpierr)
    call mpi_comm_size(comm,nprocs,mpierr)

    !$omp parallel
      nthread = omp_get_num_threads()
    !$omp end parallel
    call hdfinit()
    call setup()
    call out_init()
    
    if(wshape.eq.1) then
        call init_ts_gausian()
    else if(wshape.eq.2) then
        call init_ts_dgausian()
    else if(wshape.eq.3) then 
        call init_ts_sin()
    else if(wshape.eq.4) then 
        call init_ts_RCP()
    else if(wshape.eq.5) then 
        call init_ts_LCP()
    endif

    if(mode.eq.1) then
        call init_far()
    else if(mode.eq.2) then 
        call init_dir()
    endif

    call initpml()
    t=dt
    
    call mpi_barrier(comm,mpierr)
    time0 = mpi_wtime()

    do step=1,nstep+1

        if(myrank.eq.0) then
            write(*,'(a10,I5.5)')"Time step:",step-1
        endif

        call efield()

        if(wshape.eq.1) then 
            call e_add_gausian()
        else if(wshape.eq.2) then
            call e_add_dgausian()
        else if(wshape.eq.3) then
            call e_add_sinwave()
        else if(wshape.eq.4) then 
            call e_add_RCPwave()
        else if(wshape.eq.5) then
            call e_add_LCPwave()
        endif

        call epml()
        
        if(mode.eq.1) then 
            call mfarfld()
        else if(mode.eq.2) then
            call k_m()
        endif

        t=t+0.5d0*dt

        if(pls.ge.5) then 
            call velocity()
        endif 
        if((pls.le.4).and.(pls.ge.1)) then 
            call current()
        else if(pls.ge.5) then 
            call currentEOM()
        endif

        call hfield()

        if(wshape.eq.1) then 
            call h_add_gausian()
        else if(wshape.eq.2) then
            call h_add_dgausian()
        else if(wshape.eq.3) then
            call h_add_sinwave()
        else if(wshape.eq.4) then 
            call h_add_RCPwave()
        else if(wshape.eq.5) then
            call h_add_LCPwave()
        endif

        call hpml()

        if(mode.eq.1) then 
            call jfarfld()
        else if(mode.eq.2) then 
            call k_e()
        endif

        t=t+0.5d0*dt
        call out_emf(step-1)
    enddo

    if(mode.eq.1) then 
        call out_far()
        call finalize_far()
    else if(mode.eq.2) then 
        call out_dir()
        call finalize_dir()
    endif

    call mpi_barrier(comm,mpierr)
    time1 = mpi_wtime()
    if(myrank.eq.0) then
        open(50,file="memo.txt",position="append")
        write(50,*)"np:",nprocs,"nt:",nthread," time:",time1-time0
        close(50)
    endif
    call hdffinalize()
    call out_finailzie()
    call finalize()
    call mpi_finalize(mpierr)
end program fdtd2d