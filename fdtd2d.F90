program fdtd2d
    use constants
    use HDF5
    use hdfio
    use pml2d
    use pwave
    use tfsf_gausian
    use tfsf_dgausian
    use tfsf_sin
    use omp_lib
    use farfield
    use directivity
    
    implicit none
    open(30,file="sim.out")
    call hdfinit()
    call setup()
    
    if(kwave.eq.0.and.mode.eq.0) then
        call init_pwave()
        call initpml()
    endif
    if(mode.eq.1.or.mode.eq.4) then
        call initpml()
        call init_ts_gausian()
        call init_far()
    endif
    if(mode.eq.2) then
        call initpml()
        call init_ts_sin()
        call init_dir()
    endif
    if(mode.eq.3.or.mode.eq.6.or.mode.eq.7) then
        call initpml()
        call init_ts_sin()
        call init_dir()
    endif
    if(mode.eq.5) then 
        call initpml()
        call init_ts_dgausian()
        call init_far()
    endif
    if(mode.eq.8) then
        call initpml()
        call init_ts_dgausian()
        call init_far()
    endif
    t=dt
    
    if(mode.eq.0) then
        write(30,'(a21)'),"All total field mode"
        do step=1,nstep+1
            write(*,'(a10,I5.5)')"Time step:",step-1
            call h_pbc()
            call efield()
            t=t+0.5d0*dt
            call e_pbc()
            call hfield()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo
    elseif(mode.eq.1) then
        write(30,'(a17)'),"TF-SF pulse mode"
        do step=1,nstep+1
            write(*,'(a10,I5.5)')"Time step:",step-1
            call efield()
            call e_add_gausian()
            call epml()
            call mfarfld()
            t=t+0.5d0*dt
            call hfield()
            call h_add_gausian()
            call hpml()
            call jfarfld()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo
        call out_far()
    elseif(mode.eq.2) then
        write(30,'(a15)'),"TF-SF sin mode"
        do step=1,nstep+1
            write(*,'(a10,I5.5)')"Time step:",step-1
            call efield()
            call e_add_sin()
            call epml()
            call k_m()
            t=t+0.5d0*dt
            call hfield()
            call h_add_sin()
            call hpml()
            call k_e()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo
        call out_dir()
    elseif(mode.eq.3) then
        write(30,'(a16)'),"plasma sin mode"
        do step=1,nstep+1
            write(*,'(a10,I5.5)')"Time step:",step-1
            call efield()
            call e_add_sin()
            call epml()
            call k_m()
            t=t+0.5d0*dt
            call current()
            call hfield()
            call h_add_sin()
            call hpml()
            call k_e()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo
        call out_dir()
    elseif(mode.eq.4) then
        write(30,'(a18)'),"plasma gauss mode"
        do step=1,nstep+1
            write(*,'(a10,I5.5)')"Time step:",step-1
            call efield()
            call e_add_gausian()
            call epml()
            call mfarfld()
            t=t+0.5d0*dt
            call hfield()
            call current()
            call h_add_gausian()
            call hpml()
            call jfarfld()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo    
        call out_far()
    elseif(mode.eq.5) then
        write(30,'(a19)'),"plasma gauss mode2"
        do step=1,nstep+1
            write(*,'(a10,I5.5)')"Time step:",step-1
            call efield()
            call e_add_dgausian()
            call epml()
            call mfarfld()
            t=t+0.5d0*dt
            call hfield()
            call current()
            call h_add_dgausian()
            call hpml()
            call jfarfld()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo    
        call out_far()
    elseif(mode.eq.6) then
        write(30,'(a16)'),"plasma sin mode"
        do step=1,nstep+1
            write(*,'(a10,I5.5)')"Time step:",step-1
            call efield()
            call e_add_sin()
            call epml()
            call k_m()
            t=t+0.5d0*dt
            call current()
            call hfield()
            call h_add_sin()
            call hpml()
            call k_e()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo
        call out_dir()
    elseif(mode.eq.7) then
        write(30,'(a23)'),"plasma moving equation"
        do step=1,nstep+1
            write(*,'(a10,I5.5)')"Time step:",step-1
            call efield()
            call e_add_sin()
            call epml()
            call k_m()
            t=t+0.5d0*dt
            call velocity()
            call currentEOM()
            call hfield()
            call h_add_sin()
            call hpml()
            call k_e()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo
        call out_dir()
    elseif(mode.eq.8) then
        write(30,'(a23)'),"plasma sin pulse"
        do step=1,nstep+1
            write(*,'(a10,I5.5)')"Time step:",step-1
            call efield()
            call e_add_dgausian()
            call epml()
            call mfarfld()
            t=t+0.5d0*dt
            call current()
            call hfield()
            call h_add_dgausian()
            call hpml()
            call jfarfld()
            t=t+0.5d0*dt
            call out_emf(step-1)
        enddo
        call out_far()
    endif
    call hdffinalize()
    call finalize()
    close(30)
end program fdtd2d