subroutine finalize
    use constants
    implicit none
    
    deallocate(ex,ey,ez)
    deallocate(hx,hy,hz)
    deallocate(jx,jy,jz)
    deallocate(aex,aey,aez)
    deallocate(bexy,beyx,bezx,bezy)
    deallocate(amx,amy,amz)
    deallocate(bmxy,bmyx,bmzx,bmzy)
    deallocate(epsd,sgmed,mud,sgmmd)
    call mpi_type_free(edge,mpierr)
    call mpi_type_free(pml_edge,mpierr)
    if(pls.ge.1) then
        deallocate(nd)
        deallocate(ajx,ajy,ajz)
        deallocate(ajex,ajey,ajez)
        if(pls.ge.5) then
            deallocate(vx,vy,vz)
            deallocate(avx,avy,avz)
        endif
    endif
end subroutine