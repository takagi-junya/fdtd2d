subroutine finalize
    use constants
    implicit none
    
    deallocate(ex,ey,ez)
    deallocate(hx,hy,hz)
    deallocate(jx,jy,jz)
    deallocate(nd)
    deallocate(aex,aey,aez)
    deallocate(bexy,beyx,bezx,bezy)
    deallocate(amx,amy,amz)
    deallocate(bmxy,bmyx,bmzx,bmzy)
    deallocate(ajx,ajy,ajz)
    deallocate(ajex,ajey,ajez)
    deallocate(epsd,sgmed,mud,sgmmd)
    
    if(pls.ge.5) then
        deallocate(vx,vy,vz)
        deallocate(avx,avy,avz)
    endif
end subroutine