module hdfio
    use HDF5
    use constants
    implicit none
    character(len=20) :: filename(20)
    character(len=10):: groupname(20)
    character(len=10) :: datasetname(10)
    character(len=5) :: tag
    integer(HID_T) :: file_id(20)
    integer(HID_T) :: group_id(20)
    integer(HID_T) :: dataset_id(20)
    integer(HID_T) :: dataspace_id(20)
    integer(HID_T) :: plist_id(20)
    integer :: istat1(20),istat2(20),error(20)
    integer :: rank=2
    integer :: rank2 =1
    integer(kind=4) :: hdferr
    integer(HSIZE_T) :: dims(2)
    integer(HSIZE_T) :: dim2(1)
    integer(HSIZE_T),dimension(2) :: dimsf2d,dimsfi2d,chunk_dims2d

contains

    !HDFの初期化
    subroutine hdfinit()
        use HDF5
        implicit none
        call h5open_f(hdferr)
        call h5eset_auto_f(0,hdferr)
    end subroutine hdfinit

    !HDFの終了
    subroutine hdffinalize()
        use HDF5
        implicit none
        call h5close_f(hdferr)
    end subroutine hdffinalize

    !HDFファイルを開く
    subroutine hdfopen(filename,groupname,file_id,group_id,acc)
        use HDF5
        implicit none
        character(len=*),intent(in) ::  filename,groupname
        integer(kind=HID_T),intent(out) :: file_id,group_id
        integer(kind=4),intent(in) :: acc 
        if(acc==0) then
            call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdferr)
            call h5gcreate_f(file_id,groupname,group_id,hdferr)
        else if(acc==1) then
            call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdferr)
            call h5gopen_f(file_id,groupname,group_id,hdferr)
        else 
            call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdferr)
            call h5gopen_f(file_id,groupname,group_id,hdferr)
        endif
    end subroutine hdfopen

    subroutine h5popen(filename,groupname,file_id,group_id,plist_id,acc,comm,info,err)
        use HDF5
        implicit none 
        integer :: hdferr
        integer,intent(in) :: acc,comm,info
        integer,intent(inout) :: err
        integer(HID_T),intent(inout) :: file_id,group_id,plist_id 
        character(len=5),intent(in) :: filename
        character(len=2),intent(in) :: groupname 
        call h5pcreate_f(h5p_file_access_f,plist_id,err)
        call h5pset_fapl_mpio_f(plist_id,comm,info,err)
        if(acc.eq.0) then
            call h5fcreate_f(filename,h5f_acc_trunc_f,file_id,err,access_prp=plist_id)
            call h5gcreate_f(file_id,groupname,group_id,err)
        else if(acc.eq.1) then
            call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,err,access_prp=plist_id)
            call h5gopen_f(file_id,groupname,group_id,err)
        endif
        call h5pclose_f(plist_id,err)
    end subroutine

    subroutine h5pclose(file_id,group_id,plist_id,istate)
        use HDF5
        implicit none
        integer :: err
        integer(HID_T),intent(inout) :: file_id,group_id,plist_id
        integer(kind=4),intent(in) :: istate 
        call h5gclose_f(group_id,err)
        call h5fclose_f(file_id,err)
    end subroutine

    !HDFファイルを閉じる
    subroutine hdfclose(file_id,group_id,istate)
        use HDF5
        implicit none
        integer(HID_T),intent(in) :: file_id,group_id
        integer(kind=4),intent(out) :: istate
            call h5gclose_f(group_id,istate)
            call h5fclose_f(file_id,istate)
        return
    end subroutine hdfclose

    !1次元配列の書き込み
    subroutine wrt1d(file_id,group_id,datasetname,dim,data,istat,status)
        use HDF5
        implicit none
        integer(kind=4),parameter :: rank=1

        integer(kind=HID_T),intent(in) :: file_id,group_id
        integer(kind=HSIZE_T),intent(in) :: dim(rank)
        integer(kind=4),intent(out) :: istat,status
        real(kind=8),intent(in) :: data(dim(1))
        character(len=*),intent(in) :: datasetname
        integer(kind=HID_T) dataspace_id,dataset_id

        istat=-1
        call h5screate_simple_f(rank,dim,dataspace_id,status)
        call h5dcreate_f(group_id,datasetname,H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,status)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,data,dim,istat)
        call h5dclose_f(dataset_id,status)
        call h5sclose_f(dataspace_id,status)

        return
    end subroutine wrt1d

    !2次元配列の書き込み
    subroutine wrt2d(file_id,group_id,datasetname,dim,data,istat,status)
        use HDF5
        implicit none
        integer(kind=4),parameter :: rank=2

        integer(kind=HID_T),intent(in) :: file_id,group_id
        integer(kind=HSIZE_T),intent(in) :: dim(rank)
        integer(kind=4),intent(out) :: istat,status
        real(kind=8),intent(in) :: data(dim(1),dim(2))
        character(len=*),intent(in) :: datasetname
        integer(kind=HID_T) dataspace_id,dataset_id

        istat=-1
        call h5screate_simple_f(rank,dim,dataspace_id,status)
        call h5dcreate_f(group_id,datasetname,H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,status)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,data,dim,istat)
        call h5dclose_f(dataset_id,status)
        call h5sclose_f(dataspace_id,status)

        return
    end subroutine wrt2d

    subroutine wrt2p(group_id,plist_id,datasetname,dims,chunk_dims,coords,data)
        use HDF5
        use mpi
        implicit none
        integer(kind=4),parameter :: rank=2
        integer,intent(in) :: coords(:)
        integer(kind=HID_T),intent(inout) :: group_id,plist_id
        integer(kind=HSIZE_T),intent(in) :: dims(2),chunk_dims(2)
        real(kind=8),intent(in) :: data(:,:)
        integer(kind=HID_T) :: dataspace,dataset_id
        character(len=*),intent(in) :: datasetname

        integer :: err
        integer(kind=4):: status
        integer(kind=HID_T) :: filespace,dataspace_id,memspace
        integer(HSIZE_T),dimension(2) :: stride,count,block
        integer(HSSIZE_T),dimension(2) :: offset

        call h5screate_simple_f(rank,dims,dataspace_id,status)
        call h5screate_simple_f(rank,chunk_dims,memspace,status)

        call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,status)
        call h5pset_chunk_f(plist_id,rank,chunk_dims,status)
        call h5dcreate_f(group_id,datasetname,H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,status)
        call h5sclose_f(dataspace_id,status)
        
        stride(1) = 1
        stride(2) = 1 
        count(1) = 1
        count(2) = 1 
        block(1) = chunk_dims(1)
        block(2) = chunk_dims(2)
        offset(1) = chunk_dims(1)*coords(1)
        offset(2) = chunk_dims(2)*coords(2)
        
        call h5dget_space_f(dataset_id,dataspace_id,status)
        call h5sselect_hyperslab_f(dataspace_id,H5S_SELECT_SET_F,offset,count,status,stride,block)

        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,status)
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,status)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,data,dims,status,file_space_id=dataspace_id,mem_space_id=memspace,xfer_prp=plist_id)

        call h5sclose_f(dataspace_id,status)
        call h5sclose_f(memspace,status)
        call h5pclose_f(plist_id,status)
        call h5dclose_f(dataset_id,status)
        call h5sclose_f(dataspace_id,status)
    end subroutine

    !3次元配列の書き込み
    subroutine wrt3d(file_id,group_id,datasetname,dim,data,istat,status)
        use HDF5
        implicit none
        integer(kind=4),parameter :: rank=3

        integer(kind=HID_T),intent(in) :: file_id,group_id
        integer(kind=HSIZE_T),intent(in) :: dim(rank)
        integer(kind=4),intent(out) :: istat,status
        real(kind=8),intent(in) :: data(dim(1),dim(2),dim(3))
        character(len=*),intent(in) :: datasetname
        integer(kind=HID_T) dataspace_id,dataset_id

        istat=-1
        call h5screate_simple_f(rank,dim,dataspace_id,status)
        call h5dcreate_f(group_id,datasetname,H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,status)
        call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,data,dim,istat)
        call h5dclose_f(dataset_id,status)
        call h5sclose_f(dataspace_id,status)

        return
    end subroutine wrt3d

end module hdfio

