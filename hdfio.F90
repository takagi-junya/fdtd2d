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
    integer :: istat1(20),istat2(20),error(20)
    integer :: rank=2
    integer :: rank2 =1
    integer(HSIZE_T) :: dims(2)
    integer(HSIZE_T) :: dim2(1)

contains

    !HDFの初期化
    subroutine hdfinit()
        use HDF5
        implicit none
        integer(kind=4) :: hdferr
        call h5open_f(hdferr)
        call h5eset_auto_f(0,hdferr)
    end subroutine hdfinit

    !HDFの終了
    subroutine hdffinalize()
        use HDF5
        implicit none
        integer(kind=4) hdferr
        call h5close_f(hdferr)
    end subroutine hdffinalize

    !HDFファイルを開く
    subroutine hdfopen(filename,groupname,file_id,group_id,acc)
        use HDF5
        implicit none
        character(len=*),intent(in) ::  filename,groupname
        integer(kind=HID_T),intent(out) :: file_id,group_id
        integer(kind=4),intent(in) :: acc 
        integer(kind=4) :: hdferr
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

