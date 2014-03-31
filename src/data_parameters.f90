module data_parameters
      integer,parameter :: &
      ireals=selected_real_kind(15,200), &
      iintegers=selected_int_kind(12)
      integer :: mpiint_dummy
      integer,parameter :: mpiint = kind(mpiint_dummy)
      real(ireals),parameter :: pi=3.141592653589793, nil=-9999._ireals
      real(ireals),parameter :: zero=0._ireals, one=1._ireals
      integer(iintegers) ,parameter :: i0=0,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9,i10=10,i11=11,inil=-9999_iintegers

end module
