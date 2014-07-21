
subroutine update_dir_stream_1_2(bmc,p,S,N)
    class(t_boxmc_1_2) :: bmc
    type(photon),intent(in) :: p
    real(ireals),intent(inout) :: S(:)
    integer(iintegers),intent(inout) :: N(:)

    select case (p%side)
    case(1:2)
      S(1) = S(1)+p%weight
      N(1) = N(1)+i1
      return
    case default
      print *,'Dont know what to do with this p%side'
      call print_photon(p)
    end select
    select type(bmc)
    end select

end subroutine
subroutine update_diff_stream_1_2(bmc,p,S,N)
    class(t_boxmc_1_2) :: bmc
    type(photon),intent(in) :: p
    real(ireals),intent(inout) :: S(:)
    integer(iintegers),intent(inout) :: N(:)

    !         _______1_______
    !        |           10  |
    !      5 |            8  | 6
    !        |               |
    !        |   9           |
    !        |   7           |
    !      3 |_______________| 4
    !                2

    if(p%side.eq.1) then
      S(1) = S(1)+p%weight
      N(1) = N(1)+i1
      return
    else if(p%side.eq.2) then
      S(2) = S(2)+p%weight
      N(2) = N(2)+i1
      return
    else
      print *,'Couldnt find a stream on to which I can put the photon weight on?!'
      call print_photon(p)
      call exit
    endif
    select type(bmc)
    end select

end subroutine

subroutine init_dir_photon_1_2(bmc,p,src,direct,initial_dir,dx,dy,dz)
    class(t_boxmc_1_2) :: bmc
    type(photon),intent(inout) :: p
    real(ireals),intent(in) :: dx,dy,dz,initial_dir(3)
    integer(iintegers),intent(in) :: src
    logical,intent(in) :: direct

    p%alive = .False.

    if(src.eq.1) then
      p%loc = (/L(dx)     , L(dy)     ,    dz  /)
    else
      print *,'Dont know what to do with source spec:',src
      call exit
    endif

    p%weight=one
    p%dx   = dx
    p%dy   = dy
    p%dz   = dz
    p%alive = .True.
    p%direct= direct
    p%side = int(nil)
    p%src  = src
    p%dir = initial_dir

    select type(bmc)
    end select

end subroutine

subroutine init_diff_photon_1_2(bmc,p,src,dx,dy,dz)
    class(t_boxmc_1_2) :: bmc
    type(photon),intent(inout) :: p
    real(ireals),intent(in) :: dx,dy,dz
    integer(iintegers),intent(in) :: src

    real(ireals) :: phi,mu

    if(src.eq.1) then
      p%loc = (/L(dx), L(dy),    dz  /)
      mu = interv_R( deg2mu(90._ireals),deg2mu(180._ireals) ) ; phi= deg2rad( R()*360._ireals )
    else if(src.eq.2) then
      p%loc = (/L(dx), L(dy),    zero/)
      mu = interv_R( deg2mu(0._ireals),deg2mu(90._ireals) ) ; phi= deg2rad( R()*360._ireals )
    else
      print *,'Dont know what to do with source spec:',src
      call exit
    endif

    p%dir = (/cos(phi)*sin(acos(mu)) , sin(phi)*sin(acos(mu)) , mu  /)
    p%weight=one
    p%dx   = dx
    p%dy   = dy
    p%dz   = dz
    p%alive = .True.
    p%direct= .False.
    p%side = int(nil)
    p%src  = src
    if(.not. approx( norm(p%dir) ,one ) ) stop 'init_diff_photon_1_2 :: direction was not normed?!'

    select type(bmc)
    end select

end subroutine

subroutine intersect_distance_1_2(bmc,p,max_dist)
    class(t_boxmc_1_2) :: bmc
    type(photon),intent(inout) :: p
    real(ireals),intent(out) :: max_dist

    real(ireals) :: x,y
    !        real(ireals),parameter :: eps=1e-3

    !crossing with bottom and top plane:
    if(p%dir(3).ge.zero) then
      max_dist = hit_plane(p,[zero,zero,p%dz ],[zero,zero,one])
      p%side=1 
      x = p%loc(1)+p%dir(1)*max_dist
      y = p%loc(2)+p%dir(2)*max_dist
      return
    endif
    if(p%dir(3).le.zero) then
      max_dist = hit_plane(p,[zero,zero,zero ],[zero,zero,one])
      p%side=2
      x = p%loc(1)+p%dir(1)*max_dist
      y = p%loc(2)+p%dir(2)*max_dist
      return
    endif

    print *,'should actually not be here at the end of crossings in intersect distance! - however, please check if distance makes sense?:',max_dist
    call print_photon(p)

    select type(bmc)
    end select

end subroutine