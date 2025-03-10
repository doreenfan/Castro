! This subroutine cannot be tiled
subroutine ca_compute_lamborder(Er, Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3, &
                                kap, kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
                                lam, lam_l1, lam_l2, lam_l3, lam_h1, lam_h2, lam_h3, &
                                dx, ngrow, limiter, filter_T, S)

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDlambda
  use filter_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3, &
       kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
       lam_l1, lam_l2, lam_l3, lam_h1, lam_h2, lam_h3
  integer, intent(in) :: ngrow, limiter, filter_T, S
  real(rt)        , intent(in) :: dx(3)
  real(rt)        , intent(in) :: kap(kap_l1:kap_h1, kap_l2:kap_h2, kap_l3:kap_h3, 0:ngroups-1)
  real(rt)        , intent(in) :: Er(Er_l1:Er_h1, Er_l2:Er_h2, Er_l3:Er_h3, 0:ngroups-1)
  real(rt)        , intent(out) :: lam(lam_l1:lam_h1, lam_l2:lam_h2, lam_l3:lam_h3, 0:ngroups-1)

  integer :: i, j, k, reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3, g
  real(rt)         :: r, r1, r2, r3

  real(rt)        , allocatable :: lamfil(:,:,:)

  reg_l1 = lam_l1 + ngrow
  reg_l2 = lam_l2 + ngrow
  reg_l3 = lam_l3 + ngrow
  reg_h1 = lam_h1 - ngrow
  reg_h2 = lam_h2 - ngrow
  reg_h3 = lam_h3 - ngrow

  if (filter_T .gt. 0) then
     allocate(lamfil(reg_l1:reg_h1,lam_l2:lam_h2,lam_l3:lam_h3))
  end if

  do g = 0, ngroups-1

  do k=lam_l3, lam_h3
     do j=lam_l2, lam_h2
        do i=lam_l1, lam_h1

           lam(i,j,k,g) = -1.e50_rt

           if (Er(i,j,k,g) .eq. -1.e0_rt) then
              cycle
           end if

           if (Er(i-1,j,k,g) .eq. -1.e0_rt) then
              r1 = (Er(i+1,j,k,g) - Er(i,j,k,g)) / (dx(1))
           else if (Er(i+1,j,k,g) .eq. -1.e0_rt) then
              r1 = (Er(i,j,k,g) - Er(i-1,j,k,g)) / (dx(1))
           else
              r1 = (Er(i+1,j,k,g) - Er(i-1,j,k,g)) / (2.e0_rt*dx(1))
           end if

           if (Er(i,j-1,k,g) .eq. -1.e0_rt) then
              r2 = (Er(i,j+1,k,g) - Er(i,j,k,g)) / (dx(2))
           else if (Er(i,j+1,k,g) .eq. -1.e0_rt) then
              r2 = (Er(i,j,k,g) - Er(i,j-1,k,g)) / (dx(2))
           else
              r2 = (Er(i,j+1,k,g) - Er(i,j-1,k,g)) / (2.e0_rt*dx(2))
           end if

           if (Er(i,j,k-1,g) .eq. -1.e0_rt) then
              r3 = (Er(i,j,k+1,g) - Er(i,j,k,g)) / dx(3)
           else if (Er(i,j,k+1,g) .eq. -1.e0_rt) then
              r3 = (Er(i,j,k,g) - Er(i,j,k-1,g)) / dx(3)
           else
              r3 = (Er(i,j,k+1,g) - Er(i,j,k-1,g)) / (2.e0_rt*dx(3))
           end if

           r = sqrt(r1**2 + r2**2 + r3**2)
           r = r / (kap(i,j,k,g) * max(Er(i,j,k,g), 1.e-50_rt))

           lam(i,j,k,g) = FLDlambda(r, limiter)
        end do
     end do
  end do

  ! filter
  if (filter_T .eq. 1) then
     do k=lam_l3, lam_h3
        do j=lam_l2, lam_h2
           if (Er(reg_l1,j,k,g) .eq. -1.e0_rt) then
              lamfil(:,j,k) = -1.e-50_rt
              cycle
           end if

           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff1(0) * lam(i,j,k,g) &
                   &        + ff1(1) * (lam(i-1,j,k,g)+lam(i+1,j,k,g))
           end do

           if (Er(reg_l1-1,j,k,g) .eq. -1.e0_rt) then
              i = reg_l1
              lamfil(i,j,k) = dot_product(ff1b, lam(i:i+1,j,k,g))
           end if

           if (Er(reg_h1+1,j,k,g) .eq. -1.e0_rt) then
              i = reg_h1
              lamfil(i,j,k) = dot_product(ff1b(1:0:-1), lam(i-1:i,j,k,g))
           end if
        end do
     end do

     do k=lam_l3, lam_h3
        if (Er(reg_l1,reg_l2,k,g) .eq. -1.e0_rt) then
           cycle
        end if

        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lam(i,j,k,g) = ff1(0) * lamfil(i,j,k) &
                   &       + ff1(1) * (lamfil(i,j-1,k)+lamfil(i,j+1,k))
           end do
        end do

        if (Er(reg_l1,reg_l2-1,k,g) .eq. -1.e0_rt) then
           j = reg_l2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff1b, lamfil(i,j:j+1,k))
           end do
        end if

        if (Er(reg_l1,reg_h2+1,k,g) .eq. -1.e0_rt) then
           j = reg_h2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff1b(1:0:-1), lamfil(i,j-1:j,k))
           end do
        end if
     end do

     do k=reg_l3, reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff1(0) * lam(i,j,k,g) &
                   &        + ff1(1) * (lam(i,j,k-1,g)+lam(i,j,k+1,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end do

     if (Er(reg_l1,reg_l2,reg_l3-1,g) .eq. -1.e0_rt) then
        k = reg_l3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff1b, lam(i,j,k:k+1,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end if

     if (Er(reg_l1,reg_l2,reg_h3+1,g) .eq. -1.e0_rt) then
        k = reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff1b(1:0:-1), lam(i,j,k-1:k,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end if

     lam(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3,g) = &
          lamfil(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)

  else if (filter_T .eq. 2) then

     do k=lam_l3, lam_h3
        do j=lam_l2, lam_h2
           if (Er(reg_l1,j,k,g) .eq. -1.e0_rt) then
              lamfil(:,j,k) = -1.e-50_rt
              cycle
           end if

           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff2(0,S) * lam(i,j,k,g) &
                   &        + ff2(1,S) * (lam(i-1,j,k,g)+lam(i+1,j,k,g)) &
                   &        + ff2(2,S) * (lam(i-2,j,k,g)+lam(i+2,j,k,g))
           end do

           if (Er(reg_l1-1,j,k,g) .eq. -1.e0_rt) then
              i = reg_l1
              lamfil(i,j,k) = dot_product(ff2b0, lam(i:i+2,j,k,g))

              i = reg_l1 + 1
              lamfil(i,j,k) = dot_product(ff2b1, lam(i-1:i+2,j,k,g))
           end if

           if (Er(reg_h1+1,j,k,g) .eq. -1.e0_rt) then
              i = reg_h1 - 1
              lamfil(i,j,k) = dot_product(ff2b1(2:-1:-1), lam(i-2:i+1,j,k,g))

              i = reg_h1
              lamfil(i,j,k) = dot_product(ff2b0(2:0:-1), lam(i-2:i,j,k,g))
           end if
        end do
     end do

     do k=lam_l3, lam_h3
        if (Er(reg_l1,reg_l2,k,g) .eq. -1.e0_rt) then
           cycle
        end if

        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lam(i,j,k,g) = ff2(0,S) * lamfil(i,j,k) &
                   &       + ff2(1,S) * (lamfil(i,j-1,k)+lamfil(i,j+1,k)) &
                   &       + ff2(2,S) * (lamfil(i,j-2,k)+lamfil(i,j+2,k))
           end do
        end do

        if (Er(reg_l1,reg_l2-1,k,g) .eq. -1.e0_rt) then
           j = reg_l2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff2b0, lamfil(i,j:j+2,k))
           end do

           j = reg_l2 + 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff2b1, lamfil(i,j-1:j+2,k))
           end do
        end if

        if (Er(reg_l1,reg_h2+1,k,g) .eq. -1.e0_rt) then
           j = reg_h2 - 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff2b1(2:-1:-1), lamfil(i,j-2:j+1,k))
           end do

           j = reg_h2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff2b0(2:0:-1), lamfil(i,j-2:j,k))
           end do
        end if
     end do

     do k=reg_l3, reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff2(0,S) * lam(i,j,k,g) &
                   &        + ff2(1,S) * (lam(i,j,k-1,g)+lam(i,j,k+1,g)) &
                   &        + ff2(2,S) * (lam(i,j,k-2,g)+lam(i,j,k+2,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end do

     if (Er(reg_l1,reg_l2,reg_l3-1,g) .eq. -1.e0_rt) then
        k = reg_l3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff2b0, lam(i,j,k:k+2,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff2b1, lam(i,j,k-1:k+2,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end if

     if (Er(reg_l1,reg_l2,reg_h3+1,g) .eq. -1.e0_rt) then
        k = reg_h3 - 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff2b1(2:-1:-1), lam(i,j,k-2:k+1,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff2b0(2:0:-1), lam(i,j,k-2:k,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end if

     lam(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3,g) = &
          lamfil(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)

  else if (filter_T .eq. 3) then

     do k=lam_l3, lam_h3
        do j=lam_l2, lam_h2
           if (Er(reg_l1,j,k,g) .eq. -1.e0_rt) then
              lamfil(:,j,k) = -1.e-50_rt
              cycle
           end if

           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff3(0,S) * lam(i,j,k,g) &
                   &        + ff3(1,S) * (lam(i-1,j,k,g)+lam(i+1,j,k,g)) &
                   &        + ff3(2,S) * (lam(i-2,j,k,g)+lam(i+2,j,k,g)) &
                   &        + ff3(3,S) * (lam(i-3,j,k,g)+lam(i+3,j,k,g))
           end do

           if (Er(reg_l1-1,j,k,g) .eq. -1.e0_rt) then
              i = reg_l1
              lamfil(i,j,k) = dot_product(ff3b0, lam(i:i+3,j,k,g))

              i = reg_l1 + 1
              lamfil(i,j,k) = dot_product(ff3b1, lam(i-1:i+3,j,k,g))

              i = reg_l1 + 2
              lamfil(i,j,k) = dot_product(ff3b2, lam(i-2:i+3,j,k,g))
           end if

           if (Er(reg_h1+1,j,k,g) .eq. -1.e0_rt) then
              i = reg_h1 - 2
              lamfil(i,j,k) = dot_product(ff3b2(3:-2:-1), lam(i-3:i+2,j,k,g))

              i = reg_h1 - 1
              lamfil(i,j,k) = dot_product(ff3b1(3:-1:-1), lam(i-3:i+1,j,k,g))

              i = reg_h1
              lamfil(i,j,k) = dot_product(ff3b0(3:0:-1), lam(i-3:i,j,k,g))
           end if
        end do
     end do

     do k=lam_l3, lam_h3
        if (Er(reg_l1,reg_l2,k,g) .eq. -1.e0_rt) then
           cycle
        end if

        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lam(i,j,k,g) = ff3(0,S) * lamfil(i,j,k) &
                   &       + ff3(1,S) * (lamfil(i,j-1,k)+lamfil(i,j+1,k)) &
                   &       + ff3(2,S) * (lamfil(i,j-2,k)+lamfil(i,j+2,k)) &
                   &       + ff3(3,S) * (lamfil(i,j-3,k)+lamfil(i,j+3,k))
           end do
        end do

        if (Er(reg_l1,reg_l2-1,k,g) .eq. -1.e0_rt) then
           j = reg_l2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b0, lamfil(i,j:j+3,k))
           end do

           j = reg_l2 + 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b1, lamfil(i,j-1:j+3,k))
           end do

           j = reg_l2 + 2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b2, lamfil(i,j-2:j+3,k))
           end do
        end if

        if (Er(reg_l1,reg_h2+1,k,g) .eq. -1.e0_rt) then
           j = reg_h2 - 2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b2(3:-2:-1), lamfil(i,j-3:j+2,k))
           end do

           j = reg_h2 - 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b1(3:-1:-1), lamfil(i,j-3:j+1,k))
           end do

           j = reg_h2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b0(3:0:-1), lamfil(i,j-3:j,k))
           end do
        end if
     end do

     do k=reg_l3, reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff3(0,S) * lam(i,j,k,g) &
                   &        + ff3(1,S) * (lam(i,j,k-1,g)+lam(i,j,k+1,g)) &
                   &        + ff3(2,S) * (lam(i,j,k-2,g)+lam(i,j,k+2,g)) &
                   &        + ff3(3,S) * (lam(i,j,k-3,g)+lam(i,j,k+3,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end do

     if (Er(reg_l1,reg_l2,reg_l3-1,g) .eq. -1.e0_rt) then
        k = reg_l3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b0, lam(i,j,k:k+3,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b1, lam(i,j,k-1:k+3,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 2
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b2, lam(i,j,k-2:k+3,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end if

     if (Er(reg_l1,reg_l2,reg_h3+1,g) .eq. -1.e0_rt) then
        k = reg_h3 - 2
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b2(3:-2:-1), lam(i,j,k-3:k+2,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3 - 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b1(3:-1:-1), lam(i,j,k-3:k+1,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b0(3:0:-1), lam(i,j,k-3:k,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end if

     lam(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3,g) = &
          lamfil(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)

  else if (filter_T .eq. 4) then

     do k=lam_l3, lam_h3
        do j=lam_l2, lam_h2
           if (Er(reg_l1,j,k,g) .eq. -1.e0_rt) then
              lamfil(:,j,k) = -1.e-50_rt
              cycle
           end if

           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff4(0,S) * lam(i,j,k,g) &
                   &        + ff4(1,S) * (lam(i-1,j,k,g)+lam(i+1,j,k,g)) &
                   &        + ff4(2,S) * (lam(i-2,j,k,g)+lam(i+2,j,k,g)) &
                   &        + ff4(3,S) * (lam(i-3,j,k,g)+lam(i+3,j,k,g)) &
                   &        + ff4(4,S) * (lam(i-4,j,k,g)+lam(i+4,j,k,g))
           end do

           if (Er(reg_l1-1,j,k,g) .eq. -1.e0_rt) then
              i = reg_l1
              lamfil(i,j,k) = dot_product(ff4b0, lam(i:i+4,j,k,g))

              i = reg_l1 + 1
              lamfil(i,j,k) = dot_product(ff4b1, lam(i-1:i+4,j,k,g))

              i = reg_l1 + 2
              lamfil(i,j,k) = dot_product(ff4b2, lam(i-2:i+4,j,k,g))

              i = reg_l1 + 3
              lamfil(i,j,k) = dot_product(ff4b3, lam(i-3:i+4,j,k,g))
           end if

           if (Er(reg_h1+1,j,k,g) .eq. -1.e0_rt) then
              i = reg_h1 - 3
              lamfil(i,j,k) = dot_product(ff4b3(4:-3:-1), lam(i-4:i+3,j,k,g))

              i = reg_h1 - 2
              lamfil(i,j,k) = dot_product(ff4b2(4:-2:-1), lam(i-4:i+2,j,k,g))

              i = reg_h1 - 1
              lamfil(i,j,k) = dot_product(ff4b1(4:-1:-1), lam(i-4:i+1,j,k,g))

              i = reg_h1
              lamfil(i,j,k) = dot_product(ff4b0(4:0:-1), lam(i-4:i,j,k,g))
           end if
        end do
     end do

     do k=lam_l3, lam_h3
        if (Er(reg_l1,reg_l2,k,g) .eq. -1.e0_rt) then
           cycle
        end if

        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lam(i,j,k,g) = ff4(0,S) * lamfil(i,j,k) &
                   &       + ff4(1,S) * (lamfil(i,j-1,k)+lamfil(i,j+1,k)) &
                   &       + ff4(2,S) * (lamfil(i,j-2,k)+lamfil(i,j+2,k)) &
                   &       + ff4(3,S) * (lamfil(i,j-3,k)+lamfil(i,j+3,k)) &
                   &       + ff4(4,S) * (lamfil(i,j-4,k)+lamfil(i,j+4,k))
           end do
        end do

        if (Er(reg_l1,reg_l2-1,k,g) .eq. -1.e0_rt) then
           j = reg_l2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b0, lamfil(i,j:j+4,k))
           end do

           j = reg_l2 + 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b1, lamfil(i,j-1:j+4,k))
           end do

           j = reg_l2 + 2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b2, lamfil(i,j-2:j+4,k))
           end do

           j = reg_l2 + 3
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b3, lamfil(i,j-3:j+4,k))
           end do
        end if

        if (Er(reg_l1,reg_h2+1,k,g) .eq. -1.e0_rt) then
           j = reg_h2 - 3
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b3(4:-3:-1), lamfil(i,j-4:j+3,k))
           end do

           j = reg_h2 - 2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b2(4:-2:-1), lamfil(i,j-4:j+2,k))
           end do

           j = reg_h2 - 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b1(4:-1:-1), lamfil(i,j-4:j+1,k))
           end do

           j = reg_h2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b0(4:0:-1), lamfil(i,j-4:j,k))
           end do
        end if
     end do

     do k=reg_l3, reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff4(0,S) * lam(i,j,k,g) &
                   &        + ff4(1,S) * (lam(i,j,k-1,g)+lam(i,j,k+1,g)) &
                   &        + ff4(2,S) * (lam(i,j,k-2,g)+lam(i,j,k+2,g)) &
                   &        + ff4(3,S) * (lam(i,j,k-3,g)+lam(i,j,k+3,g)) &
                   &        + ff4(4,S) * (lam(i,j,k-4,g)+lam(i,j,k+4,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end do

     if (Er(reg_l1,reg_l2,reg_l3-1,g) .eq. -1.e0_rt) then
        k = reg_l3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b0, lam(i,j,k:k+4,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b1, lam(i,j,k-1:k+4,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 2
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b2, lam(i,j,k-2:k+4,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b3, lam(i,j,k-3:k+4,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end if

     if (Er(reg_l1,reg_l2,reg_h3+1,g) .eq. -1.e0_rt) then
        k = reg_h3 - 3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b3(4:-3:-1), lam(i,j,k-4:k+3,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3 - 2
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b2(4:-2:-1), lam(i,j,k-4:k+2,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3 - 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b1(4:-1:-1), lam(i,j,k-4:k+1,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b0(4:0:-1), lam(i,j,k-4:k,g))
              lamfil(i,j,k) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i,j,k)))
           end do
        end do
     end if

     lam(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3,g) = &
          lamfil(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)

  end if

  ! lo-x lo-y lo-z
  do k=lam_l3,reg_l3-1
     do j=lam_l2,reg_l2-1
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_l1,reg_l2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! reg-x lo-y lo-z
  do k=lam_l3,reg_l3-1
     do j=lam_l2,reg_l2-1
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(i,reg_l2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! hi-x lo-y lo-z
  do k=lam_l3,reg_l3-1
     do j=lam_l2,reg_l2-1
        do i=reg_h1+1,lam_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_h1,reg_l2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! lo-x reg-y lo-z
  do k=lam_l3,reg_l3-1
     do j=reg_l2,reg_h2
        do i=lam_l1,reg_l1-1
           lam(i,j,k,g) = lam(reg_l1,j,reg_l3,g)
        end do
     end do
  end do

  ! reg-x reg-y lo-z side
  do k=lam_l3,reg_l3-1
     do j=reg_l2,reg_h2
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(i,j,reg_l3,g)
           end if
        end do
     end do
  end do

  ! hi-x reg-y lo-z
  do k=lam_l3,reg_l3-1
     do j=reg_l2,reg_h2
        do i=reg_h1+1,lam_h1
           lam(i,j,k,g) = lam(reg_h1,j,reg_l3,g)
        end do
     end do
  end do

  ! lo-x hi-y lo-z
  do k=lam_l3,reg_l3-1
     do j=reg_h2+1,lam_h2
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_l1,reg_h2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! reg-x hi-y lo-z
  do k=lam_l3,reg_l3-1
     do j=reg_h2+1,lam_h2
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(i,reg_h2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! hi-x hi-y lo-z
  do k=lam_l3,reg_l3-1
     do j=reg_h2+1,lam_h2
        do i=reg_h1+1,lam_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_h1,reg_h2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! lo-x lo-y reg-z
  do k=reg_l3,reg_h3
     do j=lam_l2,reg_l2-1
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_l1,reg_l2,k,g)
           end if
        end do
     end do
  end do

  ! reg-x lo-y reg-z
  do k=reg_l3,reg_h3
     do j=lam_l2,reg_l2-1
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(i,reg_l2,k,g)
           end if
        end do
     end do
  end do

  ! hi-x lo-y reg-z
  do k=reg_l3,reg_h3
     do j=lam_l2,reg_l2-1
        do i=reg_h1+1,lam_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_h1,reg_l2,k,g)
           end if
        end do
     end do
  end do

  ! lo-x reg-y reg-z
  do k=reg_l3,reg_h3
     do j=reg_l2,reg_h2
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_l1,j,k,g)
           end if
        end do
     end do
  end do

  ! hi-x reg-y reg-z
  do k=reg_l3,reg_h3
     do j=reg_l2,reg_h2
        do i=reg_h1+1,lam_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_h1,j,k,g)
           end if
        end do
     end do
  end do

  ! lo-x hi-y reg-z
  do k=reg_l3,reg_h3
     do j=reg_h2+1,lam_h2
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_l1,reg_h2,k,g)
           end if
        end do
     end do
  end do

  ! reg-x hi-y reg-z
  do k=reg_l3,reg_h3
     do j=reg_h2+1,lam_h2
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(i,reg_h2,k,g)
           end if
        end do
     end do
  end do

  ! hi-x hi-y reg-z
  do k=reg_l3,reg_h3
     do j=reg_h2+1,lam_h2
        do i=reg_h1+1,lam_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_h1,reg_h2,k,g)
           end if
        end do
     end do
  end do

  ! lo-x lo-y hi-z
  do k=reg_h3+1,lam_h3
     do j=lam_l2,reg_l2-1
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_l1,reg_l2,reg_h3,g)
           end if
        end do
     end do
  end do

  ! reg-x lo-y hi-z
  do k=reg_h3+1,lam_h3
     do j=lam_l2,reg_l2-1
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(i,reg_l2,reg_h3,g)
           end if
        end do
     end do
  end do

  ! hi-x lo-y hi-z
  do k=reg_h3+1,lam_h3
     do j=lam_l2,reg_l2-1
        do i=reg_h1+1,lam_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_h1,reg_l2,reg_h3,g)
           end if
        end do
     end do
  end do

  ! lo-x reg-y hi-z
  do k=reg_h3+1,lam_h3
     do j=reg_l2,reg_h2
        do i=lam_l1,reg_l1-1
           lam(i,j,k,g) = lam(reg_l1,j,reg_h3,g)
        end do
     end do
  end do

  ! reg-x reg-y hi-z
  do k=reg_h3+1,lam_h3
     do j=reg_l2,reg_h2
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(i,j,reg_h3,g)
           end if
        end do
     end do
  end do

  ! hi-x reg-y hi-z
  do k=reg_h3+1,lam_h3
     do j=reg_l2,reg_h2
        do i=reg_h1+1,lam_h1
           lam(i,j,k,g) = lam(reg_h1,j,reg_h3,g)
        end do
     end do
  end do

  ! lo-x hi-y hi-z
  do k=reg_h3+1,lam_h3
     do j=reg_h2+1,lam_h2
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_l1,reg_h2,reg_h3,g)
           end if
        end do
     end do
  end do

  ! reg-x hi-y hi-z
  do k=reg_h3+1,lam_h3
     do j=reg_h2+1,lam_h2
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(i,reg_h2,reg_h3,g)
           end if
        end do
     end do
  end do

  ! hi-x hi-y hi-z
  do k=reg_h3+1,lam_h3
     do j=reg_h2+1,lam_h2
        do i=reg_h1+1,lam_h1
           if (Er(i,j,k,g).eq.-1.e0_rt) then
              lam(i,j,k,g) = lam(reg_h1,reg_h2,reg_h3,g)
           end if
        end do
     end do
  end do

  end do ! do g=

  if (filter_T .gt. 0) then
     deallocate(lamfil)
  end if

  return
end subroutine ca_compute_lamborder

subroutine ca_correct_dterm(  &
                            dfx, dfx_l1, dfx_l2, dfx_l3, dfx_h1, dfx_h2, dfx_h3, &
                            dfy, dfy_l1, dfy_l2, dfy_l3, dfy_h1, dfy_h2, dfy_h3, &
                            dfz, dfz_l1, dfz_l2, dfz_l3, dfz_h1, dfz_h2, dfz_h3, &
                            re, rc) bind(C, name="ca_correct_dterm")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: dfx_l1, dfx_l2, dfx_l3, dfx_h1, dfx_h2, dfx_h3
  integer, intent(in) :: dfy_l1, dfy_l2, dfy_l3, dfy_h1, dfy_h2, dfy_h3
  integer, intent(in) :: dfz_l1, dfz_l2, dfz_l3, dfz_h1, dfz_h2, dfz_h3
  real(rt)        , intent(inout) :: dfx(dfx_l1:dfx_h1,dfx_l2:dfx_h2,dfx_l3:dfx_h3)
  real(rt)        , intent(inout) :: dfy(dfy_l1:dfy_h1,dfy_l2:dfy_h2,dfy_l3:dfy_h3)
  real(rt)        , intent(inout) :: dfz(dfz_l1:dfz_h1,dfz_l2:dfz_h2,dfz_l3:dfz_h3)
  real(rt)        , intent(in) :: re(1), rc(1)

end subroutine ca_correct_dterm
