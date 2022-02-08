!  Oct. 2014    minjung.kim@yale.edu
!  This subroutine do FFTs for each bands and update Pmtrx
!
!  Spin is not yet implemented

!subroutine calc_Prrp( psi_v,  psi_N, Bmtrx, Cmtrx, FFTsize, sys, Uvec )
subroutine calc_Brrp( psi_v,  psi_N, Bmtrx, FFTsize, sys, Uvec )

   use constant
   use electronic_structure
   use gw_structure
   
   implicit none
   type(wfn), intent(in) :: psi_v,  psi_N
   type(rank2_mtrx), intent(inout) :: Bmtrx
   !type(rank2_mtrx), intent(inout) :: Cmtrx
   !type(rank2_mtrx)  :: C
   integer, intent(in) :: FFTsize(3)
   type(sysinfo), intent(in) :: sys
   integer, intent(in) :: Uvec(3)

   ! work variables
   complex(dp), dimension(:), allocatable :: a_v, a_c, a_n
   integer :: ndata
   integer :: nv, nc, nb, nbmax, nbmin
   integer :: iv, ic, ib
   integer :: utest
   
   ! Change it ------------------------------
   nv = int ( sum( psi_v%occ ) )
   !nb = int ( sizeof( psi_c%eig(:) ) )/ 8 
   nb = 52
   nc = nb - nv
   !nbmax = 6
   !nbmin = 2
   !-----------------------------------------

   ndata = fftsize(1)*fftsize(2)*fftsize(3)
   utest = abs( Uvec(1) ) + abs( Uvec(2) ) + abs( Uvec(3) )


   !print*, 'number of bands, nv, nc, nbmax:', nb, nv, nc, nbmax
   !allocate( a_v(ndata), a_n(ndata) )
   allocate( a_v(ndata))

!   call update_B( a_v, a_c, a_n, ndata, psi_N%eig(ib), psi_v%eig(iv), psi_c%eig(nb-nc+ic), Bmtrx%C , sys%nk )    
   ! no spin

   ! loop for valence bands

!   do ib = nbmin, nbmax
      !Crrp = B*conjg( a_v(1:ndata) )
!   a_n(1:ndata) = psi_N%cg(1:ndata,ib)
!   if (utest .eq. 0 ) then 
!   continue

!   if (utest .eq. 0 ) then 
!   continue
!   else
!      call modify_Uproc_wfn ( ndata, a_v, Uvec, sys, FFTsize )
!   endif

!   call update_C( a_v, a_c, a_n, ndata, psi_N%eig(ib), psi_v%eig(iv), psi_c%eig(nb-nc+ic), Bmtrx%C , sys%nk )    

     Bmtrx%C = 0.0d0
     do iv = 1, nv

        a_v(1:ndata) = psi_v%cg(1:ndata,iv)
  
        ! change wavefunction for U-process
        if (utest .eq. 0 ) then 
           continue
        else
           call modify_Uproc_wfn ( ndata, a_v, Uvec, sys, FFTsize )
        endif

        call update_B( a_v, ndata, Bmtrx%C)

     enddo

   deallocate(a_v)


   contains 
   subroutine update_B( a_v, ndata, B )

      complex(dp), dimension(ndata), intent(inout) :: a_v
      integer, intent(in) :: ndata
      complex(dp), dimension(ndata,ndata), intent(inout) :: B
   
      ! work variables
      integer :: i, j
    !  real(dp) :: fact

! for n=l 
!for  n =B,
!      f_cv(1:ndata) = a_v(1:ndata) * conjg( a_c(1:ndata) )

!      fact = 2.d0/dble(nk)/(Ev-Ec)
!      diagonal element of the  X
     
     Bmtrx%C = 0.0d0
      do j = 1, ndata
         do i = 1, ndata

            B(i,j) = a_v(i) * conjg( a_v(j) ) + B(i,j)
            
         enddo
      enddo
   end subroutine update_B
   



   subroutine modify_Uproc_wfn( ndata, a_c, Uvec, sys, FFTsize )

      integer, intent(in) :: ndata
      complex(dp), dimension(ndata), intent(inout) :: a_c
      integer, dimension(3), intent(in) :: Uvec
      type(sysinfo), intent(in) :: sys
      integer, dimension(3), intent(in) :: FFTsize
      
      ! work variables
      complex(dp), allocatable :: fftbox(:,:,:), fact(:,:,:)
      integer :: i, j, k, ii
      integer, allocatable :: idx(:,:)
      real(dp) :: a(3,3), b(3,3)
      real(dp) :: rijk(3), G0(3), phase
      
      allocate( idx(3,ndata) )
      allocate( fftbox(fftsize(1),fftsize(2),fftsize(3)) )
      allocate( fact(fftsize(1),fftsize(2),fftsize(3)) )
      
      ! set fftbox index
      call set_3Dbox_index( ndata, FFTsize, idx )
            
      call put_into_fftbox( ndata, a_c, idx, FFTsize, fftbox )
      
      
      a(1:3,1:3) = sys%avec(1:3,1:3)
      b(1:3,1:3) = sys%bvec(1:3,1:3)*2.d0*pi/sys%alat
   
      ! calculate factor to be multiplied to the fftbox
      do k = 1, FFTsize(3)
         do j = 1, FFTsize(2)
            do i = 1, FFTsize(1) 
               do ii = 1, 3
                  rijk(ii) = a(ii,1)*(i-1)/FFTsize(1) + a(ii,2)*(j-1)/FFTsize(2) + &
                             a(ii,3)*(k-1)/FFTsize(3)
                  G0(ii) = b(ii,1)*Uvec(1) + b(ii,2)*Uvec(2) + b(ii,3)*Uvec(3)
               enddo
               G0 = -1.d0*G0
               phase = dot_product( G0, rijk )
               fact(i,j,k) = cmplx( cos(phase), sin(phase) )
            enddo
         enddo
      enddo
      fftbox(1:fftsize(1),1:fftsize(2),1:fftsize(3)) =    &
            fact(1:fftsize(1),1:fftsize(2),1:fftsize(3))*    &
            fftbox(1:fftsize(1),1:fftsize(2),1:fftsize(3)) 
         
      
      call box_to_array( FFTsize, fftbox, a_c )
   
      deallocate( idx, fftbox, fact )
   end subroutine modify_Uproc_wfn
   
   
end subroutine! - End of Calc_Crrp subroutine



   
!###########################################################################
!
!                 FFT : P(r,r') -> P(G,G')
!                 #FFT = 2*ndata
!
subroutine B_r_to_g( iq, sign_fft, B, FFTsize )

   use constant
   use gw_structure
   implicit none
   integer, intent(in) :: iq, sign_fft
   integer, intent(in) :: FFTsize(3)
   type(rank2_mtrx), intent(inout) :: B
   !type(rank2_mtrx), intent(inout) :: C
   
   ! work variables
   complex(dp), allocatable :: a_r(:), fftbox(:,:,:)
   integer :: icol, irow, colrow
   integer :: ndata
   integer, allocatable :: idx(:,:)
   integer, parameter :: one = 1, minus_one = -1
   
   ndata = FFTsize(1)*FFTsize(2)*FFTsize(3)
   allocate( idx(3, ndata) )
   allocate( a_r(ndata), fftbox(FFTsize(1),FFTsize(2),FFTsize(3)) )
   
   ! set index
   call set_3Dbox_index( ndata, FFTsize, idx )
   
   
   ! Time to FFT
   ! standard: -1, and +1
!   do colrow = 1, 2
      
   do icol = 1, ndata
      
         a_r( 1:ndata ) = B%C( 1:ndata, icol )
         ! put 1D array to 3D box
         call put_into_fftbox( ndata, a_r, idx, FFTsize, fftbox ) 
         ! FFT WARNING: SIGN CHANGE
         !if ( colrow .eq. 1 ) then
         call do_fft( fftbox, FFTsize, sign_fft)
          !  if (iq .ne. 1) call do_fft( fftbox, FFTsize, minus_one)
        ! elseif ( colrow .eq. 2 ) then  
         !   if (iq .eq. 1) call do_fft( fftbox, FFTsize, minus_one)
         !   if (iq .ne. 1) call do_fft( fftbox, FFTsize, one)
         !endif
         ! save FFTed values to a_r ( but its contents are g-space values )
         call box_to_array( FFTsize, fftbox, a_r )
          ! exchange column to FFTed values
         B%C( 1:ndata, icol ) = a_r( 1:ndata )
    enddo

    deallocate(idx)
    deallocate(a_r)
      
!   enddo

end subroutine
