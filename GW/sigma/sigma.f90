!  Sigma_x main program
!
!  Sep. 2014  minjung.kim@yale.edu
!  July 2015  subhasish.mandal@yale.edu

program sigma

   use constant
   use usrinput
   use electronic_structure
   use gw_structure

   implicit none
   type(input) :: inp
   type(sysinfo) :: sys, sys_shft
   type(wfstruc) :: psi, psi_R, psi_shft, psi_R_shft
   type(wfn) :: psi_v,  psi_N 
   ! k is in cartesian unit, k_ is in crystal coordinate
   type(kptinfo) :: k, q, k_, & 
                    k_shft
   !type(polarizability) ::B
   !type(rank2_mtrx) :: Bmtrx, Cmtrx, Dmtrx
   type(rank2_mtrx) :: Brrp, Crrp
   type(rank2_mtrx) :: Drrp, Crgp
   !real(dp)  :: coul
   !real(dp) :: coulb
   !real(dp), dimension(ndata) :: coulb
   !real(dp), allocatable ::  sig_x, sig_sum
   real(dp) ::  sig_x, sig_sum
   real(dp) ::  Ha2ev
   type(gspace) :: g_pol     ! Polarizability g space
   !real(dp) :: qvec ! q vector
   real(dp) :: qvec(3) ! q vector
   ! determine FFT grid. This number is multiplied to the minimum cutoff radius
   integer :: FFTsize(3)
   integer :: nk, nq, iq, ik,  g
   integer :: icm, idm
   integer :: ikq ! kpt index for k+q vector
   !integer :: nfft ! FFTsize(1)*FFTsize(2)*FFTsize(3)
   integer :: ndata ! number of g vector in epsilon matrix   
   !integer, allocatable :: ndata, nfft ! number of g vector in epsilon matrix   
   integer, allocatable :: gidx(:) ! link G vectors between polarizability and epsilon
   integer :: Uvec(3) ! Umklapp vector
   !integer :: gvec(3,:) ! g vectors
   integer, allocatable :: gvec(:,:) ! g vectors
   integer:: i,j, idata, ii
   !complex(dp) :: Dgg
   complex(dp), dimension(:), allocatable ::a_n
   real(dp), dimension(:), allocatable ::coulb
   integer :: nv, nc, nb, nbmax, nbmin
   integer :: iv, ic, ib
!*****************************************
integer, allocatable :: gv(:,:), gpv(:,:)
integer :: ig, igp, test, kk
real(dp) :: g2, gp2
!*****************************************
!*******************
!******Band info*******
   nv = 4
   nb = 52
   nc = nb - nv
   nbmax = 1
   nbmin = 1
!********************
!********************
   Ha2ev=27.211382543519
!********************
!********************
   ! Read input values
   call read_input( inp )
   ! Read wavefunction data and save into the memory
   call read_wfn( inp%wfname, sys, psi, k )

   ! check_dim is for testing.
   ! call check_dim( psi )

   ! This routine gives you the size of the FFT box
   call set_FFTsize( psi, FFTsize )

   nk = k%nk
   sys%nk = nk
   nq = nk  ! nq = nk
   q%nk = nq
   print*,'check Nq', nq, nk
   !stop

   ! retrieve k vectors in reciprocal lattice basis from cartesian basis
   ! (unit is 2*pi/a for both)
   ! relevant for QE output
   print *, 'Calling cartesian_to_crystal'
   call cartesian_to_crystal( sys, k, k_ )

   ! Let's get k and q vector information
   print *, 'get_qvec'
   call get_qvec( k_, q )
   !print *, 


   !!print*, 'Calculate sigma at q:', q%vec(1:3,iq)
   ! Do FFT for all k points and bands
   ! psi (in gspace) will be removed from the memory
   !here psi==psi_G as input and psi_R as output
   print *, 'FFT_all_bands'
   call FFT_all_bands( psi, psi_R, FFTsize, sys )

   ! Polarizability matrix construction
  ! B%nq = nq
  ! call create_Pmtrx_struc( B, FFTsize )
  ! print *, 'Allocating a pile'
  !allocate(coulb(ndata))
   !allocate(coulb(FFTsize(1)*FFTsize(2)*FFTsize(3)))

   ndata = FFTsize(1)*FFTsize(2)*FFTsize(3)
   allocate(a_n(FFTsize(1)*FFTsize(2)*FFTsize(3)))
   allocate(Brrp%C(FFTsize(1)*FFTsize(2)*FFTsize(3),FFTsize(1)*FFTsize(2)*FFTsize(3)) )
   allocate(Crrp%C(FFTsize(1)*FFTsize(2)*FFTsize(3),FFTsize(1)*FFTsize(2)*FFTsize(3)) )
   allocate(Drrp%C(FFTsize(1)*FFTsize(2)*FFTsize(3),FFTsize(1)*FFTsize(2)*FFTsize(3)) )
   !allocate( g_pol%gvec(3,ndata) , gidx(ndata) )
!!   allocate(gvec(3,ndata))
!!   allocate( gvec(3,ndata) )
   allocate(coulb(ndata))
   allocate( gvec(3,ndata) , gidx(ndata) )
!!   g_pol%ng = ndata
!!
!***************************************************************** 
   !sig_sum = sig_x
  ! CASE2: q \= 0
!   do iq = 2, nq

!           do ik = 1, nk

                    ! let's calculate which k vector is the same as k+q
!           call get_kq_index( k_, q, ik, iq, ikq, Uvec )


!           enddo
!     enddo
       ! remove psi_R 
!   do ik = 1, nk 
!        deallocate( psi_R%wk(ik)%cg ) 
!   enddo



  do ik = 1, nk

    do  ib = nbmin, nbmax


       do iq = 2, nq

!       do g = 1, ndata

         call fftidx_to_gidx( FFTsize, ndata, gvec )

         ! let's calculate which k vector is the same as k+q
         call get_kq_index( k_, q, ik, iq, ikq, Uvec )

!!         call calc_sig( Drrp, sig_x, ndata, sys, nq )  
         !  Getting B-matrix.


         call calc_Brrp(psi_R%wk(ikq), psi_R%wk(ikq), Brrp, FFTsize, sys, Uvec )
         ndata = FFTsize(1)*FFTsize(2)*FFTsize(3)
         a_n(1:ndata) = psi_R%wk(ik)%cg(1:ndata,ib)

         do icm =1, ndata
          Crrp%C(icm, 1:ndata) = Brrp%C(icm, 1:ndata)*(a_n(1:ndata))
         enddo



!         print*, 'FFT (Crrp -> Crgp)'
         call B_r_to_g( iq, 1, Crrp, FFTsize )

         Crrp%C = Crrp%C * sys%vol / ndata
         !Now let's compute D
         do idm = 1, ndata
           Drrp%C(1:ndata,idm) = Crrp%C(1:ndata,idm)*conjg(a_n(1:ndata))
         enddo
         Drrp%C(1:ndata,1:ndata) = Transpose( Drrp%C(1:ndata,1:ndata) )
!!
         call B_r_to_g( iq, -1, Drrp, FFTsize )

         Drrp%C = Drrp%C * sys%vol / ndata


         !print *, 'Now calculating sig'
!         call calc_sig( Drrp, sig_x, ndata, sys, nq, coulb, g ) ! qvec , gvec, sys) 
   !      call calc_coulb( sys, ndata, qvec, gvec, coulb, nq )
!         print*, 'Which band you are on : ', ib, coulb
!       enddo
!    print*, 'Which k-point you are on : ', k%vec(1:3,ik),  ik
    !!print*, 'Which band you are on : ', ib, coulb
    !print*, 'Which gvec ? ', %gs(iq)%gvec(:,i)
!    print*, "sigma:", sig_x*Ha2ev



        call calc_coulb( sys, ndata, qvec, gvec, coulb, nq )

!    enddo

        !print*, 'Which k-point you are on : ', ik, ndata
!  enddo
   ! remove psi_R 
 !  do ik = 1, nk 
 !     deallocate( psi_R%wk(ikq)%cg ) 
 !  enddo


   !print*, "sigma:", sig_sum
!   print*, "ndata:", ndata
!   print*, "ndata:", nfft
 ! print*, "nk:", nk, nq



 !*****************************************
print*, 'print out coulmb'
open(1,file='COUL_1.dat',form='formatted',status='unknown')
!write(1,*) '1728'
write(1,*) 'Which k-point you are on : ', ik
write(1,*) '******************************* '
write(1,*) '******************************* '
write(1,*) 'start for in =    ', ib


write(1,*) ' Vcoulb n1   ig   ga  gb  gc      '




do j = 1, ndata 
 write(1, '(f15.8, 510i8)')coulb(j), ib, j, gvec(:,j)!
enddo
!!******************************************

enddo
enddo
enddo
close(1)



deallocate(a_n)
deallocate(Brrp%C )
deallocate(Crrp%C)
deallocate(Drrp%C)


!enddo 



end program
