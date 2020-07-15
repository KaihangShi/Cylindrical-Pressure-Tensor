! ==========================================================
! This subroutine is used to set up Ewald sum parameters
! Modified for cylindrical pressure tensor on 9/26/2019
! ==========================================================

Subroutine set_ewld(ewld_prec,q_o,q_h,n_mol_sites,n_mol_tot,iframe)

Use global

IMPLICIT NONE

! Passed 
Integer :: iframe, n_mol_sites, n_mol_tot
Double Precision :: q_o, q_h
  Double Precision :: ewld_prec

! Local
INTEGER :: n_k,ksq,kx,ky,kz,k, imol
INTEGER :: itype,ichrg,isite
DOUBLE PRECISION :: rksq,rkx,rky,rkz
DOUBLE PRECISION :: b
DOUBLE PRECISION,DIMENSION(3) :: two_Pibox
Double Precision :: dzis


DOUBLE PRECISION :: testeng


two_Pibox(:) = two_Pi/box(:,iframe)


! Using heuristic equation to calculate the kmax compatible with LAMMPS (version 12Dec18)
! ewld_prec should be equal to the 'xxx' in 'kspace_modify force xxx' in LAMMPS setting
! reference: https://lammps.sandia.gov/threads/msg18268.html
kxmax = CEILING( alpha * box(1,iframe)/Pi * dSQRT(-log(ewld_prec)))
kymax = CEILING( alpha * box(2,iframe)/Pi * dSQRT(-log(ewld_prec)))
kzmax = CEILING( alpha * box(3,iframe)/Pi * dSQRT(-log(ewld_prec)))

! Calculate ksq_max
ksq_max = MAX(kxmax,kymax,kzmax)**2 + 2


! --- Update on 9/26/2019: now we fix cutoff and adjust kmax
! Take ewald_fix_kmax scheme follows the same one in Towhee
! Calculate corresponding convergence parameter in the real space
!alpha = kalp/MINVAL(box(:,iframe))

! Set cut_off raidus in real space to half of the box length
!rcelect = MINVAL(box(:,iframe))/2.0d0
!rcelectsq = rcelect**2

! Setup reciprocal vectors
! Calculate the denominator for the k-space exponential
b = 1.0d0/(4.0d0*alpha**2)


! Calculate reciprocal vector 
Do imol = 1, n_mol_tot
   
  ! Loop over the sites on the molecule
  DO isite = 1, n_mol_sites

    ! Calculate the reciprocal vector components for k = 0
    eikx(imol,isite,0) = (1.0d0, 0.0d0)
    eiky(imol,isite,0) = (1.0d0, 0.0d0) 
    eikz(imol,isite,0) = (1.0d0, 0.0d0)

    ! Calculate the reciprocal vector components for k = 1
    eikx(imol,isite,1) = &
      & dCMPLX(dCOS(two_Pibox(1)*rx_s(isite,imol,iframe)),dSIN(two_Pibox(1)*rx_s(isite,imol,iframe)))

    eiky(imol,isite,1) = &
      & dCMPLX(dCOS(two_Pibox(2)*ry_s(isite,imol,iframe)),dSIN(two_Pibox(2)*ry_s(isite,imol,iframe)))

    eikz(imol,isite,1) = &
      & dCMPLX(dCOS(two_Pibox(3)*rz_s(isite,imol,iframe)),dSIN(two_Pibox(3)*rz_s(isite,imol,iframe)))

    ! Calculate the reciprocal vector components for k = -1 
    eikx(imol,isite,-1) = dCONJG(eikx(imol,isite,1))
    eiky(imol,isite,-1) = dCONJG(eiky(imol,isite,1))
    eikz(imol,isite,-1) = dCONJG(eikz(imol,isite,1))

  ENDDO

End Do

! Calculate the remaining reciprocal vector components by recursion
DO imol = 1,n_mol_tot

  DO isite = 1,n_mol_sites

    ! x-component 
    Do k = 2, kxmax
      eikx(imol,isite,k) = eikx(imol,isite,k-1)*eikx(imol,isite,1) 
      eikx(imol,isite,-k) = dCONJG(eikx(imol,isite,k))
    End Do

    ! y-component
    Do k = 2, kymax
      eiky(imol,isite,k) = eiky(imol,isite,k-1)*eiky(imol,isite,1)
      eiky(imol,isite,-k) = dCONJG(eiky(imol,isite,k))
    End Do
    
    ! z-component
    Do k = 2, kzmax
      eikz(imol,isite,k) = eikz(imol,isite,k-1)*eikz(imol,isite,1)
      eikz(imol,isite,-k) = dCONJG(eikz(imol,isite,k))
    End Do
      
  ENDDO
ENDDO




! Initialize the number of vectors
n_k = 0

! Loop over all the vectors
! vectors follow Stanford Ewald notes 
DO kx = -kxmax,kxmax
  rkx = two_Pibox(1)*DBLE(kx)
  DO ky = -kymax,kymax
    rky = two_Pibox(2)*DBLE(ky)
    DO kz = -kzmax,kzmax
      rkz = two_Pibox(3)*DBLE(kz)

      ! Calculate the square of the vector
      ksq = kx**2 + ky**2 + kz**2

      ! Check the bounds of the vector, the range is a sphere
      IF((ksq .LT. ksq_max) .AND. (ksq .NE. 0)) THEN

        ! Update the number of vectors
        n_k = n_k + 1

        If (n_k .gt. maxk) Then
          Write(5,*) 'FATAL ERROR: TOO SMALL OF maxk'
          STOP
        End If

        ! Calculate rksq
        rksq = rkx**2 + rky**2 + rkz**2

        ! Calculate the the prefactor of the first term 
        k_vec(1,n_k) = dEXP(-b*rksq)/rksq * (1.0d0 - 2.0*rkz**2*(1.0/rksq + b))
        k_vec(2,n_k) = dEXP(-b*rksq)/rksq * rkz


        ! Initialize the structure factor
        skewld(n_k) = (0.0d0, 0.0d0)
        skewlds(n_k) = (0.0d0, 0.0d0)

        ! Loop over the molecules
        DO imol = 1,n_mol_tot

          ! Loop over imol's sites
          DO isite = 1, n_mol_sites

            dzis = rz_s(isite,imol,iframe) - rz(imol,iframe)
            dzis = dzis - dNINT(dzis/box(3,iframe))*box(3,iframe)

            If (isite .eq. 1) Then

              ! Oxygen atom in water
              skewld(n_k) = skewld(n_k) + q_o*&
                  & eikx(imol,isite,kx)*eiky(imol,isite,ky)*eikz(imol,isite,kz)

              skewlds(n_k) = skewlds(n_k) + q_o*dzis*&
                  & eikx(imol,isite,kx)*eiky(imol,isite,ky)*eikz(imol,isite,kz)

            Else 
              ! Hydrogen atom in water
              skewld(n_k) = skewld(n_k) + q_h*&
                  & eikx(imol,isite,kx)*eiky(imol,isite,ky)*eikz(imol,isite,kz)

              skewlds(n_k) = skewlds(n_k) + q_h*dzis*&
                  & eikx(imol,isite,kx)*eiky(imol,isite,ky)*eikz(imol,isite,kz)

            End If
          ! End loop over sites
          ENDDO
        
        ! End loop over molecules
        ENDDO

      ENDIF

    ENDDO
  ENDDO
ENDDO


Return

End Subroutine 
