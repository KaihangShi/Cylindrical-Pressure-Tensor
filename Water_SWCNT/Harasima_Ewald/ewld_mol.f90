! ==========================================================
! This subroutine is used to calculate Fourier space 
! contribution to the axial pressure in cylindrical coordinates
! following Harasima definition
! 
! Created on 9/24/2019 by Kaihang Shi
!
! ==========================================================

	  Subroutine ewld_mol(q_o,q_h,n_mol_sites,n_mol_tot,imol,iframe,ptzf)

	  Use global

	  IMPLICIT NONE

	  ! Passed
	  Double Precision :: q_o, q_h
	  Integer :: imol, iframe, n_mol_sites, n_mol_tot
	  Double Precision :: ptzf

	  ! Local
	  INTEGER :: n_k,ksq,kx,ky,kz
	  INTEGER :: isite, jmol, jsite
	  !DOUBLE PRECISION,DIMENSION(3) :: two_Pibox
	  Complex*16 :: img
	  Complex*16, Dimension(n_mol_sites) :: eikis
	  Double Precision :: ptzc_1, ptzc_2
	  Double Precision :: dzis
	  Double Precision :: skimol

	  
	  ! Initialize variables
	  ptzf = 0.0d0
	  ptzc_1 = 0.0d0
	  ptzc_2 = 0.0d0
	  img = (0.0,1.0)
	  !two_Pibox(:) = two_Pi/box(:,iframe)


      
	  ! reciprocal contribution to the pressure
	  ! Initialze the vector
	  n_k = 0

	  DO kx = -kxmax,kxmax
	    DO  ky = -kymax,kymax
	      DO kz = -kzmax,kzmax

	        ! Calculate the square of the vector
	        ksq = kx**2 + ky**2 + kz**2

	        ! Check the bounds of the vector
	        IF((ksq .LT. ksq_max) .AND. (ksq .NE. 0)) THEN

	          ! Update the number of vectors
	          n_k = n_k + 1

	          ! Calculate exp(-ih*r_ia)
	          Do isite = 1, n_mol_sites

	          	eikis(isite) = eikx(imol,isite,-kx)*eiky(imol,isite,-ky)*eikz(imol,isite,-kz)

	          End Do

	          ! --- first part of the Fourier contribution ---! 
	          skimol = 0.0d0
	          ! Loop over imol's sites
	          Do isite = 1, n_mol_sites

	          	If (isite .eq. 1) Then
	          		! Oxygen atom in water
	          		skimol = skimol + q_o * REALPART( eikis(isite)*skewld(n_k) )
	          		
	          	Else
	          		! Hydrogen atom in water
	          		skimol = skimol + q_h * REALPART( eikis(isite)*skewld(n_k) )

	          	End If
	          ! End Loop over sites
	          End Do
	    	  
	    	  ! Multiply the prefactor
	          ptzc_1 = ptzc_1 + skimol * k_vec(1,n_k)


	          ! --- second part of the Fourier contribution ---! 
	          skimol = 0.0d0 
	          Do isite = 1, n_mol_sites

	          	dzis = rz_s(isite,imol,iframe) - rz(imol,iframe)
	          	dzis = dzis - dNINT(dzis/box(3,iframe))*box(3,iframe)
	          
          		If (isite .eq. 1) Then
          			! Oxygen atom in imol
          			skimol = skimol + q_o * REALPART( img*dzis*eikis(isite)*skewld(n_k) - img*eikis(isite)*skewlds(n_k) )

          		Else 

          			skimol = skimol + q_h * REALPART( img*dzis*eikis(isite)*skewld(n_k) - img*eikis(isite)*skewlds(n_k) )
          			
          		End If

          	  ! End isite
	          End Do

	          ! Update
	          ptzc_2 = ptzc_2 + skimol * k_vec(2,n_k)



	        ENDIF

	      ENDDO
	    ENDDO
	  ENDDO


	  ! Update total contribution from reciprocal space
	  ptzf = (- ptzc_1 + ptzc_2)



	  Return
	  
	  End Subroutine 















