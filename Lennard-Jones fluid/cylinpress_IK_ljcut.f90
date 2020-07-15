! ==========================================================================
!  This version is for calculating Irving-Kirkwood contour local pressure tensor (full)
!  in cylindrical coordinates
!  Input: XYZ coordiantes from MD simulation
! 
!  Note: 
!  1) No periodic boundary condition was applied for pressure calculation 
!     so rden_cut should be chosen properly
!  2) No impulsive contribution to the pressure in MD simulation
!  3) Cylindrical geometry is defined at the center of the box and 
!     axial direction is aligned with the z-axis of the box
!  
!  Author: Kaihang Shi
!  Last update: September 9, 2019
! ==========================================================================

PROGRAM virialpress_cylinder

IMPLICIT NONE

! ---------------- Define variables ---------------------------!
! Control parameters 
Integer, Parameter :: n_frame = 49999               ! Total simulation frames for sampling
! System parameter 
Double Precision, Parameter :: temp = 100.0d0       ! Temperature [Kelvin]
Integer, Parameter :: n_mol_types = 1               ! Number of molecule types  
Integer, Parameter :: n_mol_tot = 4000              ! Total number of LJ particles in the system
! Intermolecular parameter 
Double Precision, Parameter :: sigma = 3.405d0      ! LJ size paramter (sigma) [Angstrom] 
Double Precision, Parameter :: epsilonkb = 119.8d0  ! LJ energy parameter (epsilon/kB) [Kelvin]
Double Precision, Parameter :: r_cut = 10.215d0     ! Cutoff radius [Angstrom]
Double Precision, Parameter :: mol_mass = 39.948d0  ! Mass [g/mol]
! Parameters for radial density and pressure calculation
Double Precision, Parameter :: rden_cut = 28.0d0   ! Cutoff radius for density/pressure calculation [A]
Integer, Parameter :: rden_bins = 1000             ! Density/pressure profile resolution
Double Precision, Parameter :: zlo = -10.0         ! lower bound in axial direction for pressure/density calculation [A]
Double Precision, Parameter :: zhi = 10.0          ! upper bound in axial direction for pressure/density calculation [A]

! Box dimension
Double Precision, Dimension(:,:), Allocatable :: box
! Read in parameter
Character(len=128):: dummy
Character(len=128) :: site_name
Double Precision, Dimension(:,:), Allocatable :: boxangle
Integer :: frame_id
Integer :: n_sites_tot
! Coordinates
Double Precision, Dimension(:,:), Allocatable :: rx_s, ry_s, rz_s
! delta value in equation for Harasima contour use
Double Precision :: delrr
! Statistics
! 1-fluid-wall; 2-fluid-fluid; 3-total
Double Precision, Dimension(:,:), Allocatable :: virialpress_cylin_pnr
Double Precision, Dimension(:,:), Allocatable :: virialpress_cylin_ptt
Double Precision, Dimension(:,:), Allocatable :: virialpress_cylin_ptz
Double Precision, Dimension(:,:), Allocatable :: virialpnravg
Double Precision, Dimension(:,:), Allocatable :: virialpttavg
Double Precision, Dimension(:,:), Allocatable :: virialptzavg
Double Precision, Dimension(:), Allocatable :: prkin
! r-density statistics (Density profile of each molecule type in the cylindrical system)
! Cutoff radius for density and cylindrical pressure tensor calculation
! limit in radial direction
Double Precision :: rden_lim
! limit in axial direction
Double Precision :: zlo_lim, zhi_lim
! r-density length in r direction for each bin
Double Precision :: rden_dr, rden_drsq
! Statistics 
Double Precision, Dimension(:,:), Allocatable :: rdenavg
Double Precision, Dimension(:,:,:), Allocatable :: rdenblk

!Variables for pressure calculation
INTEGER :: imol, jmol, ibin, iframe, idirec, ierr, isol
INTEGER :: isite, jsite, itype, jtype, isitetype, jsitetype
!DOUBLE PRECISION :: rxi, ryi, rzi, rij, rxj, ryj, rzj
DOUBLE PRECISION :: rxis, ryis, rzis, rxjs, ryjs, rzjs
DOUBLE PRECISION :: sr3, sr6, sr12, phitz
DOUBLE PRECISION :: pnrc, ptzc, dudr, pttc
DOUBLE PRECISION :: rijs, xijs, yijs, zijs
DOUBLE PRECISION :: xijssq, yijssq, zijssq
Double Precision :: clri, clrj, clrij, posr, clrk, zlk
Double Precision :: clr, clrsq, rden_dvol
Double Precision :: aa, bb, cc, dcrt
Double Precision, Dimension(2) :: alpk

! Constants
Double Precision, Parameter :: Pi = 3.141592653589d0
Double Precision, parameter :: two_Pi = 2.0d0*Pi
Double Precision, Parameter :: Kb = 1.3806488d-23   ! [J/K] 
Double Precision, Parameter :: Na = 6.02214129d23  ! [1/mol] */
Double Precision, Parameter :: h = 6.62606957d-34   ! [J*s] = [kg*m^2/s]
Double Precision, Parameter :: R = Kb*Na            ! [J/(mol*K)] */
Double Precision, Parameter :: EETOK = 167100.9558738398d0 ! EETOK = 1/(4*pi*epsilon*kb) * e^2 * 10^10
Double Precision, Parameter :: PVTOK = 0.7242971565d-02   ! PVTOK = 1.0d5*1.0d-30/kb , [bar]*[A^3] to [K]
Double Precision, Parameter :: PCOEFF = Kb*1.0d25        ! [K]/[A]^3 to [bar]

! -------------------------------------------------------------------------------------------------! 

! -------- Allocate variables ----------!
Allocate(box(3,n_frame), STAT=ierr)
Allocate(boxangle(3,n_frame), STAT=ierr)
Allocate(rx_s(n_mol_tot,n_frame), STAT=ierr)
Allocate(ry_s(n_mol_tot,n_frame), STAT=ierr)
Allocate(rz_s(n_mol_tot,n_frame), STAT=ierr)
Allocate(virialpress_cylin_pnr(2,rden_bins), STAT=ierr) 
Allocate(virialpress_cylin_ptt(2,rden_bins), STAT=ierr) 
Allocate(virialpress_cylin_ptz(2,rden_bins), STAT=ierr) 
Allocate(virialpnravg(3,rden_bins), STAT=ierr) 
Allocate(virialpttavg(3,rden_bins), STAT=ierr) 
Allocate(virialptzavg(3,rden_bins), STAT=ierr) 
Allocate(prkin(rden_bins), STAT=ierr) 
Allocate(rdenavg(rden_bins,n_mol_types), STAT=ierr) 
Allocate(rdenblk(rden_bins,n_mol_types,n_frame), STAT=ierr) 

! -- Initialize parameter -----!
virialpress_cylin_ptt = 0.0d0
virialpress_cylin_ptz = 0.0d0
virialpress_cylin_pnr = 0.0d0
virialpnravg = 0.0d0
virialpttavg = 0.0d0
virialptzavg = 0.0d0
prkin = 0.0d0
rdenavg = 0.0d0
rdenblk = 0.0d0
rden_dr = rden_cut/DBLE(rden_bins)
rden_drsq = rden_dr**2
delrr = rden_dr/2.0d0
rden_lim = rden_cut + r_cut
zlo_lim = zlo - r_cut
zhi_lim = zhi + r_cut

! --------- Write header ------------!
Write(*,'(A)') '==========================================================  '
Write(*,'(A)') ' This program is used to calculate the pressure tensor for  ' 
Write(*,'(A)') ' cylindrical geometry from virial route                     '
Write(*,'(A)') ' using the Irving-Kirkwood (IK) definition of contour'
Write(*,'(A)') '                                                            '  
Write(*,'(A)') ' This version for simple LJ site molecule with LJ simple cut'
Write(*,'(A)') ' force. No impulsive contribution for MD simulation.        '
Write(*,'(A)') ' ========================================================== '
Write(*,'(A)') '                                                            '

! ---- Read in coordinates from my XYZ format file ------!
Open(1,File='../coor.xyz',Status='Old',Access='Sequential',Action= 'Read')
! Loop over all frames 
Do iframe = 1, n_frame

  ! Read in header
  Read(1,*) (box(idirec,iframe),idirec=1,3), (boxangle(idirec,iframe),idirec=1,3), frame_id
  Read(1,*) n_sites_tot

  ! Make sure this run is sensible
  If (n_sites_tot .NE. n_mol_tot) Then
    Write(*,*) 'n_sites_tot is not equal n_mol_tot! Check code!'
    STOP
  End If
  If (rden_lim .GT. Min(box(1,iframe),box(2,iframe))/2.0d0) Then
    Write(*,*) 'rden_cut or r_cut is too large for this box size!'
    STOP
  End If

  Read(1,'(A)') dummy

  ! Start reading in xyz coordinates
  Do isite = 1, n_sites_tot
    Read(1,*) site_name, rx_s(isite,iframe), ry_s(isite,iframe), rz_s(isite,iframe)
  End Do

  ! Read in space according to the 'coor.xyz' format 
  Read(1,'(A)') dummy

End Do 
! Close file
CLOSE(1)
Write(*,*) 'Finished reading coordinates and box information from coor.xyz file!'
! unit=6 is screen
Call FLUSH(6)




! --------------- Start postprocessing ----------------------!
Do iframe = 1, n_frame
  
  ! ----------- Sampling r-density ------------!
  ! Loop over all sites/molecules
  Do isite = 1, n_sites_tot

    ! Only sampling over particles within [zlo,zhi]
    if((rz_s(isite,iframe) .ge. zlo) .AND. (rz_s(isite,iframe) .le. zhi)) Then

      ! Calculate r-distance of isite in cylindrical coordiantes
      clrsq = rx_s(isite,iframe)**2 + ry_s(isite,iframe)**2
      clr = dSQRT(clrsq)

      If (clr .LT. rden_cut) Then
        ! Calculate ibin number
        ibin = FLOOR(clr/rden_dr) + 1

        ! Accumulate number 
        rdenblk(ibin,1,iframe) = rdenblk(ibin,1,iframe) + 1.0d0
      End If

  ENDIF

  ! End loop over all molecules 
  End Do

  ! Loop over bins
  Do ibin = 1, rden_bins
    ! Loop over molecule types
    Do itype = 1, n_mol_types

      ! Calculate volume for each bin (changing in NPT ensemble)
      rden_dvol = Pi*DBLE(2*ibin-1)*(zhi-zlo)*rden_drsq

      ! Convert to number density (1/A^3)
      rdenavg(ibin,itype) = rdenavg(ibin,itype) + rdenblk(ibin,itype,iframe)/rden_dvol

    End Do          
  End Do  
  

  ! ----- Sampling cylindrical pressure tensor -------!   
  ! Loop over sites/molecules
  DO isite = 1, n_sites_tot

    rxis = rx_s(isite,iframe) 
    ryis = ry_s(isite,iframe) 
    rzis = rz_s(isite,iframe)

    ! Only sampling particles within [zlo,zhi]
    if((rzis .lt. zlo_lim) .or. (rzis .gt. zhi_lim)) CYCLE

    ! Calculate R-distance of site i in cylindrical coordiantes
    clri = dSQRT(rxis**2 + ryis**2)

    ! Determine if carrying on (rden_lim = rden_cut + r_cut)
    If (clri .GT. rden_lim) CYCLE

    !Fluid-fluid interaction, avoiding double count of i,j and j,i so just
    !loop by using i<j     
    IF (isite .LT. n_sites_tot) THEN
      DO jsite = isite+1, n_sites_tot

        rxjs = rx_s(jsite,iframe) 
        ryjs = ry_s(jsite,iframe) 
        rzjs = rz_s(jsite,iframe)

        ! Only sampling particles within [zlo,zhi]
        if((rzjs .lt. zlo_lim) .or. (rzjs .gt. zhi_lim)) CYCLE
        if((rzis .gt. zhi) .and. (rzjs .gt. zhi)) CYCLE
        if((rzis .lt. zlo) .and. (rzjs .lt. zlo)) CYCLE

        ! Calculate R-distance of site j 
        clrj = dSQRT(rxjs**2+ryjs**2)

        ! Determine if carrying on
        If (clrj .GT. rden_lim) CYCLE

        !Calculate vector between sites
        xijs = rxjs - rxis
        yijs = ryjs - ryis
        zijs = rzjs - rzis

        ! Apply minimum image convention
        ! Only apply to z-direction
        !zijs = zijs - dNINT(zijs/box(3,iframe))*box(3,iframe)

        !Square the values
        xijssq=xijs*xijs
        yijssq=yijs*yijs
        zijssq=zijs*zijs

        !Calculate distance between i and j sites
        rijs=dSQRT(xijssq+yijssq+zijssq)

        ! Check with force cutoff
        If (rijs .GT. r_cut) CYCLE

        ! Calculate coefficient in ax^2+bx+c form for R=R_l use
        aa = xijssq + yijssq
        bb = 2.0*(rxis*xijs + ryis*yijs)

        !Calculate 12-6LJ force (based on sites)
        sr3=(sigma/rijs)**3
        sr6=sr3*sr3
        sr12=sr6*sr6
        phitz = 24.0*epsilonkb/rijs
        dudr = phitz*(sr6 - 2.0d0*sr12)

        ! == Irving-Kirkwood definition ==
        !Loop through the bins
        DO ibin= 1, rden_bins

          ! Get r-distance of ibin in cylindrical system
          posr = (DBLE(ibin)-0.5d0)*rden_dr

          ! Calculate alpha_k
          cc = rxis**2 + ryis**2 - posr**2
          ! Discriminant 
          dcrt = bb**2 - 4.0*aa*cc
          if(dcrt .LE. 0.0d0) CYCLE
          ! two solutions for alpha
          alpk(1) = (-bb + dSQRT(dcrt))/(2.0*aa)
          alpk(2) = (-bb - dSQRT(dcrt))/(2.0*aa)

          if((alpk(1) .LT. 0.0d0 ) .and. (alpk(2) .GT. 1.0d0)) CYCLE
          if((alpk(2) .LT. 0.0d0 ) .and. (alpk(1) .GT. 1.0d0)) CYCLE

          ! -------- P_RR & P_thetatheta ---------
          ! i, j both in [zlo,zhi]
          If ((rzis .ge. zlo) .and. (rzis .le. zhi) .and. (rzjs .ge. zlo) .and. (rzjs .le. zhi)) Then

            ! Loop over alpha 
            Do isol = 1,2 

              ! contribute
              if((alpk(isol) .ge. 0.0d0) .and. (alpk(isol) .le. 1.0d0)) Then
                ! R_l at alpk
                clrk = dSQRT( (rxis + alpk(isol)*xijs)**2 + (ryis + alpk(isol)*yijs)**2 )

                ! Normal pressure component
                pnrc = dABS(xijs * (rxis + alpk(isol)*xijs)/clrk + yijs * (ryis + alpk(isol)*yijs)/clrk) * dudr/rijs

                !Update normal pressure from the contribution of fluid-fluid interactions (2)
                virialpress_cylin_pnr(2,ibin) = virialpress_cylin_pnr(2,ibin) + pnrc

                ! azimuthal pressure component
                pttc = (-xijs * (ryis + alpk(isol)*yijs)/clrk + yijs * (rxis + alpk(isol)*xijs)/clrk)**2/ &
                  & dABS(xijs * (rxis + alpk(isol)*xijs)/clrk + yijs * (ryis + alpk(isol)*yijs)/clrk) * &
                  & dudr/rijs

                !Update normal pressure from the contribution of fluid-fluid interactions (2)
                virialpress_cylin_ptt(2,ibin) = virialpress_cylin_ptt(2,ibin) + pttc

              End if
            ! End looping over alpha
            End Do

          ! i and j either one in [zlo,zhi]
          Else

            ! Loop over alpha 
            Do isol = 1,2 

              ! contribute
              if((alpk(isol) .ge. 0.0d0) .and. (alpk(isol) .le. 1.0d0)) Then
                ! R_l at alpk
                clrk = dSQRT( (rxis + alpk(isol)*xijs)**2 + (ryis + alpk(isol)*yijs)**2 )

                ! Equally contribute half
                pnrc = 0.5d0*dABS(xijs * (rxis + alpk(isol)*xijs)/clrk + yijs * (ryis + alpk(isol)*yijs)/clrk) * dudr/rijs

                !Update normal pressure from the contribution of fluid-fluid interactions (2)
                virialpress_cylin_pnr(2,ibin) = virialpress_cylin_pnr(2,ibin) + pnrc

                ! Normal pressure component
                pttc = 0.5d0*(-xijs * (ryis + alpk(isol)*yijs)/clrk + yijs * (rxis + alpk(isol)*xijs)/clrk)**2/ &
                  & dABS(xijs * (rxis + alpk(isol)*xijs)/clrk + yijs * (ryis + alpk(isol)*yijs)/clrk) * &
                  & dudr/rijs

                !Update normal pressure from the contribution of fluid-fluid interactions (2)
                virialpress_cylin_ptt(2,ibin) = virialpress_cylin_ptt(2,ibin) + pttc

              End if
            ! End looping over alpha
            End Do
          End If 

          ! --------- P_zz ----------
          ! Loop over solutions
          Do isol = 1,2 

            ! contribute
            if((alpk(isol) .ge. 0.0d0) .and. (alpk(isol) .le. 1.0d0)) Then

              ! Calculate z-intercept with R-surface
              zlk = rzis + alpk(isol)*zijs

              ! Sample within [zlo,zhi]
              if ((zlk .ge. zlo) .and. (zlk .le. zhi)) Then

                ! R_l at alpk
                clrk = dSQRT( (rxis + alpk(isol)*xijs)**2 + (ryis + alpk(isol)*yijs)**2 )
                ! axial pressure component
                ptzc = zijssq/dABS(xijs * (rxis + alpk(isol)*xijs)/clrk + yijs * (ryis + alpk(isol)*yijs)/clrk) * dudr/rijs
                ! Update pzz
                virialpress_cylin_ptz(2,ibin) = virialpress_cylin_ptz(2,ibin) + ptzc

              End If
            End if
          End Do

        !End bin cycle
        ENDDO

      !End loop over jmol sites 
      ENDDO       
    !End IF imol < n_mol_tot  
    ENDIF    
  !End loop over imol sites
  ENDDO  

  If(MOD(iframe,1000) .eq. 0) Then
    Write(*,'(A,I10)') 'Finished sampling of frame #', iframe
    ! unit=6 is screen
    Call FLUSH(6)
  End If

! End loop over frames
End do




! ------------- Average statistics and write statistics to file ------------------!
! --------------- Radial density ---------------!
! Open file
OPEN(2,FILE='r-density.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

! Loop over bins
Do ibin = 1, rden_bins
  ! Loop over molecule types
  Do itype = 1, n_mol_types

    rdenavg(ibin,itype) = rdenavg(ibin,itype)/DBLE(n_frame)
      ! Convert unit from [1/A^3] to [g/ml]
      rdenavg(ibin,itype) = (mol_mass/(Na*1.0d-24))*rdenavg(ibin,itype)

  End Do
End Do

! Write data to file
! loop over molecule types
Do itype = 1, n_mol_types
  ! Write molecule info
  Write(2,'(2A)') 'Molecule name: ', site_name
  Write(2,*) ' R             R-rho [g/ml]           R-rho [1/A^3]'

  ! Loop over bins
  Do ibin = 1, rden_bins

    ! Write to file
    Write(2,'(F8.4,7X,E15.7,8X,E15.7)')  (DBLE(ibin)-0.5d0)*rden_dr, rdenavg(ibin,itype), &
                      & (rdenavg(ibin,itype)/mol_mass)*Na*1.0d-24
    
  End Do
End Do

! Close file
CLOSE(2)

! ------------ Cylindrical pressure tensor ----------------!
Do ibin = 1, rden_bins

  ! Pressure (Kinetic part) in unit of [K/A^3]
  prkin(ibin) = (rdenavg(ibin,1)/mol_mass)*Na*1.0d-24*temp

  ! Get r-distance of ibin in cylindrical system
  posr = (DBLE(ibin)-0.5d0)*rden_dr

  ! Irving-Kirkwood route
  ! Fluid-fluid contribution
  virialpnravg(2,ibin) = virialpress_cylin_pnr(2,ibin)/DBLE(n_frame)
  virialpttavg(2,ibin) = virialpress_cylin_ptt(2,ibin)/DBLE(n_frame)
  virialptzavg(2,ibin) = virialpress_cylin_ptz(2,ibin)/DBLE(n_frame)

  ! Calculate final pressure (configurational part) in unit of [K/A^3]
  virialpnravg(2,ibin) = -1.0d0/(two_Pi*posr*(zhi-zlo))*virialpnravg(2,ibin)
  virialpttavg(2,ibin) = -1.0d0/(two_Pi*posr*(zhi-zlo))*virialpttavg(2,ibin)
  virialptzavg(2,ibin) = -1.0d0/(two_Pi*posr*(zhi-zlo))*virialptzavg(2,ibin)

  ! Calculate totoal pressure tensor in unit of [K/A^3]
  virialpnravg(3,ibin) = virialpnravg(2,ibin)
  virialpttavg(3,ibin) = virialpttavg(2,ibin)
  virialptzavg(3,ibin) = virialptzavg(2,ibin)


End Do

! Open file
OPEN(3,FILE='press_cylinIK.txt', STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

! Write file head
Write(3,*) 'Cylindrical pressure tensor from virial route using Irving-Kirkwood definition'
Write(3,*) 'Unit: pressure in [bar] and length in [Angstrom]'
Write(3,'(A)') &
        & '  R               Pkin          Pnr(ff)             Pnr(tot)          Ptz(ff)&
        &         Ptz(tot)         Ptt(ff)          Ptt(tot)'

! Write data to file
! Loop over bins
Do ibin = 1, rden_bins

  ! Write to file
  Write(3,'(F8.4,F16.4,F16.4,3X,F16.4,3X,F16.4,F16.4,F16.4,F16.4)')  &
    & (DBLE(ibin)-0.5d0)*rden_dr, prkin(ibin)*PCOEFF, &
    & virialpnravg(2,ibin)*PCOEFF, &
    & (virialpnravg(3,ibin)+prkin(ibin))*PCOEFF, &
    & virialptzavg(2,ibin)*PCOEFF, &
    & (virialptzavg(3,ibin)+prkin(ibin))*PCOEFF, &
    & virialpttavg(2,ibin)*PCOEFF, &
    & (virialpttavg(3,ibin)+prkin(ibin))*PCOEFF

End Do

! Close file
CLOSE(3)




! ---------- Deallocate variables to free space -------------!
Deallocate(box)
Deallocate(boxangle)
Deallocate(rx_s)
Deallocate(ry_s)
Deallocate(rz_s)
Deallocate(virialpress_cylin_pnr) 
Deallocate(virialpress_cylin_ptt) 
Deallocate(virialpress_cylin_ptz) 
Deallocate(virialpnravg) 
Deallocate(virialpttavg) 
Deallocate(virialptzavg) 
Deallocate(prkin) 
Deallocate(rdenavg) 
Deallocate(rdenblk) 
  


END PROGRAM virialpress_cylinder
