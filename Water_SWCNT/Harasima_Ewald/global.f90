! ==========================================================
! This module is used to set global variables for 
! cylindrical pressure tensor calculation
! Created on 9/24/2019
! ==========================================================


	  Module global

	  IMPLICIT NONE


	  ! Simulation box info
	  !================================================!
	  ! Box size dimension
	  Double Precision, Dimension(:,:), Allocatable :: box

	  ! Box angle (shape)
	  Double Precision, Dimension(:,:), Allocatable :: boxangle

	  ! Box volume
	  Double Precision, Dimension(:), Allocatable :: volbox


	  ! Molecular parameters
	  ! ===============================================!
	  ! Coordinates
      Double Precision, Dimension(:,:,:), Allocatable :: rx_s, ry_s, rz_s
      Real, Dimension(:), Allocatable :: rxs_temp, rys_temp, rzs_temp
      Double Precision, Dimension(:,:), Allocatable :: rx, ry, rz

      ! Coordinates of the nanotube (assuming the nanotube is fixed in the box) 
      Real, Dimension(:), Allocatable :: rx_snt, ry_snt, rz_snt
      Double Precision, save :: rx_nt, ry_nt, rz_nt



      ! Ewald Method parameters
      ! ===============================================!
	  ! Cut-off radius in real space 
	  Double Precision, save :: rcelect, rcelectsq

	  ! Real space convergence parameter
	  Double Precision, save :: alpha

	  ! Max number of reciprocal vector
	  Integer, save :: kxmax, kymax, kzmax

	  ! Maximum number of reciprocal vector 
	  Integer, save :: ksq_max
	  Integer, save :: maxk 

	  ! Reciprocal prefactor 
	  Double Precision, Dimension(:,:), Allocatable :: k_vec

	  ! Reciprocal component vectors 
	  Complex*16, Dimension(:,:,:), Allocatable :: eikx, eiky, eikz

	  ! Reciprocal structure factor
	  Complex*16, Dimension(:), Allocatable :: skewld, skewlds



	  ! Pressure statistics
	  !================================================!
	  ! 1st rank: 1-fluid-wall; 2-fluid-fluid
	  ! 2nd rank: 1 - real space ; 2-Fourier; 3 - Lennard-Jones
	  ! 3rd rank: bin id
	  !Double Precision, Dimension(:,:,:), Allocatable :: virialpress_cylin_pnr
	  !Double Precision, Dimension(:,:,:), Allocatable :: virialpress_cylin_ptt
	  Double Precision, Dimension(:,:,:), Allocatable :: virialpress_cylin_ptz
	  ! 1st rank: 1-fluid-wall; 2-fluid-fluid; 3-total
	  !Double Precision, Dimension(:,:), Allocatable :: virialpnravg
	  !Double Precision, Dimension(:,:), Allocatable :: virialpttavg
	  Double Precision, Dimension(:,:,:), Allocatable :: virialptzavg
	  Double Precision, Dimension(:), Allocatable :: prkin


	  ! Constants [CODATA 2010 recommendation]
	  !================================================!
	  Double Precision, Parameter :: Pi = 3.141592653589d0
	  Double Precision, Parameter :: sqrtPi = dSQRT(Pi)
	  Double Precision, parameter :: two_Pi = 2.0d0*Pi
	  Double Precision, Parameter :: Kb = 1.3806488d-23   ! [J/K] 
	  Double Precision, Parameter :: Na = 6.02214129d23  ! [1/mol] */
	  Double Precision, Parameter :: h = 6.62606957d-34   ! [J*s] = [kg*m^2/s]
	  Double Precision, Parameter :: R = Kb*Na            ! [J/(mol*K)] */
	  Double Precision, Parameter :: EETOK = 167100.9558738398d0 ! EETOK = 1/(4*pi*epsilon*kb) * e^2 * 10^10
	  Double Precision, Parameter :: PVTOK = 0.7242971565d-02   ! PVTOK = 1.0d5*1.0d-30/kb , [bar]*[A^3] to [K]
	  Double Precision, Parameter :: PCOEFF = Kb*1.0d25        ! [K]/[A]^3 to [bar]
	  
	  


	  End Module global





      



