
!--------------------------< Parameters for NLNDynamics >--------------------------!
!                                                                                  !    
! Contains: Basic mathematical and physical constants, parameters, mesh creation   !                                                                             
!                                                                                  !
! Last revision:    02/03/2026                                                     !
!                                                                                  !
!----------------------------------------------------------------------------------!
    
MODULE Parameters
    IMPLICIT NONE
    INTEGER, PARAMETER :: idk = 8
    INTEGER, PARAMETER :: qdk = 16
    
    ! Mathematical constants:
    REAL(KIND = idk), PARAMETER :: pi = 3.141592653589793238462643383279502d0
    REAL(KIND = idk), PARAMETER :: euler = 2.718281828459045235d0
    
    ! Physical constants:
    REAL(KIND = idk), PARAMETER :: &
        phys_a0 = 0.52917706d0, &   ! Bohr radius [angstroms]
        phys_h0 = 27.2113961d0, &   ! 1 Hartree in [eV]
        phys_k = 8.617385d-5, &     ! Boltzman konstant in [eV/K]
        phys_me = 9.1093897d-31,&   ! electron mass [kg]
        phys_u = 1.6605402d-27, &   ! (unified) atomic mass unit [kg]
        hbar = 1.0d0, &             ! Planck constant [au]
        cmtoht = 4.556d-6           ! cm^-1 to hartree
    
    
    ! Computational settings
    LOGICAL, PARAMETER                  :: model_output = .TRUE.
    LOGICAL, PARAMETER                  :: VE_computation = .TRUE.

    
    
    ! Mesh settings
    REAL(KIND = idk), PARAMETER :: Rmin = 1.0d-10         ! in [au]
    REAL(KIND = idk), PARAMETER :: Rmax = 10.0d0          ! in [au]
    INTEGER, PARAMETER          :: Rmp = 3000
    
    ! Energy mesh settings
    REAL(KIND = idk), PARAMETER :: Emin = 0.1d-4/phys_h0       ! [eV] to [au]
    REAL(KIND = idk), PARAMETER :: Emax = 1.49996d0/phys_h0    ! [eV] to [au]
    INTEGER, PARAMETER          :: ep = 30000

    ! Energy model grid settings
    REAL(KIND = idk), PARAMETER :: EModelmin = -15.0d0/phys_h0       ! [eV] to [au]
    REAL(KIND = idk), PARAMETER :: EModelmax = 15.0d0/phys_h0        ! [eV] to [au]
    INTEGER, PARAMETER          :: eModelp = 5000

    ! R model grid settings
    REAL(KIND = idk), PARAMETER :: RModelmin = 0.7d0        ! in [au]
    REAL(KIND = idk), PARAMETER :: RModelmax = 10.0d0       ! in [au]
    INTEGER, PARAMETER          :: RModelnv = 10

    ! Other parameters
    INTEGER, PARAMETER          :: n_eig = 80   ! number of vib. eigenstates used in level-shift operator expansion
    INTEGER, PARAMETER          :: n_fdvr = 200 ! size of FDVR basis to compute eigenstates
    INTEGER, PARAMETER          :: n_print = 6  ! number of printed final states
    
    CONTAINS
    
    SUBROUTINE make_mesh(xmin, xmax, mp, x, dx)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(OUT) :: dx
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: x(:)
        INTEGER, INTENT(IN) :: mp
        REAL(KIND = idk), INTENT(IN) :: xmin, xmax 
        INTEGER :: i
        
        IF (ALLOCATED(x)) DEALLOCATE(x)
        ALLOCATE(x(mp))
        dx = (xmax - xmin) / (mp-1.0d0)
        DO i = 1,mp
            x(i) = xmin + (i-1)*dx
        END DO
        
    END SUBROUTINE


    SUBROUTINE make_mesh_ndyn(xmin, xmax, xcore1, xcore2, mp, x)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: xmin, xmax, xcore1, xcore2
        INTEGER, INTENT(IN) :: mp
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: x(:)
        
        INTEGER :: mp1, mp2, mp3, i, steps2
        REAL(KIND = idk) :: dx1, dx2, dx3
        
        steps2 = MAX(5, NINT((0.60d0 * REAL(mp, idk)) / 5.0d0) * 5)
        mp2 = steps2 + 1
        
        mp1 = MAX(2, NINT(0.15d0 * REAL(mp, idk)))
        
        mp3 = mp - mp1 - mp2 + 2
        
        IF (ALLOCATED(x)) DEALLOCATE(x)
        ALLOCATE(x(mp))
        
        dx1 = (xcore1 - xmin) / REAL(mp1 - 1, idk)
        dx2 = (xcore2 - xcore1) / REAL(mp2 - 1, idk)
        dx3 = (xmax - xcore2) / REAL(mp3 - 1, idk)
        
        DO i = 1, mp1
            x(i) = xmin + REAL(i - 1, idk) * dx1
        END DO
        
        DO i = 2, mp2
            x(mp1 + i - 1) = xcore1 + REAL(i - 1, idk) * dx2
        END DO

        DO i = 2, mp3
            x(mp1 + mp2 - 2 + i) = xcore2 + REAL(i - 1, idk) * dx3
        END DO
        
    END SUBROUTINE make_mesh_ndyn
    
    
    SUBROUTINE CONSOLE(message)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: message
        CHARACTER(LEN=10) :: timestr
        INTEGER :: h, m, s
        
        CALL DATE_AND_TIME(TIME=timestr)

        READ(timestr(1:2),*) h
        READ(timestr(3:4),*) m
        READ(timestr(5:6),*) s

        WRITE(*,'(A,I2.2,A,I2.2,A,I2.2,A,1X,A)') '[', h, ':', m, ':', s, ']:', TRIM(message)
    END SUBROUTINE CONSOLE
    
    FUNCTION STR(k) RESULT(str_out)
        INTEGER, INTENT(IN) :: k
        CHARACTER(LEN=20)   :: str_out
        WRITE(str_out, '(I0)') k 
        str_out = ADJUSTL(str_out)    
    END FUNCTION STR


    FUNCTION cutoff_energy(l) RESULT(Ecut)
        INTEGER, INTENT(IN) :: l
        REAL(KIND = idk) :: Ecut
        IF (l == 0) THEN
            Ecut = -1.0d0
        ELSE
            Ecut = -1.0d0 / (2.0d0 * REAL(l, KIND=idk)**2)
        END IF
    END FUNCTION cutoff_energy
    
    
    
    
    SUBROUTINE parameters_write(unit)
        INTEGER, INTENT(IN)  :: unit
        
        WRITE(unit,*) "=============================================================="
        WRITE(unit,*) "                  Compuational Settings                       "
        WRITE(unit,*) "=============================================================="
        
        WRITE(unit,*) ""
        WRITE(unit,'(A35,1X,L1)') "model_output:", model_output
        WRITE(unit,'(A35,1X,L1)') "VE_computation:", VE_computation
        WRITE(unit,*) ""

        
        WRITE(unit,*) "=============================================================="
        WRITE(unit,*) "                 Computation Parameters                       "
        WRITE(unit,*) "=============================================================="
        
        WRITE(unit,*) ""
        WRITE(unit,*) "[x-Grid]"
        WRITE(unit,'(A25,1X,E0.6)') "Rmin:", Rmin
        WRITE(unit,'(A25,1X,F0.6)') "Rmax:", Rmax
        WRITE(unit,'(A25,1X,I0)')   "Rmp:", Rmp 
        WRITE(unit,*) ""
        
        WRITE(unit,*) ""
        WRITE(unit,*) "[E-Grid]"
        WRITE(unit,'(A25,1X,F0.6)') "Emin:", Emin * phys_h0
        WRITE(unit,'(A25,1X,F0.6)') "Emax:", Emax * phys_h0
        WRITE(unit,'(A25,1X,I0)')   "ep:", ep 
        WRITE(unit,*) ""     

        WRITE(unit,*) ""
        WRITE(unit,*) "[E-Grid]"
        WRITE(unit,'(A25,1X,F0.6)') "EModelmin:", EModelmin * phys_h0
        WRITE(unit,'(A25,1X,F0.6)') "EModelmax:", EModelmax * phys_h0
        WRITE(unit,'(A25,1X,I0)')   "eModelp:", eModelp 
        WRITE(unit,*) ""     

        WRITE(unit,*) ""
        WRITE(unit,*) "[Other Parameters]"
        WRITE(unit,'(A25,1X,I0)')   "n_eig:", n_eig
        WRITE(unit,'(A25,1X,I0)')   "n_fdvr:", n_fdvr 
        WRITE(unit,'(A25,1X,I0)')   "n_print:", n_print 
        WRITE(unit,*) ""   
        
        
        WRITE(unit,*) "=============================================================="
    
    END SUBROUTINE parameters_write
    
    
    SUBROUTINE parameters_o()
        INTEGER :: param_u
        INTEGER :: console_u = 6
        
        CALL parameters_write(console_u)
        
        OPEN(newunit = param_u, file = 'parameters.txt', status = 'replace', form = 'formatted')
        CALL parameters_write(param_u)
        CALL CONSOLE('Parameters have been succesfully written to parameters.txt.')
    
    
    END SUBROUTINE



    
    
END MODULE Parameters