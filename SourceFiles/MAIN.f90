!-------------------------------------< MAIN >-------------------------------------!
!                                                                                  !    
! Contains: Main program for model and result data for NLVD                        !
!                                                                                  !
! Last revision:    03/03/2026                                                    !                  
!                                                                                  !
!----------------------------------------------------------------------------------!



PROGRAM MAIN
    USE Parameters
    USE RadialSE
    USE VibrationalExcitation
    USE RCModel2
    USE OMP_LIB

    IMPLICIT NONE
    REAL(KIND = idk) :: dR, dE, dME, dMR
    REAL(KIND = idk), ALLOCATABLE :: Rgrid(:), RModelgrid(:), Egrid(:), EModelgrid(:)
    COMPLEX(KIND = idk), ALLOCATABLE :: Flvlshiftmodel(:)
    INTEGER :: i, j, status
    INTEGER :: iounit

    REAL(KIND = idk), ALLOCATABLE :: eigE_V(:), eigF_V(:,:), FeigF_V(:,:)
    REAL(KIND = idk), ALLOCATABLE :: cs_VE(:,:), cs_VE_tot(:), cs_VE_bg(:,:), cs_VE_tot_bg(:), cs_DR(:), cs_DR_wave(:), cs_DR_opt(:)



    ! Interface for R->R function
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.0d0)) FUNCTION real_function_interface(x)
            IMPLICIT NONE
            REAL(KIND = KIND(1.0d0)), INTENT(IN) :: x
        END FUNCTION real_function_interface
    END INTERFACE 

    ! Interface for R->C function
    ABSTRACT INTERFACE
        COMPLEX(KIND = KIND(1.0d0)) FUNCTION complex_function_interface(x)
            IMPLICIT NONE
            REAL(KIND = KIND(1.0d0)), INTENT(IN) :: x
        END FUNCTION complex_function_interface
    END INTERFACE


    PROCEDURE(real_function_interface), POINTER :: V0_ptr, Vd_ptr, g_ptr, f_ptr, VLCP_ptr
    PROCEDURE(complex_function_interface), POINTER :: Flvlshift_ptr

    V0_ptr => V0
    Vd_ptr => Vd
    VLCP_ptr => VLCP
    g_ptr => g
    f_ptr => f
    Flvlshift_ptr => F_lvlshift



    print*, '=============================================================='
    print*, '            NONLOCAL NUCLEAR DYNAMICS COMPUTATION             '
    print*, '=============================================================='

    ! parameters output
    CALL parameters_o()

    ! Create directories
    INQUIRE(DIRECTORY="DATA", EXIST=status)
    IF (.NOT. status) THEN
        CALL SYSTEM("mkdir -p DATA")
        CALL CONSOLE('Directory "DATA" created successfully.')
    ELSE IF (status == -1) THEN
        CALL CONSOLE('Directory "DATA" already exists.')
    END IF

    ! Create directories
    INQUIRE(DIRECTORY="MODEL", EXIST=status)
    IF (.NOT. status) THEN
        CALL SYSTEM("mkdir -p MODEL")
        CALL CONSOLE('Directory "MODEL" created successfully.')
    ELSE IF (status == -1) THEN
        CALL CONSOLE('Directory "MODEL" already exists.')
    END IF

    ! model initialization
    CALL initiate_model()
    CALL CONSOLE('Model successfully initiated.')

    ! Computational coordinate grid
    CALL make_mesh(Rmin, Rmax, Rmp, Rgrid, dR)
    CALL CONSOLE('Coordinate computational grid created successfully.')

    ! Energy computational grid
    CALL make_mesh(Emin, Emax, ep, Egrid, dE)
    CALL CONSOLE('Energy computational grid created successfully.')


    !=================================================================================
    !                                Model Output
    !=================================================================================
    IF (model_output) THEN

        ! Computational coordinate grid
        CALL make_mesh(RModelmin, RModelmax, RModelnv , RModelgrid, dMR)
        CALL CONSOLE('Coordinate computational grid created successfully.')

        ! Energy model grid
        CALL make_mesh(EModelmin, EModelmax, eModelp, EModelgrid, dME)
        CALL CONSOLE('Energy grid for model testing created successfully.')

        ! Printing V0(x) to file
        OPEN(NEWUNIT=iounit, FILE='MODEL/V0.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A20)') '#', 'R [a0]','V0(R) [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(Rgrid)
            WRITE(iounit, '(1E20.12)', advance='no') Rgrid(i)
            WRITE(iounit, '(1E20.12)', advance='no') V0_ptr(Rgrid(i)) * phys_h0
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/V0.txt')

        ! Printing Vd(x) to file
        OPEN(NEWUNIT=iounit, FILE='MODEL/Vd.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A20)') '#', 'R [a0]','Vd(R) [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(Rgrid)
            WRITE(iounit, '(1E20.12)', advance='no') Rgrid(i)
            WRITE(iounit, '(1E20.12)', advance='no') Vd_ptr(Rgrid(i)) * phys_h0
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/Vd.txt')

        ! Printing VLPC(x) to file
        OPEN(NEWUNIT=iounit, FILE='MODEL/VLCP.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A20, A20)') '#', 'R [a0]','VLCP(R) [eV]','GammaLCP(R) [eV]'
        WRITE(iounit,*)
        DO i = 1, SIZE(Rgrid)
            WRITE(iounit, '(1E20.12)', advance='no') Rgrid(i)
            WRITE(iounit, '(1E20.12)', advance='no') VLCP(Rgrid(i)) * phys_h0
            WRITE(iounit, '(1E20.12)', advance='no') GammaLCP(Rgrid(i)) * phys_h0
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/VLCP.txt')

        ! Printing g(x) to file
        OPEN(NEWUNIT=iounit, FILE='MODEL/g.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A20)') '#', 'R [a0]','g(R) [a.u.]'
        WRITE(iounit,*)
        DO i = 1, SIZE(Rgrid)
            WRITE(iounit, '(1E20.12)', advance='no') Rgrid(i)
            WRITE(iounit, '(1E20.12)', advance='no') g(Rgrid(i))
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/g.txt')


        CALL CONSOLE('Computing complex level-shift operator for output...')
        ALLOCATE(Flvlshiftmodel(SIZE(EModelgrid)))
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
        DO i = 1, SIZE(EModelgrid)
            Flvlshiftmodel(i) = Flvlshift_ptr(EModelgrid(i))
        END DO
        !$OMP END PARALLEL DO
        CALL CONSOLE('Computation of complex level-shift operator finished successfully.')

        ! Printing Delta(E,R) to file
        OPEN(NEWUNIT=iounit, FILE='MODEL/delta.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A20)') '#', 'E [eV]', 'Delta(R,E) [eV] ...'
        WRITE(iounit,*)
        DO i = 1, SIZE(EModelgrid)
            WRITE(iounit, '(1E20.12)', advance='no') EModelgrid(i) * phys_h0
            DO j = 1, SIZE(RModelgrid)
                WRITE(iounit, '(1E20.12)', advance='no') REAL(Flvlshiftmodel(i)) * phys_h0 * (g_ptr(RModelgrid(j))) ** 2
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/delta.txt')

        ! Printing Gamma(E,R) to file
        OPEN(NEWUNIT=iounit, FILE='MODEL/gamma.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A)') '#', 'E [eV]', 'Gamma(R,E) [eV] ...'
        WRITE(iounit,*)
        DO i = 1, SIZE(EModelgrid)
            WRITE(iounit, '(1E20.12)', advance='no') EModelgrid(i) * phys_h0
            DO j = 1, SIZE(RModelgrid)
                WRITE(iounit, '(1E20.12)', advance='no') -2.0d0 * AIMAG(Flvlshiftmodel(i)) * (g_ptr(RModelgrid(j)))**2 * phys_h0
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/gamma.txt')

        DEALLOCATE(Flvlshiftmodel)



        ! Calculation of eigenstates in V0(x)
        CALL CONSOLE('Computing eigenstates in V0 potential...')
        IF (ALLOCATED(eigE_V)) DEALLOCATE(eigE_V)
        IF (ALLOCATED(eigF_V)) DEALLOCATE(eigF_V)
        IF (ALLOCATED(FeigF_V)) DEALLOCATE(FeigF_V)
        CALL Four_DVR(Rmin, Rmax, V0_ptr, mass, n_fdvr, n_fdvr, eigE_V, FeigF_V)
        CALL Inverse_F(Rgrid, FeigF_V, eigF_V)
        IF (ALLOCATED(FeigF_V)) DEALLOCATE(FeigF_V)
        CALL CONSOLE('Eigenstates in V0 potential computed successfully.')
        
        OPEN(NEWUNIT=iounit, FILE='MODEL/V0_eigE.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A20)') '#', 'n', 'E_n [eV]'
        DO i = 1, n_FDVR
            WRITE(iounit, '(I0)', advance='no') i
            WRITE(iounit, '(1E20.12)', advance='no') eigE_V(i)*phys_h0
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/V0_eigE.txt')
        
        OPEN(NEWUNIT=iounit, FILE='MODEL/V0_eigF.txt', STATUS='replace', ACTION='write')    
        WRITE(iounit, '(A3, A17, A)') '#', 'R [a0]', 'nu_i [u.a.]'

        DO j = 1, SIZE(Rgrid)
            WRITE(iounit, '(1E20.12)', advance='no') Rgrid(j)
            DO i = 1, n_print+1
                WRITE(iounit, '(1E20.12)', advance='no') eigF_V(i,j)
            END DO
            WRITE(iounit,*) 
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/V0_eigF.txt')
        
        IF (ALLOCATED(eigE_V)) DEALLOCATE(eigE_V)
        IF (ALLOCATED(eigF_V)) DEALLOCATE(eigF_V)


        ! Calculation of eigenstates in Vd(x)
        CALL CONSOLE('Computing eigenstates in Vd potential...')
        IF (ALLOCATED(eigE_V)) DEALLOCATE(eigE_V)
        IF (ALLOCATED(eigF_V)) DEALLOCATE(eigF_V)
        IF (ALLOCATED(FeigF_V)) DEALLOCATE(FeigF_V)
        CALL Four_DVR(Rmin, Rmax, Vd_ptr, mass, n_fdvr, n_fdvr, eigE_V, FeigF_V)
        CALL Inverse_F(Rgrid, FeigF_V, eigF_V)
        IF (ALLOCATED(FeigF_V)) DEALLOCATE(FeigF_V)
        CALL CONSOLE('Eigenstates in Vd potential computed successfully.')
        
        OPEN(NEWUNIT=iounit, FILE='MODEL/Vd_eigE.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A20)') '#', 'n', 'E_n [eV]'
        DO i = 1, n_FDVR
            WRITE(iounit, '(I0)', advance='no') i
            WRITE(iounit, '(1E20.12)', advance='no') eigE_V(i)*phys_h0
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/Vd_eigE.txt')
        
        OPEN(NEWUNIT=iounit, FILE='MODEL/Vd_eigF.txt', STATUS='replace', ACTION='write')    
        WRITE(iounit, '(A3, A17, A)') '#', 'R [a0]', 'nu_i [u.a.]'

        DO j = 1, SIZE(Rgrid)
            WRITE(iounit, '(1E20.12)', advance='no') Rgrid(j)
            DO i = 1, n_print+1
                WRITE(iounit, '(1E20.12)', advance='no') eigF_V(i,j)
            END DO
            WRITE(iounit,*) 
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/Vd_eigF.txt')
        
        IF (ALLOCATED(eigE_V)) DEALLOCATE(eigE_V)
        IF (ALLOCATED(eigF_V)) DEALLOCATE(eigF_V)

        
        ! Calculation of eigenstates in VLCP(x)
        CALL CONSOLE('Computing eigenstates in VLCP potential...')
        IF (ALLOCATED(eigE_V)) DEALLOCATE(eigE_V)
        IF (ALLOCATED(eigF_V)) DEALLOCATE(eigF_V)
        IF (ALLOCATED(FeigF_V)) DEALLOCATE(FeigF_V)
        CALL Four_DVR(Rmin, Rmax, VLCP_ptr, mass, n_fdvr, n_fdvr, eigE_V, FeigF_V)
        CALL Inverse_F(Rgrid, FeigF_V, eigF_V)
        IF (ALLOCATED(FeigF_V)) DEALLOCATE(FeigF_V)
        CALL CONSOLE('Eigenstates in VLCP potential computed successfully.')
        
        OPEN(NEWUNIT=iounit, FILE='MODEL/VLCP_eigE.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A1, A2, A20)') '#', 'n', 'E_n [eV]'
        DO i = 1, n_FDVR
            WRITE(iounit, '(I0)', advance='no') i
            WRITE(iounit, '(1E20.12)', advance='no') eigE_V(i)*phys_h0
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/VLCP_eigE.txt')
        
        OPEN(NEWUNIT=iounit, FILE='MODEL/VLCP_eigF.txt', STATUS='replace', ACTION='write')    
        WRITE(iounit, '(A3, A17, A)') '#', 'R [a0]', 'nu_i [u.a.]'

        DO j = 1, SIZE(Rgrid)
            WRITE(iounit, '(1E20.12)', advance='no') Rgrid(j)
            DO i = 1, n_print+1
                WRITE(iounit, '(1E20.12)', advance='no') eigF_V(i,j)
            END DO
            WRITE(iounit,*) 
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/VLCP_eigF.txt')
        
        IF (ALLOCATED(eigE_V)) DEALLOCATE(eigE_V)
        IF (ALLOCATED(eigF_V)) DEALLOCATE(eigF_V)


    END IF


    !=================================================================================
    !                      Vibrational Excitation & DA / DR
    !=================================================================================
    IF (VE_computation) THEN
        CALL CONSOLE('Running computation of VE & DA/DR...')
        CALL vib_exc(Rgrid, Egrid, mass, V0_ptr, n_fdvr, Vd_ptr, 0, n_eig, g_ptr, f_ptr, Flvlshift_ptr, cs_VE, cs_VE_bg, cs_VE_tot, cs_VE_tot_bg, cs_DR, cs_DR_wave, cs_DR_opt)
        CALL CONSOLE('Computation of VE & DA/DR successfully finished.')

        OPEN(NEWUNIT=iounit, FILE='DATA/cs_VE.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A)') '#', 'e [eV]', 'cs_VE [au] ...'
        DO j = 1, SIZE(Egrid)
            WRITE(iounit, '(1E20.12)', advance='no') Egrid(j) * phys_h0
            DO i = 1, n_print+1
                WRITE(iounit, '(1E20.12)', advance='no') cs_VE(j,i)
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/cs_VE.txt')

        
        OPEN(NEWUNIT=iounit, FILE='DATA/cs_VE_bg.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A)') '#', 'e [eV]', 'cs_VE_bg [au] ...'
        DO j = 1, SIZE(Egrid)
            WRITE(iounit, '(1E20.12)', advance='no') Egrid(j) * phys_h0
            DO i = 1, n_print+1
                WRITE(iounit, '(1E20.12)', advance='no') cs_VE_bg(j,i)
            END DO
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/cs_VE_bg.txt')

        
        OPEN(NEWUNIT=iounit, FILE='DATA/cs_VE_tot.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A20, A20)') '#', 'e [eV]', 'cs_VE_tot [au]', 'cs_VE_tot_bg [au]'
        DO j = 1, SIZE(Egrid)
            WRITE(iounit, '(1E20.12)', advance='no') Egrid(j) * phys_h0
            WRITE(iounit, '(1E20.12)', advance='no') cs_VE_tot(j)
            WRITE(iounit, '(1E20.12)', advance='no') cs_VE_tot_bg(j)
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/cs_VE_tot.txt')

        
        OPEN(NEWUNIT=iounit, FILE='DATA/cs_DR.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A20, A20, A20)') '#', 'e [eV]', 'cs_DR [au]', 'cs_DR_wave [au]', 'cs_DR_opt [au]'
        DO j = 1, SIZE(Egrid)
            WRITE(iounit, '(1E20.12)', advance='no') Egrid(j) * phys_h0
            WRITE(iounit, '(1E20.12)', advance='no') cs_DR(j)
            WRITE(iounit, '(1E20.12)', advance='no') cs_DR_wave(j)
            WRITE(iounit, '(1E20.12)', advance='no') cs_DR_opt(j)
            WRITE(iounit,*)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to DATA/cs_DR.txt')

        IF(ALLOCATED(cs_VE)) DEALLOCATE(cs_VE)
        IF(ALLOCATED(cs_VE_tot)) DEALLOCATE(cs_VE_tot)
        IF(ALLOCATED(cs_VE_bg)) DEALLOCATE(cs_VE_bg)
        IF(ALLOCATED(cs_VE_tot_bg)) DEALLOCATE(cs_VE_tot_bg)
        IF(ALLOCATED(cs_DR)) DEALLOCATE(cs_DR)
        IF(ALLOCATED(cs_DR_wave)) DEALLOCATE(cs_DR_wave)
        IF(ALLOCATED(cs_DR_opt)) DEALLOCATE(cs_DR_opt)

    END IF

    IF(ALLOCATED(Rgrid)) DEALLOCATE(Rgrid)
    IF(ALLOCATED(Egrid)) DEALLOCATE(Egrid)
    IF(ALLOCATED(RModelgrid)) DEALLOCATE(RModelgrid)
    IF(ALLOCATED(EModelgrid)) DEALLOCATE(EModelgrid)


END PROGRAM MAIN