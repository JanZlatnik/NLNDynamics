
!--------------------------< Radial Schrodinger Equation >--------------------------!
!                                                                                   !    
! Contains: Subroutine to calculate the eigenstates and eigenenergies in            !
!            a radial potential V(x)                                                !                                                                             
!                                                                                   !
! Last revision:    18/03/2024                                                      !
!                                                                                   !
!-----------------------------------------------------------------------------------!

MODULE RadialSE
    USE Parameters
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eig, tridiag_H_matrix, Four_DVR, Four_DVR_fixed, Inverse_F
    
    ! Interface for R->R potential V(x)
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.0d0)) FUNCTION potential_eig_interface(x)
            IMPLICIT NONE
            REAL(KIND = KIND(1.0d0)), INTENT(IN) :: x
        END FUNCTION potential_eig_interface
    END INTERFACE

    CONTAINS 
    
    ! Construction of tridiagonal matrix of the Hamiltonian with a potential V(x)
    SUBROUTINE tridiag_H_matrix(x, potential_eig, m, diag, side_diag)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: m
        REAL(KIND = idk), INTENT(IN) :: x(:)
        PROCEDURE(potential_eig_interface), POINTER, INTENT(IN) :: potential_eig
        REAL(KIND = idk), INTENT(OUT) :: diag(:), side_diag(:)
        INTEGER :: i
        REAL(KIND = idk) :: h

        h = ABS(x(2) - x(1))

        DO i = 1, SIZE(diag)
            diag(i) = hbar**2/m * 1/h**2 + potential_eig(x(i))
        END DO
        DO i = 1, SIZE(side_diag)
            side_diag(i) = -hbar**2/(2*m) * 1/h**2
        END DO

    END SUBROUTINE tridiag_H_matrix

    ! Calculation of eigenenergiers and eigenstates of a tridiagonal matrix
    SUBROUTINE eig_tridiag(diag,side_diag,eigenvalues,eigenvectors)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(INOUT) :: diag(:)
        REAL(KIND = idk), INTENT(INOUT) :: side_diag(:)
        REAL(KIND = idk), INTENT(OUT) :: eigenvalues(:)
        REAL(KIND = idk), INTENT(OUT) :: eigenvectors(:,:)
        INTEGER :: info, lwork, liwork
        REAL(KIND = idk), ALLOCATABLE :: work(:)
        INTEGER, ALLOCATABLE :: iwork(:)
        CHARACTER(LEN=1) :: compz
        REAL(KIND = idk) :: vl, vu
        INTEGER :: il, iu
        

        compz = 'I'
        lwork = max(1, 1 + 4*SIZE(diag) + SIZE(diag)**2)
        liwork = max(1, 3 + 5*SIZE(diag))

        ALLOCATE(work(lwork))
        ALLOCATE(iwork(liwork))

        CALL dstedc(compz, SIZE(diag), diag, side_diag, eigenvectors, SIZE(diag), work, lwork, iwork, liwork, info)

        IF (info /= 0) THEN
            PRINT*, 'Error: DSTEDC failed with INFO = ', info
        END IF

        DEALLOCATE(work)
        DEALLOCATE(iwork)

        eigenvalues = diag
        eigenvectors = TRANSPOSE(eigenvectors)

    END SUBROUTINE eig_tridiag

    ! Given a potential V(x) this subroutine returns the first n eigenvectors and eigenenergies
    SUBROUTINE eig(x, potential_eig, m, n, eigenenergies, eigenvectors)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: m
        INTEGER, INTENT(IN) :: n
        REAL(KIND = idk), INTENT(IN) :: x(:)
        PROCEDURE(potential_eig_interface), POINTER, INTENT(IN) :: potential_eig
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: eigenenergies(:)
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: eigenvectors(:,:)
        REAL(KIND = idk), ALLOCATABLE :: diag(:), side_diag(:), full_eigenenergies(:), full_eigenvectors(:,:)

        ALLOCATE(diag(SIZE(x)))
        ALLOCATE(side_diag(SIZE(x)-1))
        ALLOCATE(full_eigenenergies(SIZE(x)))
        ALLOCATE(full_eigenvectors(SIZE(x),SIZE(x)))
        IF (ALLOCATED(eigenenergies)) DEALLOCATE(eigenenergies)
        ALLOCATE(eigenenergies(n))
        IF (ALLOCATED(eigenvectors)) DEALLOCATE(eigenvectors)
        ALLOCATE(eigenvectors(n,SIZE(x)))

        CALL tridiag_H_matrix(x, potential_eig, m, diag, side_diag)
        CALL eig_tridiag(diag, side_diag, full_eigenenergies, full_eigenvectors)

        eigenenergies = full_eigenenergies(1:n)
        eigenvectors = full_eigenvectors(1:n,:)/SQRT(ABS(x(2) - x(1)))

        DEALLOCATE(diag)
        DEALLOCATE(side_diag)
        DEALLOCATE(full_eigenenergies)
        DEALLOCATE(full_eigenvectors)

    END SUBROUTINE eig
    
    
    ! FourierDVR method for finding eigenstates and eigenenergies in a pontential V(x)
    SUBROUTINE Four_DVR(x0, xN, V, mass, Nmax, n, eigenergies, F_eigenvectors)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: mass, x0, xN
        INTEGER, INTENT(IN) :: Nmax, n
        PROCEDURE(potential_eig_interface), POINTER, INTENT(IN) :: V
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: eigenergies(:),  F_eigenvectors(:,:)
        REAL(KIND = idk), ALLOCATABLE :: H_matrix(:,:), all_eigenergies(:), X_optor(:,:), x_eig(:)
        REAL(KIND = idk) :: L
        INTEGER :: i, j, k
        INTEGER :: info, lwork
        CHARACTER(len=1) :: jobz, uplo
        REAL(KIND = idk), ALLOCATABLE :: work(:)
        
        L = ABS(xN - x0)
        
        IF (n>Nmax) THEN
            PRINT *, "Four_DVR failed: n > Nmax"
        ELSE

            ALLOCATE(X_optor(Nmax,Nmax), x_eig(Nmax))
            
            X_optor = 0.0d0
            
            DO i = 2, Nmax
                DO j = (MOD(i, 2) + 1), i, 2
                    X_optor(i,j) = - 8.0d0 * L * i * j / (pi**2 * (i**2 - j**2)**2)
                    X_optor(j,i) = X_optor(i,j)
                END DO
            END DO
            
            jobz = 'V'  
            uplo = 'U'  
            lwork = -1
            ALLOCATE(work(1), STAT=info)
            CALL DSYEV(jobz, uplo, Nmax, X_optor, Nmax, x_eig, work, lwork, info)
            IF (info /= 0) THEN
                PRINT *, "DSYEV failed with info = ", info
                RETURN
            END IF

            lwork = INT(work(1))
            DEALLOCATE(work)
            ALLOCATE(work(lwork)) 
            CALL DSYEV(jobz, uplo, Nmax, X_optor, Nmax, x_eig, work, lwork, info)
            IF (info /= 0) THEN
                PRINT *, "DSYEV failed with info = ", info
                RETURN
            END IF
            DEALLOCATE(work)
            
            x_eig = x_eig + (x0+xN)/2
            
            ALLOCATE(H_matrix(Nmax,Nmax), all_eigenergies(Nmax))
            
            H_matrix = 0.0d0
            DO i = 1, Nmax
                DO j = 1, i
                    DO k = 1, Nmax
                        H_matrix(i,j) = H_matrix(i,j) + X_optor(i,k) * V(x_eig(k)) * X_optor(j,k)
                    END DO
                    H_matrix(j,i) = H_matrix(i,j)
                END DO
                H_matrix(i,i) = H_matrix(i,i) + i**2 * pi**2 * hbar**2 / (2 * mass * L**2)
            END DO
            
            DEALLOCATE(X_optor, x_eig)
            IF (ALLOCATED(F_eigenvectors)) DEALLOCATE(F_eigenvectors)
            IF (ALLOCATED(eigenergies)) DEALLOCATE(eigenergies)
            
            jobz = 'V'  
            uplo = 'U'  
            lwork = -1
            ALLOCATE(work(1), STAT=info)
            CALL DSYEV(jobz, uplo, Nmax, H_matrix, Nmax, all_eigenergies, work, lwork, info)
            IF (info /= 0) THEN
                PRINT *, "DSYEV failed with info = ", info
                RETURN
            END IF

            lwork = INT(work(1))
            DEALLOCATE(work)
            ALLOCATE(work(lwork), STAT=info) 
            CALL DSYEV(jobz, uplo, Nmax, H_matrix, Nmax, all_eigenergies, work, lwork, info)
            DEALLOCATE(work)
            IF (info == 0) THEN
                ALLOCATE(eigenergies(n),F_eigenvectors(n,Nmax))
                eigenergies = all_eigenergies(1:n)
                F_eigenvectors = TRANSPOSE(H_matrix(:,1:n))
            ELSE
                PRINT *, "DSYEV failed with info = ", info
            END IF
            
            DEALLOCATE(all_eigenergies,H_matrix)
            

        END IF

    END SUBROUTINE Four_DVR
    
    
    ! FourierDVR method with fixed potential values on DVR grid
    SUBROUTINE Four_DVR_fixed(x_values, V_values, mass, Nmax, n, eigenergies, F_eigenvectors)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: mass, x_values(:) 
        INTEGER, INTENT(IN) :: Nmax, n
        REAL(KIND = idk), INTENT(IN) :: V_values(Nmax)
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: eigenergies(:),  F_eigenvectors(:,:)
        REAL(KIND = idk), ALLOCATABLE :: H_matrix(:,:), all_eigenergies(:), X_optor(:,:), x_eig(:), V_eig(:)
        REAL(KIND = idk) :: L, x0, xN
        INTEGER :: i, j, k
        INTEGER :: info, lwork
        CHARACTER(len=1) :: jobz, uplo
        REAL(KIND = idk), ALLOCATABLE :: work(:)
        
        x0 = x_values(1)
        xN = x_values(size(x_values))
        L = ABS(xN - x0)
        
        IF (n>Nmax) THEN
            PRINT *, "Four_DVR failed: n > Nmax"
            RETURN
        END IF

            ALLOCATE(X_optor(Nmax,Nmax), x_eig(Nmax))
            
            X_optor = 0.0d0
            
            DO i = 2, Nmax
                DO j = (MOD(i, 2) + 1), i, 2
                    X_optor(i,j) = - 8.0d0 * L * i * j / (pi**2 * (i**2 - j**2)**2)
                    X_optor(j,i) = X_optor(i,j)
                END DO
            END DO
            
            jobz = 'V'  
            uplo = 'U'  
            lwork = -1
            ALLOCATE(work(1), STAT=info)
            CALL DSYEV(jobz, uplo, Nmax, X_optor, Nmax, x_eig, work, lwork, info)
            IF (info /= 0) THEN
                PRINT *, "DSYEV failed with info = ", info
                RETURN
            END IF

            lwork = INT(work(1))
            DEALLOCATE(work)
            ALLOCATE(work(lwork)) 
            CALL DSYEV(jobz, uplo, Nmax, X_optor, Nmax, x_eig, work, lwork, info)
            IF (info /= 0) THEN
                PRINT *, "DSYEV failed X with info = ", info
                RETURN
            END IF
            DEALLOCATE(work)
            
            x_eig = x_eig + (x0+xN)/2
            
            ALLOCATE(H_matrix(Nmax,Nmax), all_eigenergies(Nmax))
            ALLOCATE(V_eig(Nmax))
            CALL interpolate_V(x_values, V_values, size(x_values), x_eig, Nmax, V_eig)
            
            H_matrix = 0.0d0
    
        DO i = 1, Nmax
            DO j = 1, i
                DO k = 1, Nmax
                    H_matrix(i,j) = H_matrix(i,j) + X_optor(i,k) * V_eig(k) * X_optor(j,k)
                END DO
                H_matrix(j,i) = H_matrix(i,j)
            END DO
            H_matrix(i,i) = H_matrix(i,i) + i**2 * pi**2 * hbar**2 / (2 * mass * L**2)
        END DO
    
        DEALLOCATE(X_optor, x_eig, V_eig)
        IF (ALLOCATED(F_eigenvectors)) DEALLOCATE(F_eigenvectors)
        IF (ALLOCATED(eigenergies)) DEALLOCATE(eigenergies)
            
        jobz = 'V'  
        uplo = 'U'  
        lwork = -1
        ALLOCATE(work(1), STAT=info)
        CALL DSYEV(jobz, uplo, Nmax, H_matrix, Nmax, all_eigenergies, work, lwork, info)
        IF (info /= 0) THEN
            PRINT *, "DSYEV failed with info = ", info
            RETURN
        END IF

        lwork = INT(work(1))
        DEALLOCATE(work)
        ALLOCATE(work(lwork), STAT=info) 
        CALL DSYEV(jobz, uplo, Nmax, H_matrix, Nmax, all_eigenergies, work, lwork, info)
        DEALLOCATE(work)
        IF (info == 0) THEN
            ALLOCATE(eigenergies(n),F_eigenvectors(n,Nmax))
            eigenergies = all_eigenergies(1:n)
            F_eigenvectors = TRANSPOSE(H_matrix(:,1:n))
        ELSE
            PRINT *, "DSYEV V failed with info = ", info
        END IF
            
        DEALLOCATE(all_eigenergies,H_matrix)
    
    END SUBROUTINE Four_DVR_fixed
    
    ! Computes wavefunctions on a grid x(:) given its sine coefficients
    SUBROUTINE Inverse_F(x_mesh, F_coefficient, wavefunctions)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: x_mesh(:)
        REAL(KIND = idk), ALLOCATABLE, INTENT(IN) :: F_coefficient(:,:)
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: wavefunctions(:,:)
        REAL(KIND = idk) :: L, x0, xN
        INTEGER :: i, j, k, Nmax, n, n_mesh
        
        n = SIZE(F_coefficient,1)
        n_mesh = SIZE(x_mesh)
        Nmax = SIZE(F_coefficient,2)
        x0 = x_mesh(1)
        xN = x_mesh(n_mesh)
        L = ABS(xN - x0)
        
        ALLOCATE(wavefunctions(n,n_mesh))
        wavefunctions = 0.0d0
        
        DO k = 1, n
            DO i = 1, Nmax
                DO j = 1, n_mesh
                    wavefunctions(k,j) = wavefunctions(k,j) + F_coefficient(k,i) * SQRT(2.0d0 / L) * SIN(i * pi * (x_mesh(j) - x0) / L)
                END DO
            END DO
        END DO
        
    END SUBROUTINE Inverse_F
    
    
    subroutine interpolate_V(x_values, V_values, n_x, x_eig, n_vib, V_interp)
        implicit none
        integer, intent(in) :: n_x, n_vib
        real(kind = idk), intent(in) :: x_values(n_x), V_values(n_x), x_eig(n_vib)
        real(kind = idk), intent(out) :: V_interp(n_vib)
    
        integer :: i, j
        real*8 :: dx, x_left, x_right, V_left, V_right, weight
    
        dx = x_values(2) - x_values(1)
    
        do i = 1, n_vib
            j = int((x_eig(i) - x_values(1)) / dx) + 1
            if (j < 1) j = 1
            if (j > n_x - 1) j = n_x - 1
    
            x_left = x_values(j)
            x_right = x_values(j+1)
            V_left = V_values(j)
            V_right = V_values(j+1)
    
            weight = (x_eig(i) - x_left) / (x_right - x_left)
            V_interp(i) = (1.0d0 - weight) * V_left + weight * V_right
        end do
    
    end subroutine interpolate_V



END MODULE RadialSE
