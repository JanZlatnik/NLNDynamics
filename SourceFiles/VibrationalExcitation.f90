!----------------------------< Vibrational Excitation >----------------------------!
!                                                                                  !    
! Contains: Calculation of vibrational excitation amplitudes and cross sections    !                                                   
!                                                                                  !
! Last revision:    22/03/2024                                                     !
!                                                                                  !
!----------------------------------------------------------------------------------!

MODULE VibrationalExcitation
    USE Parameters
    USE Math
    USE RadialSE
    USE Green
    USE OMP_LIB
    USE RCModel2, ONLY : DS_R0, compute_Tbg
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: vib_exc
    
    ! Abstract interface for real->real functions
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.d0)) FUNCTION real_ve_interface(x)
            IMPLICIT NONE
            REAL(KIND = KIND(1.d0)), INTENT(IN) :: x
        END FUNCTION real_ve_interface
    END INTERFACE

    ! Abstract interface for real->complex functions
    ABSTRACT INTERFACE
        COMPLEX(KIND = KIND(1.d0)) FUNCTION complex_ve_interface(x)
            IMPLICIT NONE
            REAL(KIND = KIND(1.d0)), INTENT(IN) :: x
        END FUNCTION complex_ve_interface
    END INTERFACE


    
    CONTAINS
    
    ! Subroutine that calculates cross sections of vibrational exitation for multiple given (as a field) energy of the incoming electron using F_eig
    SUBROUTINE vib_exc(x_mesh,e_i,mass,ground_pot,N_Four,exc_pot,in_state,n_states,g_Vde,f_Vde,F_trans,cs_VE,cs_VE_bg,cs_VE_tot,cs_VE_tot_bg,cs_DA,cs_DA_infty,cs_DA_inel)
        ! x_mesh        -   interval (xmin,xmax) as a mesh in [au]
        ! e_i           -   vector of energies of the incoming electron in [au]
        ! mass          -   mass of the molecule in [au]
        ! ground_pot    -   R->R function - potential of neutral molecule
        ! N_Four        -   total number of particle-in-a-box functions used in Four_DVR
        ! exc_pot       -   R->R function - potential of electronicaly excited molecule
        ! in_state      -   initial state [integer], 0 for the ground state
        ! n_states      -   total number of states of the neutral molecule potential calculated
        ! g_Vde         -   R->R function - radial component of the V_de
        ! f_Vde         -   R->R function - energy component of the V_de
        ! F_trans       -   R->C function - transform of the energy component of the V_de - F(e)
        ! cs_VE         -   cross sections for all final states
        ! cs_DA         -   cross sections for disociative attachment
        ! cs_DA_infty   -   cross sections for disociative attachment using magnitude of wavefunction at infinity
        ! cs_DA_inel    -   cross sections for disociative attachment using cs_DA = cs_Inel - cs_VE
        
        REAL(KIND = idk), INTENT(IN) :: e_i(:), mass, x_mesh(:)
        INTEGER, INTENT(IN) :: in_state, n_states, N_Four
        PROCEDURE(real_ve_interface), POINTER, INTENT(IN) :: exc_pot, g_Vde, f_Vde, ground_pot
        PROCEDURE(complex_ve_interface), POINTER, INTENT(IN) :: F_trans
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT) :: cs_VE(:,:), cs_VE_tot(:), cs_VE_bg(:,:), cs_VE_tot_bg(:), cs_DA(:), cs_DA_infty(:), cs_DA_inel(:)
        REAL(KIND = idk), ALLOCATABLE :: E_nu(:), Fnu(:,:), nu(:,:), gnu(:,:), fVec(:), Kfree(:), Kpot(:)
        COMPLEX(KIND = idk), ALLOCATABLE :: greenedgnu(:,:), G_matrix(:,:), auxVec(:), F_matrix(:,:), rhsVec(:), coeff_matrix(:,:), T_VEaux(:), Kprop(:)
        COMPLEX(KIND = idk) :: T_DAaux, auxnum
        REAL(KIND = idk) ::  E_out, delta_in, delta_out
        COMPLEX(KIND = idk) :: T_res, T_total, T_bg_vib
        COMPLEX(KIND = idk), ALLOCATABLE :: T_bg(:)
        INTEGER :: n_mesh
        REAL(KIND = idk) :: dx, E_nui, K_vec
        REAL(KIND = idk) :: V_inf
        REAL(KIND = idk), ALLOCATABLE :: E_cons(:)
        INTEGER :: i, j, k, l, i_state, nr_inf
        CHARACTER(LEN=50) :: full_message
        
        n_mesh = SIZE(x_mesh)
        dx = ABS(x_mesh(2)-x_mesh(1))
        i_state = in_state + 1
        V_inf = exc_pot(x_mesh(n_mesh))
        
        ! Computation of eigenstates of ground_pot
        CALL Four_DVR(x_mesh(1), x_mesh(n_mesh), ground_pot, mass, N_Four, n_states, E_nu, Fnu)
        CALL Inverse_F(x_mesh, Fnu, nu)
        DEALLOCATE(Fnu)
        E_nui = E_nu(i_state)
        ALLOCATE(E_cons(SIZE(e_i)))
        DO l = 1, SIZE(e_i)
            E_cons(l) = e_i(l) + E_nui
        END DO

        
        ! Construction of G_matrix
        ALLOCATE(gnu(n_states,n_mesh))
        DO i = 1, n_states
            DO j = 1, n_mesh
                gnu(i,j) = g_Vde(x_mesh(j)) * nu(i,j)               
            END DO
        END DO
        
        ALLOCATE(cs_VE(SIZE(e_i),n_states))
        ALLOCATE(cs_VE_tot(SIZE(e_i)))
        ALLOCATE(cs_VE_bg(SIZE(e_i),n_states))
        ALLOCATE(cs_VE_tot_bg(SIZE(e_i)))
        ALLOCATE(cs_DA(SIZE(e_i)))
        ALLOCATE(cs_DA_infty(SIZE(e_i)))
        ALLOCATE(cs_DA_inel(SIZE(e_i)))

        
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, greenedgnu, G_matrix, auxVec, F_matrix, rhsVec, coeff_matrix, T_VEaux, Kprop, T_DAaux, auxnum, nr_inf, i, j, k, &
        !$OMP full_message, Kfree, Kpot, fVec, K_vec, T_res, T_total, E_out, T_bg, T_bg_vib, delta_in, delta_out)
        !$OMP DO SCHEDULE(STATIC)
        DO l = 1, SIZE(e_i)

            ALLOCATE(greenedgnu(n_states,n_mesh))
            CALL apply_green(x_mesh, E_cons(l), mass, 0.0d0, 0, exc_pot, gnu, greenedgnu)
            
            ALLOCATE(G_matrix(n_states,n_states), auxVec(n_mesh))
        
            DO i = 1, n_states
                DO j = 1, n_states
                    DO k = 1, n_mesh
                        auxVec(k) = gnu(i,k)*greenedgnu(j,k)
                    END DO
                    CALL definite_integral(auxVec,dx,G_matrix(i,j))
                
                END DO
            END DO
        
            DEALLOCATE(greenedgnu,auxVec)
        
            ! Construction of F_matrix
            ALLOCATE(F_matrix(n_states,n_states))
            F_matrix = (0.0d0, 0.0d0)
            DO i = 1, n_states
                F_matrix(i,i) = F_trans(E_cons(l) - E_nu(i))
            END DO
            
            ! Solving [I-G_matrix F_matrix]T_VE/f(e_f) = G_matrix fVec
            ALLOCATE(coeff_matrix(n_states,n_states))
            coeff_matrix = (0.0d0, 0.0d0)
            DO i = 1, n_states
                coeff_matrix(i,i) = (1.0d0, 0.0d0)
            END DO
            coeff_matrix = coeff_matrix - MATMUL(G_matrix,F_matrix)
        
            ALLOCATE(fVec(n_states),rhsVec(n_states))
            fVec = 0.0d0 
            fVec(i_state) = f_Vde(e_i(l))
            rhsVec = MATMUL(G_matrix,fVec)
        
            DEALLOCATE(G_matrix,F_matrix,fVec)
        
            CALL solve_SLE(coeff_matrix,rhsVec,T_VEaux)
        
            DEALLOCATE(coeff_matrix,rhsVec)
            
            cs_VE_tot(l) = 0.0d0
            cs_VE_tot_bg(l) = 0.0d0
            
            ALLOCATE(T_bg(n_mesh))
            CALL compute_Tbg(x_mesh, e_i(l), e_i(l), T_bg, delta_in)
            DO i = 1, n_states
                E_out = E_cons(l) - E_nu(i)
                T_res = T_VEaux(i) * f_Vde(E_cons(l) - E_nu(i))
                
                CALL compute_Tbg(x_mesh, e_i(l), E_out, T_bg, delta_out)

                DO j = 1, n_mesh
                    T_bg(j) = T_bg(j) * nu(i_state,j) * nu(i,j)
                END DO
                CALL definite_integral(T_bg, dx, T_bg_vib)
                
                T_total = T_res * EXP(CMPLX(0.0d0, -(delta_out + delta_in), KIND=idk)) + T_bg_vib
                
                cs_VE(l,i) = 2.0d0 * PI**3 * (ABS(T_res))**2 / e_i(l)
                cs_VE_tot(l)= cs_VE_tot(l) + cs_VE(l,i) 
                
                cs_VE_bg(l,i) = 2.0d0 * PI**3 * (ABS(T_total))**2 / e_i(l)
                cs_VE_tot_bg(l) = cs_VE_tot_bg(l) + cs_VE_bg(l,i) 
            END DO
            DEALLOCATE(T_bg)
            
            
            IF (E_cons(l) > V_inf) THEN 
                
                K_vec = SQRT(2.0d0 * mass * (E_cons(l) - V_inf))
                     
                T_DAaux = 0.0d0 
                nr_inf = CEILING(n_mesh * 0.95d0)
                DO j = 1, n_states
                    T_DAaux = T_DAaux + T_VEaux(j)*nu(j,nr_inf)
                END DO
                T_DAaux = T_DAaux / g_Vde(x_mesh(nr_inf))
                cs_DA_infty(l) = PI**2 * K_vec / (mass * e_i(l)) * (ABS(T_DAaux))**2 
                
                ALLOCATE(Kfree(n_mesh),Kpot(n_mesh))
                DO i = 1, n_mesh
                    Kfree(i) = SIN(K_vec*x_mesh(i))
                    Kpot(i) = Kfree(i) * (exc_pot(x_mesh(i)) - V_inf)
                END DO
                Kfree = Kfree / SQRT(2/(PI * K_vec))
                Kpot = Kpot / SQRT(2/(PI * K_vec))
                
                ALLOCATE(Kprop(n_mesh))
                CALL apply_green(x_mesh, E_cons(l), mass, 0.0d0, 0, exc_pot, Kpot, Kprop)
                DO i = 1, n_mesh
                    Kprop(i) = Kprop(i) + Kfree(i)
                END DO
                DEALLOCATE(Kfree,Kpot)
                T_DAaux = 0.0d0 
                
                ALLOCATE(auxVec(n_mesh))
                DO i = 1, n_mesh
                    auxVec(i) = Kprop(i) * gnu(i_state,i)
                END DO
                CALL definite_integral(auxVec,dx,T_DAaux)
                T_DAaux = T_DAaux * f_Vde(e_i(l))
                
                DO j = 1, n_states
                    DO i = 1, n_mesh
                        auxVec(i) = Kprop(i) * gnu(j,i)
                    END DO
                    CALL definite_integral(auxVec,dx,auxnum)
                    T_DAaux = T_DAaux + F_trans(E_cons(l)-E_nu(j)) * auxnum * T_VEaux(j)
                END DO
                
                cs_DA(l) = 8 * PI * mass / (e_i(l) * K_vec**2) * (ABS(T_DAaux))**2 
                
                DEALLOCATE(auxVec,Kprop)
                
                cs_DA_inel(l) = -2.0d0 * PI**2 / e_i(l)  * AIMAG(T_VEaux(i_state)*f_Vde(E_cons(l)-E_nui)) - cs_VE_tot(l)
                
            ELSE 
                cs_DA(l) = 0.0d0
                cs_DA_infty(l) = 0.0d0
                cs_DA_inel(l) = 0.0d0
            END IF
            
            DEALLOCATE(T_VEaux)
            
            !$OMP CRITICAL
            WRITE(full_message, '(A,F7.4,A)') 'TVE calculated successfully for e =', e_i(l)*phys_h0, ' eV'
            CALL CONSOLE(full_message)
            !$OMP END CRITICAL
            
        END DO 
        !$OMP END DO
        !$OMP END PARALLEL
            
        DEALLOCATE(gnu, nu)
    
    END SUBROUTINE vib_exc


END MODULE VibrationalExcitation