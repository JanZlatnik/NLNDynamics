!-----------------------------------< RCModel2 >------------------------------------!
!                                                                                   !    
! Contains: Model functions constructed by Dr. Roman Curik                          !                                                                             
!                                                                                   !
! Last revision:    02/03/2026                                                      !
!                                                                                   !
!-----------------------------------------------------------------------------------!



MODULE RCModel2
    USE Parameters
    USE Math
    USE Green
    USE whittaker_w, ONLY: coulomb_whittaker
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: initiate_model, compute_Tbg
    PUBLIC :: V0, Vd, g, f, F_lvlshift, VLCP, GammaLCP
    PUBLIC :: mass, DS_R0


    REAL(KIND = idk), PARAMETER :: mass = 918.076337d0



    INTEGER, PARAMETER          :: DS_mp = 5000
    REAL(KIND = idk), PARAMETER :: DS_xmin = 1.0d-10 
    REAL(KIND = idk), PARAMETER :: DS_xmax = 25.0d0 
    REAL(KIND = idk), PARAMETER :: DS_xrange_SR = 13.0d0 
    REAL(KIND = idk), PARAMETER :: DS_mass = 1.0d0
    REAL(KIND = idk), PARAMETER :: DS_Z = 1.0d0
    INTEGER, PARAMETER          :: DS_l = 1
    REAL(KIND = idk), PARAMETER :: DS_R0 = 2.0d0

    REAL(KIND = idk), ALLOCATABLE :: DS_x(:)
    REAL(KIND = idk), ALLOCATABLE :: DS_dstate(:)
    REAL(KIND = idk)              :: DS_dx


    ! Interface for R->R potential V(x)
    ABSTRACT INTERFACE
        REAL(KIND = KIND(1.d0)) FUNCTION potential_interface(R)
            IMPLICIT NONE
            REAL(KIND = KIND(1.d0)), INTENT(IN) :: R
        END FUNCTION potential_interface
    END INTERFACE



    CONTAINS

    ! DS Interaction Potentials
    ! Potential V_2D_SHORT(R)
    REAL(KIND = idk) FUNCTION V_2D_SR(x,R)
        REAL(KIND = idk), INTENT(IN) :: x, R
        REAL(KIND = idk), PARAMETER :: IID_a = 0.4d0
        REAL(KIND = idk), PARAMETER :: IID_b = 2.0d0
        REAL(KIND = idk), PARAMETER :: IID_c = 1.5d0
        REAL(KIND = idk), PARAMETER :: IID_d = 0.72d0
        REAL(KIND = idk), PARAMETER :: IID_rc = 5.0d0
        REAL(KIND = idk), PARAMETER :: IID_rrc = 1.5d0
        V_2D_SR = IID_a * EXP(-(x-IID_rc)**2/IID_b**2) - IID_d * EXP(-x**2/4.0d0) * TANH( (R-IID_rrc)/IID_c )
    END FUNCTION V_2D_SR


    ! Potential V_2D(R)
    REAL(KIND = idk) FUNCTION V_2D(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        V_2D = - DS_Z / R + 0.5d0 * REAL(DS_l*(DS_l+1),KIND=idk) / R**2 + V_2D_SR(R, DS_R0)
    END FUNCTION V_2D


    ! Potential V_2D(R)
    REAL(KIND = idk) FUNCTION V_IID_asymptotic(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R
        REAL(KIND = idk), PARAMETER :: IID_a = 0.4d0
        REAL(KIND = idk), PARAMETER :: IID_b = 2.0d0
        REAL(KIND = idk), PARAMETER :: IID_c = 1.5d0
        REAL(KIND = idk), PARAMETER :: IID_d = 0.72d0
        REAL(KIND = idk), PARAMETER :: IID_rc = 5.0d0
        REAL(KIND = idk), PARAMETER :: IID_rrc = 1.5d0
        V_IID_asymptotic = - 1.0d0 / R + 0.5d0 * REAL(DS_l*(DS_l+1),KIND=idk) / R**2 + IID_a * EXP(-(R-IID_rc)**2/IID_b**2) - IID_d * EXP(-R**2/4.0d0)
    END FUNCTION V_IID_asymptotic








    ! V0 potential
    REAL(KIND = idk) FUNCTION V0(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R

        REAL(KIND = idk), PARAMETER :: D0 = 0.1027d0
        REAL(KIND = idk), PARAMETER :: alpha0 = 0.69
        REAL(KIND = idk), PARAMETER :: Re = 2.0d0

        V0 = EXP(-2.0d0 * alpha0 * (R-Re)) - 2.0d0 * EXP(-alpha0 * (R-Re))
        V0 = V0 * D0

    END FUNCTION V0


    ! Distcrete state energy
    REAL(KIND = idk) FUNCTION DSE(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R

        REAL(KIND = idk), PARAMETER :: DSE_const = 1.267048079198696d-01
        REAL(KIND = idk), PARAMETER :: DSE_A = -2.517156579364969d-01
        REAL(KIND = idk), PARAMETER :: DSE_alpha = 6.667134381849124d-01
        REAL(KIND = idk), PARAMETER :: DSE_Rc = 1.497541302278889d+00
        REAL(KIND = idk), PARAMETER :: DSE_B = -1.425252530184318d-03
        REAL(KIND = idk), PARAMETER :: DSE_beta = 3.645617517412014d-01
        REAL(KIND = idk), PARAMETER :: DSE_Rw = -2.545455471249359d-01

        DSE = DSE_const + DSE_A * TANH(DSE_alpha * (R - DSE_Rc)) + DSE_B * exp(-DSE_beta * (R - DSE_Rw)**2)

    END FUNCTION DSE


    ! Vd potential
    REAL(KIND = idk) FUNCTION Vd(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R

        Vd = V0(R) + DSE(R)

    END FUNCTION Vd


    ! LCP potential
    REAL(KIND = idk) FUNCTION VLCP(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R

        REAL(KIND = idk), PARAMETER :: LCP_Delta_A = 6.551393198713883d-02
        REAL(KIND = idk), PARAMETER :: LCP_Delta_alpha = 1.031154246905986d+00
        REAL(KIND = idk), PARAMETER :: LCP_Delta_Rc = 1.163254381883272d+00

        REAL(KIND = idk) :: Delta

        Delta = LCP_Delta_A * (TANH(LCP_Delta_alpha * (R - LCP_Delta_Rc)) - 1.0_idk)
        
        VLCP = Vd(R) + Delta

    END FUNCTION VLCP


    ! LCP width (Imaginary part)
    REAL(KIND = idk) FUNCTION GammaLCP(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R

        ! E_res parameters
        REAL(KIND = idk), PARAMETER :: LCP_E_res_const = 3.989745802406164d-02
        REAL(KIND = idk), PARAMETER :: LCP_E_res_A = -1.654452789870694d-01
        REAL(KIND = idk), PARAMETER :: LCP_E_res_alpha = 6.912177018334650d-01
        REAL(KIND = idk), PARAMETER :: LCP_E_res_Rc = 1.849820690923383d+00
        REAL(KIND = idk), PARAMETER :: LCP_E_res_B = 4.488097282210504d-03
        REAL(KIND = idk), PARAMETER :: LCP_E_res_beta = 3.206490451536629d-01

        ! Gamma parameters
        REAL(KIND = idk), PARAMETER :: LCP_Gamma_A = -7.443818187208139d-03
        REAL(KIND = idk), PARAMETER :: LCP_Gamma_B = 7.608497664604953d-02
        REAL(KIND = idk), PARAMETER :: LCP_Gamma_beta = 9.339636444257674d-01
        REAL(KIND = idk), PARAMETER :: LCP_Gamma_gamma_exp = -8.520190202347913d-02

        REAL(KIND = idk) :: E_res_H, E_res_eV

        E_res_H = LCP_E_res_const + LCP_E_res_A * TANH(LCP_E_res_alpha * (R - LCP_E_res_Rc)) + LCP_E_res_B * EXP(-LCP_E_res_beta * R)

        E_res_eV = E_res_H * phys_h0

        IF (E_res_eV > 0.0_idk) THEN
            GammaLCP = (LCP_Gamma_A + LCP_Gamma_B * EXP(-LCP_Gamma_beta * R)) * (E_res_eV**LCP_Gamma_gamma_exp)
        ELSE
            GammaLCP = 0.0_idk
        END IF

    END FUNCTION GammaLCP







    ! g(R) function
    REAL(KIND = idk) FUNCTION g(R)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: R

        REAL(KIND = idk), PARAMETER :: g_alpha = 1.145617054714233d+00
        REAL(KIND = idk), PARAMETER :: g_Rc = 1.782092739061193d+00
        REAL(KIND = idk), PARAMETER :: g_A_g = 7.424888821139934d+00
        REAL(KIND = idk), PARAMETER :: g_beta = 3.082516174633375d-01
        REAL(KIND = idk), PARAMETER :: g_Rw = -6.162091666248357d-01
        REAL(KIND = idk), PARAMETER :: g_norm = 6.039593778762647d-01

        g = 1.0_idk - TANH(g_alpha * (R - g_Rc)) + g_A_g * EXP(-g_beta * (R - g_Rw)**2)
        g = SQRT(g_norm * g)

    END FUNCTION g


    ! F_trans function
    COMPLEX(KIND = idk) FUNCTION F_lvlshift(e)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: e

        PROCEDURE(potential_interface), POINTER :: DSV_ptr
        COMPLEX(KIND = idk), ALLOCATABLE :: green_dstate(:)
        COMPLEX(KIND = idk) :: greendot

        DSV_ptr => V_2D

        ALLOCATE(green_dstate(DS_mp))
        CALL apply_green(DS_x, e, DS_mass, DS_Z, DS_l, DSV_ptr, DS_dstate, green_dstate)
        green_dstate(:) = green_dstate(:) * DS_dstate(:)
        CALL definite_integral(green_dstate, DS_dx, greendot)
        DEALLOCATE(green_dstate)

        F_lvlshift = e - DSE(DS_R0) - 1.0d0 / greendot


    END FUNCTION F_lvlshift



    ! f function
    REAL(KIND = idk) FUNCTION f(e)
        IMPLICIT NONE
        REAL(KIND = idk), INTENT(IN) :: e

        IF (e > 0.0d0) THEN
            f = SQRT(-AIMAG(F_lvlshift(e))/pi)
        ELSE
            f = 0.0d0
        END IF

    END FUNCTION f




    ! Background scattering
    SUBROUTINE compute_Tbg(R, Ein, Eout, Tbg, phaseout)
        REAL(KIND = idk), INTENT(IN) :: R(:), Ein, Eout
        COMPLEX(KIND = idk), INTENT(OUT) :: Tbg(:)
        REAL(KIND = idk), INTENT(OUT) :: phaseout

        INTEGER :: i, j
        PROCEDURE(potential_interface), POINTER :: DSV_ptr
        REAL(KIND = idk), ALLOCATABLE :: freestate_in(:), freestate_out(:)
        REAL(KIND = idk), ALLOCATABLE :: mat(:,:)
        COMPLEX(KIND = idk), ALLOCATABLE :: greened_mat(:,:), green_dstate(:), scattered_state(:)
        COMPLEX(KIND = idk), ALLOCATABLE :: phi_bg_R0(:), integrand(:), tmp_vec(:)
        COMPLEX(KIND = idk) :: greendot, scattdot, Vde_out_R0
        COMPLEX(KIND = idk) :: term1, term2, term2_base
        REAL(KIND = idk) :: R_val

        IF (Eout <= 0.0d0) THEN
            Tbg = (0.0d0, 0.0d0)
            phaseout = 0.0d0
            RETURN
        END IF

        ALLOCATE(freestate_in(DS_mp), freestate_out(DS_mp))
        ALLOCATE(mat(2,DS_mp), greened_mat(2,DS_mp))

        DO i = 1, DS_mp
            freestate_in(i)  = general_free_solution(DS_x(i), Ein, DS_mass, DS_Z, DS_l)
            freestate_out(i) = general_free_solution(DS_x(i), Eout, DS_mass, DS_Z, DS_l)
            mat(2,i) = freestate_out(i) * V_2D_SR(DS_x(i), DS_R0)
        END DO
        mat(1,:) = DS_dstate(:)

        DSV_ptr => V_2D
        CALL apply_green(DS_x, Eout, DS_mass, DS_Z, DS_l, DSV_ptr, mat, greened_mat)
        DEALLOCATE(mat)

        ALLOCATE(green_dstate(DS_mp), scattered_state(DS_mp), tmp_vec(DS_mp))
        DO i = 1, DS_mp
            green_dstate(i) = greened_mat(1,i) * DS_dstate(i)
        END DO
        CALL definite_integral(green_dstate, DS_dx, greendot)
        greendot = 1.0d0 / greendot
        DEALLOCATE(green_dstate)

        scattered_state = freestate_out + greened_mat(2,:)
        DEALLOCATE(freestate_out)
        DO i = 1, DS_mp
            tmp_vec(i) = scattered_state(i) * DS_dstate(i)
        END DO
        CALL definite_integral(tmp_vec, DS_dx, scattdot)

        Vde_out_R0 = scattdot * greendot
        phaseout = ATAN2(AIMAG(Vde_out_R0), REAL(Vde_out_R0))

        ALLOCATE(phi_bg_R0(DS_mp))
        DO i = 1, DS_mp
            phi_bg_R0(i) = scattered_state(i) - Vde_out_R0 * greened_mat(1,i)
        END DO
        DEALLOCATE(greened_mat, scattered_state)

        DO i = 1, DS_mp
            tmp_vec(i) = DS_dstate(i) * freestate_in(i)
        END DO
        CALL definite_integral(tmp_vec, DS_dx, term2_base)
        DEALLOCATE(tmp_vec)

        ALLOCATE(integrand(DS_mp))
        DO j = 1, SIZE(R)
            R_val = R(j)

            DO i = 1, DS_mp
                integrand(i) = phi_bg_R0(i) * V_2D_SR(DS_x(i),R_val) * freestate_in(i)
            END DO
            CALL definite_integral(integrand, DS_dx, term1)

            term2 = Vde_out_R0 * g(R_val) * term2_base

            Tbg(j) = term1 - term2
        END DO

        DEALLOCATE(freestate_in, phi_bg_R0, integrand)


    END SUBROUTINE compute_Tbg








    ! Initiation of DS grid and discrete state
    SUBROUTINE initiate_model()
        IMPLICIT NONE
        
        REAL(KIND = idk) :: E_asympt
        INTEGER :: i, N, iounit
        PROCEDURE(potential_interface), POINTER :: Vasym_ptr
        REAL(KIND = idk), ALLOCATABLE :: psi_tmp(:), normstate(:)
        
        REAL(KIND = idk) :: norm, norm_tail, k_num, eta, R_match, w_match, w_deriv
        INTEGER :: sf_match, sf_curr
        REAL(KIND = idk) :: A_scale_sq, r_curr, w_curr, term, term_prev, dr_adaptive
        REAL(KIND = idk) :: R_turn, E_kin_local, V_coul, local_k
        LOGICAL :: tail_converged
        REAL(KIND = idk), PARAMETER :: eps = 1.0d-30
        CHARACTER(LEN=256) :: message

        IF (.NOT. ALLOCATED(DS_x)) ALLOCATE(DS_x(DS_mp))
        CALL make_mesh(DS_xmin, DS_xmax, DS_mp, DS_x, DS_dx)
        N = SIZE(DS_x)

        E_asympt = DSE(1.0d10)

        Vasym_ptr => V_IID_asymptotic
        CALL calc_H_wavefunction(DS_x, E_asympt, DS_mass, DS_l, Vasym_ptr, psi_tmp)

        ALLOCATE(normstate(N))
        DO i = 1, N
            normstate(i) = psi_tmp(i)**2
        END DO
        CALL definite_integral(normstate, DS_dx, norm)
        DEALLOCATE(normstate)

        k_num = SQRT(2.0d0 * DS_mass * ABS(E_asympt)) / hbar
        eta = - DS_Z * DS_mass / (hbar**2 * k_num)
        R_match = DS_x(N)
        
        CALL coulomb_whittaker(eta, DS_l, k_num * R_match, w_match, w_deriv, sf_match)
        A_scale_sq = (psi_tmp(N) / w_match) ** 2
        
        R_turn = DS_Z / ABS(E_asympt)
        norm_tail = 0.0d0
        r_curr = R_match
        tail_converged = .FALSE.
        term_prev = w_match ** 2
        dr_adaptive = DS_dx
        
        DO WHILE (.NOT. tail_converged)
            V_coul = - DS_Z / r_curr
            E_kin_local = E_asympt - V_coul

            IF (E_kin_local > 0.0d0) THEN
                local_k = SQRT(2.0d0 * DS_mass * E_kin_local) / hbar
                IF(local_k > 1.0d-10) THEN
                    dr_adaptive = 2.0d0 * PI / local_k / 1000.0d0
                ELSE
                    dr_adaptive = dr_adaptive * 1.1d0
                END IF
            ELSE
                dr_adaptive = dr_adaptive * 1.05d0
            END IF
                
            dr_adaptive = MAX(dr_adaptive, DS_dx)
            r_curr = r_curr + dr_adaptive

            CALL coulomb_whittaker(eta, DS_l, k_num * r_curr, w_curr, w_deriv, sf_curr)
            term = w_curr ** 2 * (10.0d0) ** (2*(sf_curr - sf_match))
            norm_tail = norm_tail + 0.5d0 * (term + term_prev) * dr_adaptive
            term_prev = term
            
            IF (r_curr > MAX(R_turn, R_match) ) THEN 
                IF (term < eps * norm_tail) tail_converged = .TRUE.
            END IF
            
            IF (r_curr > 100.0d0 * MAX(R_turn, R_match)) THEN
                CALL CONSOLE('[ERROR]: Normalization of dstate tail did not converge!')
                STOP
            END IF
        END DO
        
        norm_tail = norm_tail * A_scale_sq
        norm = SQRT(ABS(norm + norm_tail))
        
        IF (ALLOCATED(DS_dstate)) DEALLOCATE(DS_dstate)
        ALLOCATE(DS_dstate(N))
        DO i = 1, N
            DS_dstate(i) = psi_tmp(i) / norm
        END DO
        DEALLOCATE(psi_tmp)

        WRITE(message, '(A,G0.8,A)') ' -> dstate initialized directly at energy E = ', E_asympt*phys_h0, ' eV'
        CALL CONSOLE(message)

        OPEN(NEWUNIT=iounit, FILE='MODEL/dstate.txt', STATUS='replace', ACTION='write')
        WRITE(iounit, '(A3, A17, A21)') '#', 'R [au]','dstate [a.u.]'
        DO i = 1, N
            WRITE(iounit, '(2E20.12)') DS_x(i), DS_dstate(i)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Data successfully written to MODEL/dstate.txt')

    END SUBROUTINE initiate_model




    SUBROUTINE calc_H_wavefunction(x, E, m, l, potential, psi)
        REAL(KIND = idk), INTENT(IN)                                :: x(:), E, m
        INTEGER, INTENT(IN)                                         :: l
        PROCEDURE(potential_interface), POINTER, INTENT(IN)         :: potential
        REAL(KIND = idk), ALLOCATABLE, INTENT(OUT)                  :: psi(:)
        
        INTEGER :: N, i, idx_match, sf_match, sf_curr
        REAL(KIND = idk) :: h, k_num, eta, w_match, w_curr, w_deriv, scale_base
        REAL(KIND = idk), ALLOCATABLE :: k2(:)
        
        IF(ALLOCATED(psi)) DEALLOCATE(psi)
        N = SIZE(x)
        ALLOCATE(psi(N)) 
        
        h = ABS(x(2)-x(1))
        idx_match = NINT((DS_xrange_SR - x(1)) / h) + 1
        
        IF (idx_match > N) idx_match = N
        IF (idx_match < 2) idx_match = 2 

        ALLOCATE(k2(idx_match))
        DO i = 1, idx_match
            k2(i) = 2.0d0*m/hbar**2 * (E - potential(x(i)))
        END DO
        
        psi(1) = 0.0d0
        psi(2) = h**(l+1)
        
        DO i = 2, idx_match - 1
            psi(i+1) = (2.0d0*(1.0d0-5.0d0*h**2*k2(i)/12.0d0)*psi(i) - (1.0d0+h**2*k2(i-1)/12.0d0)*psi(i-1)) / (1.0d0+h**2*k2(i+1)/12.0d0)
        END DO
        
        DEALLOCATE(k2) 
        
        k_num = SQRT(2.0d0 * m * ABS(E)) / hbar
        eta = -DS_Z * m / (hbar**2 * k_num) 
        CALL coulomb_whittaker(eta, l, k_num * x(idx_match), w_match, w_deriv, sf_match)
        scale_base = psi(idx_match) / w_match
        DO i = idx_match + 1, N
            CALL coulomb_whittaker(eta, l, k_num * x(i), w_curr, w_deriv, sf_curr)
            psi(i) = scale_base * w_curr * (10.0d0)**(sf_curr - sf_match)
        END DO
    
    END SUBROUTINE calc_H_wavefunction



END MODULE RCModel2