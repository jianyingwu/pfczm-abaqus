!**********************************************************************************************************
!                                                                                                         !      
!                                             UMAT for PF-CZM                                             !
!                                                                                                         !
!**********************************************************************************************************
*  Copyright (C) 2019 South China University of Technology, China. All rights reserved.
*  
*  This subroutine implemens the phase field regularized coheisve zone model
*
*  Status: only 2D plane stress and CPS4 elements are considered.
*  
*  Author: J. Y. Wu (jywu@scut.edu.cn) and Y. Huang
*  Date: 31 Oct. 2019
*
*
*  If you use this subroutine (research ONLY), you are also asked to cite our papers:
*  1. Wu, J. Y., 2017. A unified phase-field theory for the mechanics of damage and quasi-brittle failure. 
*     Journal of the Mechanics and Physics of Solids, 103: 72-99.
*  2. Wu, J. Y., 2018. A geometrically regularized gradient damage model with energetic equivalence. 
*     Computer Methods in Applied Mechanics and Engineering, 328: 612-637.
*  3. Wu, J. Y. and Nguyen, V.P., 2018. A length scale insensitive phase-field damage model for brittle 
*     fracture. Journal of the Mechanics and Physics of Solids, 119: 20-42.
*  4. Wu, J. Y., Huang, Y. and Nguyen, V. P., 2019. On the BFGS monolithic algorithm for the unified 
*     phase-field damage theory. Computer Methods in Applied Mechanics and Engineering, 112704.
*  5. Wu, J. Y. and Huang, Y., 2019. Comprehensive ABAQUS implementation of phase-field damage models 
*     for fracture in solids. Theoretical and Applied Fracutre Mechancis, in press.
*      
!**********************************************************************************************************
!
      module NumKind
!
!**********************************************************************************************************
        implicit none
        integer (kind(1)), parameter :: ikind = kind(1), 
     &                                  rkind = kind(0.D0), 
     &                                  lkind = kind(.true.)
!       
      end module Numkind

!**********************************************************************************************************
!
      module ModelParam
!
!**********************************************************************************************************
        use NumKind
        implicit none
        
        real (rkind) :: EA, nu, ft, Gf, lb
        real (rkind) :: p, a1, a2, a3
        real (rkind) :: De(3, 3) 
        integer (ikind) :: is
        real(rkind), parameter :: c0 = 3.1415926535897932384626433832d0
        
        ! Flag of initilization
        logical (lkind) :: bInitialized = .false.
        !logical (lkind) :: bOutput

        !logical (lkind) :: btemp
        
        contains
          !=============================================
          subroutine Initialize(props, nprops)
          !=============================================
            integer (ikind), intent (in) :: nprops
            real    (rkind), intent (in) :: props(nprops)
            !********************************************
            real(rkind) :: G0, K11, K12

            EA = props(1)       ! Young's modulus
            nu = props(2)       ! Poisson's ratio
            ft = props(3)       ! failure strength
            Gf = props(4)       ! fracture energy
            lb = props(5)       ! length scale
            is = nint(props(6)) ! softening type

            ! elastic stiffness matrix
            G0      = EA / (2.d0*(1.d0 + nu))
            K11     = EA / (1.D0 - nu * nu)
            K12     = nu * K11
            De(:,1) = (/ K11,  K12, 0.D0/)
            De(:,2) = (/ K12,  K11, 0.D0/)
            De(:,3) = (/0.D0, 0.D0,   G0/)

            ! softening parameters
            a1   =  4.d0/(c0*lb)*EA*Gf/(ft*ft)
            if      (is == 1) then  ! linear softening
              p  =  2.d0
              a2 = -0.5d0
              a3 =  0.0d0
            else if (is == 2) then  ! exponential softening
              p  =  2.5d0
              a2 =  2.0d0**(5.d0/3.d0) - 3.0d0
              a3 =  0.0d0
            else if (is == 3) then  ! blinear softening
              p  =  2.0d0
              a2 =  0.03687d0
              a3 =  20.8343d0
            else if (is == 4) then  ! concrete softening
              p  =  2.0d0
              a2 =  1.3868d0
              a3 =  0.6567d0
            else if (is == 5) then  ! hyperbolic softening
              p  =  4.0d0 
              a2 =  2.0d0**(7.d0/3.d0) - 4.5d0
              a3 =  0.0d0
            else
              write (*,*) '**error: Softening law No. ', is, 
     &                    'does not exist!'
            end if
     
            bInitialized = .true.
            
            return
          
          end subroutine Initialize

      end module ModelParam

!**********************************************************************************************************
!
      subroutine energeticFunc(omega,domega,ddomega,phi)
!
!**********************************************************************************************************
        use NumKind
        use ModelParam
        implicit none
      
        ! local varibles
        real(rkind) :: omega, domega, ddomega, phi
        real(rkind) :: fac1, dfac1, ddfac1, fac2, dfac2, ddfac2
      
        fac1    =  (1.d0 - phi)**p
        dfac1   = -p*(1.d0 - phi)**(p - 1.d0); 
        ddfac1  =  p*(p - 1.d0)*(1.d0 - phi)**(p - 2.d0)
        
        fac2    =  fac1   + a1*phi + a1*a2*phi**2.d0 + a1*a2*a3*phi**3.d0
        dfac2   =  dfac1  + a1 + 2.d0*a1*a2*phi + 3.d0*a1*a2*a3*phi**2.d0
        ddfac2  =  ddfac1 + 2.d0*a1*a2 + 6.d0*a1*a2*a3*phi
        
        omega   =  fac1/fac2        
        domega  =  (dfac1*fac2  - fac1*dfac2)/(fac2**2.d0)
        ddomega = ((ddfac1*fac2 - fac1*ddfac2)*fac2 - 2.d0*
     &             (dfac1*fac2 - fac1*dfac2)*dfac2)/(fac2**3.d0)
     
        return
      end subroutine energeticFunc
      
!**********************************************************************************************************
!
      subroutine geometricFunc(dalpha,ddalpha,phi)
!
!**********************************************************************************************************
        use NumKind
        implicit none
        
        real(rkind) :: dalpha, phi, ddalpha
        
        dalpha  = 2.d0 - 2.d0*phi
        ddalpha =-2.d0
        
        return 
      end subroutine geometricFunc  
      
      
!****************************************************************************************
      subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl,
     &                ddsddt, drplde, drpldt, stran, dstran, time, 
     &                dtime, temp, dtemp, predef, dpred, cmname, 
     &                ndi, nshr, ntens, nstatv, props, nprops, 
     &                coords, drot, pnewdt, celent, dfgrd0, dfgrd1, 
     &                noel, npt, layer, kspt, kstep, kinc)
!****************************************************************************************

        use NumKind
        use ModelParam
        implicit none
        
        integer (ikind) :: ndi, nshr, ntens, nstatv, nprops
        integer (ikind) :: noel, npt, layer
        integer (ikind) :: kspt, kstep(4), kinc
        
        real (rkind) :: stress(ntens), statev(nstatv), 
     &                  ddsdde(ntens,ntens), ddsddt(ntens),
     &                  drplde(ntens), stran(ntens), dstran(ntens), 
     &                  props(nprops)
        real (rkind) :: coords(3), drot(3), time(2), dfgrd0(3,3), 
     &                  dfgrd1(3,3)
        real (rkind) :: predef, dpred, sse, spd, scd, rpl, drpldt, 
     &                  dtime, temp, dtemp, pnewdt, celent
        character*8 cmname
        
        
      !****************************************************************************************    
        real (rkind) :: damage, omega, omega1, omega2, ps1,  
     &                  H, H0, H1, smax
        real (rkind) :: stressEff(ntens), ps(3)     
        write (*, *) 'TIME = ', time(2)

        ! Initilize
        if (.not. bInitialized) call Initialize(props, nprops)
        
      !****************************************************************************************
        !calculate the major principal stress and direction
        stressEff = matmul(De, stran+dstran)
        call SPRINC(stressEff, ps, 1, ndi, nshr)
        smax = max(ps(1), ps(2), ps(3))

        H0 = statev(1)
        H1 = 0.5d0 * max(smax, ft)**2 / EA
        H  = max(H1, H0)
        statev(1) = H        
           
        !
        damage = temp + dtemp
        call energeticFunc(omega,omega1,omega2, damage)
          
        !
        stress = omega  * stressEff
        ddsdde = omega  * De
        !ddsddt = omega1 * stressEff
        
        return

      end subroutine umat
      
      
!****************************************************************************************
      subroutine hetval(cmname, temp, time, dtime, statev, flux,
     &                  predef, dpred)
!****************************************************************************************
        use NumKind
        use ModelParam
        implicit none
                
        real (rkind) :: temp(2), time(2), flux(2), dtime
        real (rkind) :: statev(*), predef(*), dpred(*)
        character*8 cmname

        real (rkind) :: damage, omega, omega1, omega2, alpha1, alpha2, 
     &                  c1, c2, H

        damage  =  temp(1)
        call energeticFunc(omega,omega1,omega2,damage)
        call geometricFunc(alpha1,alpha2,damage)

        c1 = 0.5d0 *  c0 / (lb * Gf)
        c2 = 0.5d0 / (lb * lb)

        !
        H  = statev(1)
        flux(1) = -(c1 * omega1 * H + c2 * alpha1)
        flux(2) = -(c1 * omega2 * H + c2 * alpha2)
       
        return
      end subroutine hetval 
