!**********************************************************************************************************
!                                                                                                         !      
!                                             UEL for PF-CZM                                              !
!                                alternate minimization (staggered) solver                                !
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

        ! Flag of initilization
        logical (lkind) :: bInitialized = .false.

        ! Tolerance
        real    (rkind), parameter :: TOL = 1.0d-12 
        ! number of guass points
        integer (ikind), parameter :: ngp = 2
        
        ! geometric function parameter
        real(rkind), parameter :: c0 = 3.1415926535897932384626433832d0
        
        real (rkind) :: tc
        integer (ikind) :: isw, je

        ! 
        real(rkind) :: thk, EA, nu, Gf, ft, lb
        real(rkind) :: p, a1, a2, a3
        real(rkind) :: gp(ngp), gw(ngp)
        real(rkind) :: De(3, 3), QQ(16,16)

        !
        contains

        !===================================
          subroutine Initialize(props, nprops, istype, jelem)
        !===================================

            integer (ikind), intent (in) :: nprops, istype, jelem
            real    (rkind), intent (in) :: props(nprops)

            !********************************************
            real   (rkind) :: G0, K11, K12
            integer(ikind) :: indexq(16), i

            je = jelem

            ! material properties
            EA  = props(1) ! props(1) -- Young's modulus
            nu  = props(2) ! props(2) -- Poisson's ratio
            ft  = props(3) ! props(3) -- failure strength
            Gf  = props(4) ! props(4) -- fracture energy
            lb  = props(5) ! props(5) -- length scale            
            ! thickness
            thk = props(6) ! props(6) -- thickness
            if (thk < TOL) thk = 1.0
          
            ! elastic stiffness matrix
            G0      = EA / (2.d0*(1.d0 + nu))
            K11     = EA / (1.D0 - nu * nu)
            K12     = nu * K11
            De(:,1) = (/ K11,  K12, 0.D0/)
            De(:,2) = (/ K12,  K11, 0.D0/)
            De(:,3) = (/0.D0, 0.D0,   G0/)
            
            ! softening parameters
            a1   =  4.d0/(c0*lb)*EA*Gf/(ft*ft)
            if      (istype == 1) then  ! linear softening
              p  =  2.d0
              a2 = -0.5d0
              a3 =  0.0d0
            elseif (istype == 2) then  ! exponential softening
              p  =  2.5d0
              a2 =  2.0d0**(5.d0/3.d0) - 3.0d0
              a3 =  0.0d0
            elseif (istype == 3) then  ! blinear softening
              p  =  2.0d0
              a2 =  0.03687d0
              a3 =  20.8343d0
            elseif (istype == 4) then  ! concrete softening
              p  =  2.0d0
              a2 =  1.3868d0
              a3 =  0.6567d0
            elseif (istype == 5) then  ! hyperbolic softening
              p  =  4.0d0 
              a2 =  2.0d0**(7.d0/3.d0) - 4.5d0
              a3 =  0.0d0
            else
              write (*,*) '**error: Softening law No. ', istype, 
     &                    'does not exist!'
            endif
            
            ! integration
            gp = (/ -1.d0, 1.d0 /) / dsqrt(3.d0)
            gw = (/  1.d0, 1.d0 /)
            
            ! dof interchange
            indexq = (/ 1,2,13,9, 3,4,14,10, 5,6,15,11, 7,8,16,12 /)
            ! interchange the locations of dofs
            QQ = 0.d0
            do i = 1, 16
              QQ(indexq(i),i) = 1.d0
            enddo 
            
            bInitialized = .true.
            
            return
          end subroutine Initialize
      !========================================================================= 
      end module ModelParam


!**********************************************************************************************************
!
      module FEM
!
!**********************************************************************************************************
        use NumKind
        implicit none

        contains      
          !==================shape function and its derivative with xi and eta======================    
          subroutine shapefuc(n, dn_xieta, xi, eta)
          
            implicit none      
            real(rkind) :: n(4), dn_xieta(2, 4), xi, eta

            n(1) = 0.25d0*(1.d0 - xi)*(1.d0 - eta)
            n(2) = 0.25d0*(1.d0 + xi)*(1.d0 - eta)
            n(3) = 0.25d0*(1.d0 + xi)*(1.d0 + eta)
            n(4) = 0.25d0*(1.d0 - xi)*(1.d0 + eta)
            
            dn_xieta(1, 1) = -0.25d0*(1.d0 - eta)
            dn_xieta(1, 2) =  0.25d0*(1.d0 - eta)
            dn_xieta(1, 3) =  0.25d0*(1.d0 + eta)
            dn_xieta(1, 4) = -0.25d0*(1.d0 + eta)
            
            dn_xieta(2, 1) = -0.25d0*(1.d0 - xi)
            dn_xieta(2, 2) = -0.25d0*(1.d0 + xi)
            dn_xieta(2, 3) =  0.25d0*(1.d0 + xi)
            dn_xieta(2, 4) =  0.25d0*(1.d0 - xi)
            
            return 
          end subroutine shapefuc
          !===============traditional b matrix==============================================      
          subroutine b_matrix(nd,bd,b,det_jacb, coords,xi,eta)
          
            implicit none
            real(rkind) :: nd(4), bd(2,4), b(3,8)
            real(rkind) :: jacb(2,2), inv_jacb(2,2), coords(2, 4)
            real(rkind) :: det_jacb, xi, eta
            
            !local varibles
            real(rkind) :: n(4), dn_xieta(2,4), dn_x(4), dn_y(4)
            integer(ikind) :: i, j
             
            ! shape functions 
            call shapefuc(n,dn_xieta,xi,eta)
            nd = n
            
            ! jacob matrix
            jacb = matmul(dn_xieta, transpose(coords))            
            det_jacb = jacb(1,1)*jacb(2,2) - jacb(1,2)*jacb(2,1)
            inv_jacb(1, 1) = jacb(2, 2)
            inv_jacb(1, 2) =-jacb(1, 2)
            inv_jacb(2, 1) =-jacb(2, 1)
            inv_jacb(2, 2) = jacb(1, 1)
            inv_jacb = 1.d0/det_jacb*inv_jacb            
            
            !initialize varibles
            b = 0.d0
            do i = 1,4
              dn_x(i) = inv_jacb(1,1)*dn_xieta(1,i)
     &                + inv_jacb(1,2)*dn_xieta(2,i)
              dn_y(i) = inv_jacb(2,1)*dn_xieta(1,i)
     &                + inv_jacb(2,2)*dn_xieta(2,i)
            enddo
            
            do j = 1, 4
              b(1, 2*(j-1) + 1) = dn_x(j)
              b(2, 2*(j-1) + 2) = dn_y(j)
              b(3, 2*(j-1) + 1) = dn_y(j)
              b(3, 2*(j-1) + 2) = dn_x(j)
            enddo
            
            do j = 1,4
              bd(1,j) = dn_x(j)
              bd(2,j) = dn_y(j)
            enddo
          
            return
          end subroutine b_matrix

        !********************************************************************
        ! define the dyadic function
          function dyadic(vector1,vector2, vlen)
        !********************************************************************
            integer (ikind) :: vlen, i, j
            real    (rkind) :: vector1(vlen),vector2(vlen)
            real    (rkind) :: dyadic(vlen,vlen)
          
            do i = 1, vlen
              do j = 1, vlen
                dyadic(i,j) = vector1(i) * vector2(j)
              enddo
            enddo

            return
          end function dyadic

      end module FEM

!**********************************************************************************************************
!
      subroutine pfczm(fk,ak, coords,u,svars)     
!
!**********************************************************************************************************
        use NumKind
        use ModelParam
        use FEM
        implicit none

        real(rkind) :: fk(16), ak(16,16), u(16), coords(2,4), svars(4)
       
        ! local varibles
        real(rkind):: b(3,8), nd(4), bd(2,4)
        real(rkind):: uu(8), dd(4), rd(4), ru(8), kdd(4,4), kuu(8,8)
        real(rkind):: rr(16), kk(16,16)
        real(rkind):: strain(3), stressEff(3)
        real(rkind):: det_jacb, energy_pos, phi, omega, domega, ddomega
        real(rkind):: dalpha, ddalpha, phi_source, dphi_source, dvol
        integer(ikind):: i, j, k

        real(rkind):: savg, sdif, sdev, smax, smin
        
        do i = 1, 4
          uu(2*i - 1 ) = u(4*i  - 3)
          uu(2*i     ) = u(4*i  - 2)
          dd(i       ) = u(4*i     )
        enddo
        
        ! initialize varibles
        rd  = 0.d0
        kdd = 0.d0
        kuu = 0.d0
        do i = 1, ngp
          do j = 1, ngp      
            call b_matrix(nd,bd,b,det_jacb, coords,gp(i),gp(j))
              
            strain = matmul(b, uu)  ! strain field
            stressEff = matmul(De, strain) ! effective stress

c           max/min pricipal stress
            savg = 0.5*(stressEff(1) + stressEff(2))
            sdif = 0.5*(stressEff(1) - stressEff(2))
            sdev = sqrt(sdif*sdif + stressEff(3)*stressEff(3))
            smax = savg + sdev
            smin = savg - sdev

            ! crack driving force
            k = (i - 1) * 2 + j
            energy_pos = 0.5d0*max(smax, ft)**2/EA
            energy_pos = max(energy_pos, svars(k))
            svars(k) = energy_pos
            
            phi  = dot_product(nd,dd) ! crack phase-field
            call geometricFunc(dalpha,ddalpha,phi)
            call energeticFunc(omega,domega,ddomega,phi)

            phi_source  = domega*energy_pos  + Gf/(c0*lb)*dalpha
            dphi_source = ddomega*energy_pos + Gf/(c0*lb)*ddalpha

            ! residual for damage
            dvol =  gw(i)*gw(j)*det_jacb*thk
            rd   =  rd  + dvol*(phi_source*nd + 2.d0*lb*Gf/c0
     3           *  matmul(transpose(bd), matmul(bd, dd)))

            kuu  =  kuu + dvol*matmul(matmul(transpose(b), omega*De), b)

            kdd  =  kdd + dvol*(dphi_source*dyadic(nd, nd, 4)
     &           +  2.d0*lb*Gf/c0*matmul(transpose(bd),bd))            
          enddo
        enddo
        ru = matmul(kuu,uu) ! applies to hybrid formulation
        
        rr = 0.d0     
        kk = 0.d0
        rr(1:8 ) = ru
        rr(9:12) = rd        
        kk(1:8 , 1:8 ) = kuu
        kk(9:12, 9:12) = kdd
        do i = 13, 16
          kk(i,i) = 1.d3     ! stiffness for the dummy dofs
        enddo
      
        ! interchange the locations of various
        fk = matmul(transpose(QQ), rr)
        ak = matmul(matmul(transpose(QQ),k k), QQ)
  
        return 
      end subroutine pfczm
      
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
     
!**********************************************************************************************************
      subroutine UEL(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars,
     &               props, nprops, coords, mcrd, nnode, 
     &               u, du, v, a, jtype, time, dtime, kstep, 
     &               kinc, jelem, params, ndload, jdltyp, adlmag,
     &               predef, npredf, lflags, mlvarx, ddlmag, mdload,
     &               pnewdt, jprops,njprop,period)
!**********************************************************************************************************

        use NumKind
        use ModelParam
        implicit none

!**********************************************************************************************************
      ! interface of uel, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! variables passed in
        integer (ikind), intent (in    ) :: ndofel, mlvarx, nrhs, 
     &    nsvars, nprops, mcrd, nnode, jtype, kstep, kinc, jelem, 
     &    ndload,npredf, mdload, njprop
     
        integer (ikind), intent (in    ) :: jdltyp(mdload,*), 
     &    lflags(*), jprops(njprop)
     
        real    (rkind), intent (in    ) :: props(nprops), 
     &    coords(mcrd,nnode), u(ndofel), du(mlvarx,*), v(ndofel), 
     &    a(ndofel), time(2), params(3), adlmag(mdload,*),
     &    ddlmag(mdload,*), predef(2,npredf,nnode), dtime, period
     
        ! variables can be updated
        real    (rkind), intent (in out) :: pnewdt
  
        ! variables to be updated (the update of energy(8) is optional)
        real    (rkind), intent (in out) :: rhs(mlvarx,nrhs), 
     &    amatrx(ndofel,ndofel), svars(nsvars), energy(8)
      ! interface of uel, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !********************************************************************************************************
      !                   
      ! user coding to define rhs, amatrx, svars, energy and pnewdt (optional for the last two)
      !

         real (rkind) :: ak(16,16), fk(16)
         write (*, *) 'TIME = ', time(2)

        ! initialize parameters, etc.                     
        if (.not. bInitialized)
     2    call Initialize(props, nprops, jtype, jelem)
        
        if (jelem == je) then
          if (time(2) > tc) then
            tc  = time(2)
            isw = 1
          endif
          isw = -isw
        endif
        
        ! right now only Q4 element is implemented
        call pfczm(fk,ak, coords,u,svars)
     
        rhs = 0.d0
        amatrx = 0.d0
        rhs(1:16, 1) = rhs(1:16, 1) - fk
        amatrx = amatrx + ak
        if (isw <= 0) then
          rhs(3:15:4, 1) = -1.d3
          rhs(4:16:4, 1) =  0.d0
        endif
c
        return

      end subroutine uel
