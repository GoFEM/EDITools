c-- This file is from the edi2j program by Alan G Jones,
c-- all copyright belongs to the original author.
c-- The code was downloaded from: https://www.mtnet.info/main/source.html

c----------------------------------------------------------------------

c     SUBROUTINE MTCOMP

c---------------------------------------------------------------------

        SUBROUTINE mtcomp(xpr,frq,angle,avgs,Z,vimp,rho,phs,srh,sph)

        INTEGER Ex,Ey,Hx,Hy,Hz,Rx,Ry
        DOUBLE COMPLEX Zxx,Zxy,Zyx,Zyy,Zzx,Zzy,Temp1
        DOUBLE COMPLEX Z(3,2), cavg(7,7), swap
        DOUBLE COMPLEX Num,Den
        REAL*8 Nx2,Ny2,Nz2 ,Nphi, Ax2, Ay2, T2, Den2, D3
        REAL*8 rho(2,2), phs(2,2), srh(2,2), sph(2,2)
        REAL*8 vimp(3,2), vtr(2), angle, theta, Pi, Rtd, Rtd2
        REAL*8 ct, st, s2, c2, rotmat(2,2), avgs, frq
        REAL*8 xpr(7,7)


c---------------------------------------------------------------------

c-- set some parameters

        Pi = 3.141592
        Rtd = 180./Pi
        Rtd2 = Rtd**2

c-- set spectra matrix order (Jones order)

        Ex = 3
        Ey = 4
        Hx = 1
        Hy = 2
        Hz = 5
        Rx = 6
        Ry = 7

c-- unpack spectral matrix

        do i = 1, 7
          do j = i, 7
            if (i .eq. j) then
              cavg(i,j) = cmplx(xpr(i,j),0.0)
            else
              cavg(i,j) = (cmplx(xpr(j,i),-xpr(i,j)))
              cavg(j,i) = (cmplx(xpr(j,i),xpr(i,j)))
            end if
          end do
        end do

c-- change matrix from geotools to jones (geotools Hx Hy Hz Ex Ey Rx Ry)

c swap columns 4 & 5 (Ex and Ey)

        do i = 1, 7
          swap = cavg(i,4)
          cavg(i,4) = cavg(i,5)
          cavg(i,5) = swap
        enddo

c swap columns 3 & 5 (Hz and Ex)

        do i = 1, 7
          swap = cavg(i,3)
          cavg(i,3) = cavg(i,5)
          cavg(i,5) = swap
        enddo

c swap rows 4 & 5 (Ex and Ey)

        do i = 1, 7
          swap = cavg(4,i)
          cavg(4,i) = cavg(5,i)
          cavg(5,i) = swap
        enddo

c swap rows 3 & 5 (Hz and Ex)

        do i = 1, 7
          swap = cavg(3,i)
          cavg(3,i) = cavg(5,i)
          cavg(5,i) = swap
        enddo

c-- rotate spectral matrix

        IF (ABS(angle) .GT. 0.01) THEN
         call rotss7(angle,cavg)
        ENDIF

c-- calculate IMPEDANCE

        Den = Cavg(Hx,Rx)*Cavg(Hy,Ry)-Cavg(Hx,Ry)*Cavg(Hy,Rx)
        Num = Cavg(Ex,Rx)*Cavg(Hy,Ry)-Cavg(Hy,Rx)*Cavg(Ex,Ry)
        Zxx = Num/Den
        Num = Cavg(Hx,Rx)*Cavg(Ex,Ry)-Cavg(Ex,Rx)*Cavg(Hx,Ry)
        Zxy = Num/Den
        Num = Cavg(Ey,Rx)*Cavg(Hy,Ry)-Cavg(Hy,Rx)*Cavg(Ey,Ry)
        Zyx = Num/Den
        Num = Cavg(Hx,Rx)*Cavg(Ey,Ry)-Cavg(Ey,Rx)*Cavg(Hx,Ry)
        Zyy = Num/Den

c-- calculate TIPPER IMPEDANCE

        IF (Hz .NE. 0) THEN
          Den = Cavg(Hx,Rx)*Cavg(Hy,Ry)-Cavg(Hx,Ry)*Cavg(Hy,Rx)
          Num = Cavg(Hz,Rx)*Cavg(Hy,Ry)-Cavg(Hy,Rx)*Cavg(Hz,Ry)
          Zzx = Num/Den
          Num = Cavg(Hx,Rx)*Cavg(Hz,Ry)-Cavg(Hz,Rx)*Cavg(Hx,Ry)
          Zzy = Num/Den
        END IF

c-- save impedances and transfer fuctions

        Z(1,1) = Zxx
        Z(1,2) = Zxy
        Z(2,1) = Zyx
        Z(2,2) = Zyy
        Z(3,1) = Zzx
        Z(3,2) = Zzy

c--  calculate RHO-AMP & RHO-PHASE
c (Note that spectra assumed to be in field units!!!)

        rho(1,1) = 0.2*ABS(Zxx)**2/Frq
        rho(1,2) = 0.2*ABS(Zxy)**2/Frq
        rho(2,1) = 0.2*ABS(Zyx)**2/Frq
        rho(2,2) = 0.2*ABS(Zyy)**2/Frq
        phs(1,1) = ATAN2(imag(Zxx),real(Zxx))*Rtd
        phs(1,2) = ATAN2(imag(Zxy),real(Zxy))*Rtd
        phs(2,1) = ATAN2(imag(Zyx),real(Zyx))*Rtd
        phs(2,2) = ATAN2(imag(Zyy),real(Zyy))*Rtd

c-- calculate VARIANCE

c    JOHN STODT VARIANCE EXPRESSIONS FOR REMOTE REFERENCE.

        Nx2 =   Cavg(Ex,Ex)/Cavg(Rx,Rx)
     &          -2*REAL( Zxx*Cavg(Hx,Ex)+Zxy*Cavg(Hy,Ex)
     &          -Zxx*CONJG(Zxy)*Cavg(Hx,Hy))/Cavg(Rx,Rx)
     &          +Zxx*CONJG(Zxx)*Cavg(Hx,Hx)/Cavg(Rx,Rx)
     &          +Zxy*CONJG(Zxy)*Cavg(Hy,Hy)/Cavg(Rx,Rx)
        Nx2 = DABS(Nx2)
        Ny2 =   Cavg(Ey,Ey)/Cavg(Rx,Rx)
     &          -2*REAL( Zyx*Cavg(Hx,Ey)+Zyy*Cavg(Hy,Ey)
     &          -Zyx*CONJG(Zyy)*Cavg(Hx,Hy))/Cavg(Rx,Rx)
     &          +Zyx*CONJG(Zyx)*Cavg(Hx,Hx)/Cavg(Rx,Rx)
     &          +Zyy*CONJG(Zyy)*Cavg(Hy,Hy)/Cavg(Rx,Rx)
        Ny2 = DABS(Ny2)
        IF (Hz .NE. 0) THEN
           T2 = Cavg(Hz,Hz)/Cavg(Rx,Rx)
     &          -2*REAL( Zzx*CONJG(Cavg(Hz,Hx))+Zzy*CONJG(Cavg(Hz,Hy))
     &          -Zzx*CONJG(Zzy)*Cavg(Hx,Hy))/Cavg(Rx,Rx)
     &          +Zzx*CONJG(Zzx)*Cavg(Hx,Hx)/Cavg(Rx,Rx)
     &          +Zzy*CONJG(Zzy)*Cavg(Hy,Hy)/Cavg(Rx,Rx)
           T2 = DABS(T2)
        END IF

        Ax2=    Cavg(Hy,Ry)*CONJG(Cavg(Hy,Ry))
     &      +Cavg(Ry,Ry)/Cavg(Rx,Rx)*Cavg(Hy,Rx)*CONJG(Cavg(Hy,Rx))
     &      -2*REAL(Cavg(Rx,Ry)/Cavg(Rx,Rx)*
     &       CONJG(Cavg(Hy,Ry))*Cavg(Hy,Rx))
        Ax2 = DABS(Ax2)
        Ay2 =   Cavg(Ry,Ry)/Cavg(Rx,Rx)
     &          *Cavg(Hx,Rx)*CONJG(Cavg(Hx,Rx))
     &      +   Cavg(Hx,Ry)*CONJG(Cavg(Hx,Ry))
     &      -2*REAL(Cavg(Rx,Ry)/Cavg(Rx,Rx)
     &            *CONJG(Cavg(Hx,Ry))*Cavg(Hx,Rx))
        Ay2 = DABS(Ay2)
        Temp1=(Cavg(Hx,Rx)*Cavg(Hy,Ry)
     &         -Cavg(Hx,Ry)*Cavg(Hy,Rx))/Cavg(Rx,Rx)
        Den2 = Temp1 * CONJG(Temp1)
        D3 = avgs * Den2

        vimp(1,1) = (Nx2*Ax2/D3)
        vimp(1,2) = (Nx2*Ay2/D3)
        vimp(2,1) = (Ny2*Ax2/D3)
        vimp(2,2) = (Ny2*Ay2/D3)
        IF (Hz .NE. 0) THEN
          vimp(3,1) = (T2*Ax2/D3)
          vimp(3,2) = (T2*Ay2/D3)
        END IF

        srh(1,1) = dsqrt((0.3772*vimp(1,1))/ABS(Zxx)**2.0)
        srh(1,2) = dsqrt((0.3772*vimp(1,2))/ABS(Zxy)**2.0)
        srh(2,1) = dsqrt((0.3772*vimp(2,1))/ABS(Zyx)**2.0)
        srh(2,2) = dsqrt((0.3772*vimp(2,2))/ABS(Zyy)**2.0)


c       sph(1,1) = sqrt(8.*COS(phs(1,1)/Rtd)**4
c    &                  *vimp(1,1)*CABS(Zxx)**2/( 2*REAL(Zxx) )**4)
c       sph(1,2) = sqrt(8.*COS(phs(1,2)/Rtd)**4
c    &                  *vimp(1,2)*CABS(Zxy)**2/( 2*REAL(Zxy) )**4)
c       sph(2,1) = sqrt(8.*COS(phs(2,1)/Rtd)**4
c    &                  *vimp(2,1)*CABS(Zyx)**2/( 2*REAL(Zyx) )**4)
c       sph(2,2) = sqrt(8.*COS(phs(2,2)/Rtd)**4
c    &                  *vimp(2,2)*CABS(Zyy)**2/( 2*REAL(Zyy) )**4)

        sph(1,1) = dsqrt((rtd**2.0*0.5*vimp(1,1))/
     &                  (ABS(Zxx)**2.0) )
        sph(1,2) = dsqrt((rtd**2.0*0.5*vimp(1,2))/
     &                  (ABS(Zxy)**2.0) )
        sph(2,1) = dsqrt((rtd**2.0*0.5*vimp(2,1))/
     &                  (ABS(Zyx)**2.0) )
        sph(2,2) = dsqrt((rtd**2.0*0.5*vimp(2,2))/
     &                  (ABS(Zyy)**2.0) )

        RETURN
        END

c----------------------------------------------------------------
      SUBROUTINE rotss7( angle, ss )
c----------------------------------------------------------------
      integer j, i
      REAL*8    theta, a, b, c, d, e, angle
      DOUBLE COMPLEX ss(7,7), sx(7,7)
C
C     SUBROUTINE ROTATES THE POWER MATRIX  SS  BY ANGLE  'ANGLE'
C
C     ALL VECTORS ROTATED, including Hz
c


      THETA=ANGLE*3.141592/180.
      D = COS(THETA)
      E = SIN(THETA)
      A = D*D
      B = E*E
      C = D*E
C
      sx(1,1)=ss(1,1)*A +ss(2,2)*B +(ss(1,2)+ss(2,1))*C
      sx(1,2)=ss(1,2)*A -ss(2,1)*B +(ss(2,2)-ss(1,1))*C
      sx(1,3)=ss(1,3)*A +ss(2,4)*B +(ss(1,4)+ss(2,3))*C
      sx(1,4)=ss(1,4)*A -ss(2,3)*B +(ss(2,4)-ss(1,3))*C
      sx(1,6)=ss(1,6)*A +ss(2,7)*B +(ss(1,7)+ss(2,6))*C
      sx(1,7)=ss(1,7)*A -ss(2,6)*B +(ss(2,7)-ss(1,6))*C
      sx(2,2)=ss(2,2)*A +ss(1,1)*B -(ss(1,2)+ss(2,1))*C
      sx(2,3)=ss(2,3)*A -ss(1,4)*B +(ss(2,4)-ss(1,3))*C
      sx(2,4)=ss(2,4)*A +ss(1,3)*B -(ss(2,3)+ss(1,4))*C
      sx(2,6)=ss(2,6)*A -ss(1,7)*B +(ss(2,7)-ss(1,6))*C
      sx(2,7)=ss(2,7)*A +ss(1,6)*B -(ss(2,6)+ss(1,7))*C
      sx(3,3)=ss(3,3)*A +ss(4,4)*B +(ss(3,4)+ss(4,3))*C
      sx(3,4)=ss(3,4)*A -ss(4,3)*B +(ss(4,4)-ss(3,3))*C
      sx(3,6)=ss(3,6)*A +ss(4,7)*B +(ss(3,7)+ss(4,6))*C
      sx(3,7)=ss(3,7)*A -ss(4,6)*B +(ss(4,7)-ss(3,6))*C
      sx(4,4)=ss(4,4)*A +ss(3,3)*B -(ss(3,4)+ss(4,3))*C
      sx(4,6)=ss(4,6)*A -ss(3,7)*B +(ss(4,7)-ss(3,6))*C
      sx(4,7)=ss(4,7)*A +ss(3,6)*B -(ss(4,6)+ss(3,7))*C
      sx(6,6)=ss(6,6)*A +ss(7,7)*B +(ss(6,7)+ss(7,6))*C
      sx(6,7)=ss(6,7)*A -ss(7,6)*B +(ss(7,7)-ss(6,6))*C
      sx(7,7)=ss(7,7)*A +ss(6,6)*B -(ss(6,7)+ss(7,6))*C
c
c rotate xpower with Hz just once
c
      sx(5,5) = ss(5,5)
      sx(1,5) = ss(1,5)*d + ss(2,5)*e
      sx(2,5) =-ss(1,5)*e + ss(2,5)*d
      sx(3,5) = ss(3,5)*d + ss(4,5)*e
      sx(4,5) =-ss(3,5)*e + ss(4,5)*d
      sx(5,6) = ss(5,6)*d + ss(5,7)*e
      sx(5,7) =-ss(5,6)*e + ss(5,7)*d

      sx(2,1)=CONJG(sx(1,2))
      sx(3,1)=CONJG(sx(1,3))
      sx(3,2)=CONJG(sx(2,3))
      sx(4,1)=CONJG(sx(1,4))
      sx(4,2)=CONJG(sx(2,4))
      sx(4,3)=CONJG(sx(3,4))
      sx(5,1)=CONJG(sx(1,5))
      sx(5,2)=CONJG(sx(2,5))
      sx(5,3)=CONJG(sx(3,5))
      sx(5,4)=CONJG(sx(4,5))
      sx(6,1)=CONJG(sx(1,6))
      sx(6,2)=CONJG(sx(2,6))
      sx(6,3)=CONJG(sx(3,6))
      sx(6,4)=CONJG(sx(4,6))
      sx(6,5)=CONJG(sx(5,6))
      sx(7,1)=CONJG(sx(1,7))
      sx(7,2)=CONJG(sx(2,7))
      sx(7,3)=CONJG(sx(3,7))
      sx(7,4)=CONJG(sx(4,7))
      sx(7,5)=CONJG(sx(5,7))
      sx(7,6)=CONJG(sx(6,7))

      do j = 1, 7
        do i = 1, 7
          ss(i,j) = sx(i,j)
        enddo
      enddo

      RETURN
      END
