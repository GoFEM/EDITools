! Interface to call mtcomp code from C++

subroutine call_mtcomp(spectra, frequency, avgt, &
                       impedance, tipper,  &
                       apparent_resistivity, phase, &
                       impedance_err, tipper_err, &
                       appres_err, phase_err) bind(c)
    use iso_c_binding
    implicit none

    DOUBLE COMPLEX Z(3,2)
    REAL*8 rho(2,2), phs(2,2), srh(2,2), sph(2,2)
    REAL*8 vimp(3,2), delta, avgs, angle/0.0/
    REAL*8 temp(7,7), freq
    INTEGER i,j,idx

    REAL (c_double), dimension(49), intent(in) :: spectra
    REAL (c_double), intent(in), value :: frequency, avgt
    COMPLEX (c_double_complex), dimension(4), intent(out) :: impedance
    COMPLEX (c_double_complex), dimension(2), intent(out) :: tipper
    REAL (c_double), dimension(4), intent(out) :: apparent_resistivity
    REAL (c_double), dimension(4), intent(out) :: phase
    REAL (c_double), dimension(4), intent(out) :: impedance_err
    REAL (c_double), dimension(2), intent(out) :: tipper_err
    REAL (c_double), dimension(4), intent(out) :: appres_err
    REAL (c_double), dimension(4), intent(out) :: phase_err

    freq = frequency
    avgs = avgt

    idx = 1;
    do i=1,7
      do j=1,7
        temp(i,j) = spectra(idx)
        idx = idx + 1;
      enddo
    enddo

    call mtcomp(temp,freq,delta,avgs,Z,vimp,rho,phs,srh,sph)

    impedance(1) = Z(1,1);
    impedance(2) = Z(1,2);
    impedance(3) = Z(2,1);
    impedance(4) = Z(2,2);

    tipper(1) = Z(3,1);
    tipper(2) = Z(3,2);

    impedance_err(1) = dsqrt(vimp(1,1));
    impedance_err(2) = dsqrt(vimp(1,2));
    impedance_err(3) = dsqrt(vimp(2,1));
    impedance_err(4) = dsqrt(vimp(2,2));

    tipper_err(1) = dsqrt(vimp(3,1));
    tipper_err(2) = dsqrt(vimp(3,2));

    apparent_resistivity(1) = rho(1,1);
    apparent_resistivity(2) = rho(1,2);
    apparent_resistivity(3) = rho(2,1);
    apparent_resistivity(4) = rho(2,2);

    phase(1) = phs(1,1);
    phase(2) = phs(1,2);
    phase(3) = phs(2,1);
    phase(4) = phs(2,2);

    appres_err(1) = srh(1,1);
    appres_err(2) = srh(1,2);
    appres_err(3) = srh(2,1);
    appres_err(4) = srh(2,2);

    phase_err(1) = sph(1,1);
    phase_err(2) = sph(1,2);
    phase_err(3) = sph(2,1);
    phase_err(4) = sph(2,2);
end subroutine call_mtcomp
