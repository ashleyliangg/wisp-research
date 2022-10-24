c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
      program readdata
c
c This program reads in 2D astrosphere data in the format of the
c f4restp.dat output files of the simulations, and prepares the
c data into several 2D variable arrays that can be used for calculations.
c This version is supposed to calculate the hydrogen Ly-alpha absorption.
c
c Input data filename and the 2 output files for dtau and absorb
c are given in the command line. If f4restp.dat contains several 
c "Zones" (=time snapshots), then create
c another file with only one snapshot with the unix command:
c tail -16423 f4restp.dat > [inputfilename.dat]
c (one might have to fiddle with the number 16423 to make "Variables ="
c  the first line of the shortened file.) OR with the command:
c sed -n -e '/20000/,$p' f4restp.dat > [inputfilename.dat]
c (adjust the 20000 to the desired timezone)
c 
c Input: 
c Output: two data files - [dtaufilename] and [absorbfilename] where the dtau
c contains the tau for each angle and radial distance, the ISM tau, and the difference
c and the absorb contains the sum of all dtaus for each angle, the sum for the ISM
c (interstellar medium), and the difference (star's dtau - ISM's dtau)
c
c Authors: Ashley Liang and Hans Mueller
c
c last update 2022 October 18
c  

      implicit none

c ***********************
c Start of Variable Definitions
c ***********************

c parameters
      integer is, ie, js, je,jep1, nheader, fn    !i-start (lowest radial index), i-end, j-start (lowest angular index), j-end, j-end + 1
      parameter(is=3, js=3, ie=434, je=39)     !1 and 2 are boundary zone so never used
      parameter(jep1=je+1)                      !1 extra point for interpolations
      parameter(nheader=6)  !header lines to skip

      real*8 pi18, pi, t_ism, d_ism, v_ism
      parameter(t_ism=7500) !all ism values from inzeus parameters
      parameter(d_ism=0.18)
      parameter(v_ism=27.8)
      parameter (pi=3.14159265358979D0)
      parameter (pi18=pi/180.0D0)
      real*8 baseline, c, f0, massH, kb, df, gh
      parameter(baseline=1.05)
      parameter(c=299792458)  !speed of light in m/s
      parameter(f0=c/121.567D-9) !lyman-alpha frequency
      parameter(kb=1.380649D-23) !J/K
      parameter(massH=1.6738D-27) !mass of hydrogen in kg
      parameter(fn=200)
      parameter(df = 2.0d12/fn) !derived quantity for equidistant frequencies (Hz)
      parameter(gh = 1.02D7)

c local simple variables
      integer i,j, region, f
      real*8  x, y, th, cs, sn
     %      , n0, vx0, vy0, t0
     %      , n1, vx1, vy1, t1
     %      , n2, vx2, vy2, t2
     %      , n3, vx3, vy3, t3
      logical qr           !logical flag
      character*72 sfn, line     !file name, line

      real*8 r(ie), coldens(jep1)
     % , d(ie,jep1), t(ie,jep1)
     % , vx(ie,jep1), vy(ie,jep1), vr(ie,jep1), vt(ie,jep1)
     % , dtau(fn, ie, jep1), xf(fn), sum_tau(jep1, fn)
     % , sum_tau_ism(jep1, fn), phi_r(ie,jep1,fn)
     % , dtau_ism(fn, ie, jep1)

      logical begin, err
      real*8 base_dens, integral, total, vth, dopp_w, phi_real, phi_im, 
     & fd, base_tau
      character(len=50) filename, densfile, dtaufile, absorbfile


c ***********************
c Start of the Code
c ***********************
c  Checks number of arguments, assigns arguments to corresponding variable names
      if (command_argument_count() .ne. 3) then
        print *, 'wrong number of args: ',command_argument_count()
        STOP 1
      endif
      call get_command_argument(1, filename)
      call get_command_argument(2, dtaufile)
      call get_command_argument(3, absorbfile)

c  open the input data file
      open(15,file=filename,status='unknown')    !open file
c      do i=1,nheader       !read past the 6 header lines
c        read(15,*) line
c      enddo
      read(15,*) line      !read past header line = 1
      print *,'last read header line:',line

c  main read-in loop 2D
      do j=js,jep1    !loop over angles
        qr = (j.eq.10)
        do i=is,ie     !loop over radii
          read(15,*) x, y
     %      , n0, vx0, vy0, t0, region
     %      , n1, vx1, vy1, t1
     %      , n2, vx2, vy2, t2
     %      , n3, vx3, vy3, t3
           d(i,j) = n1   !neutral fluid 1 density
          vx(i,j) = vx1  !x-componenet of velocity
          vy(i,j) = vy1  !y-component
           t(i,j) = t1   !temperature of kelvin
          if (qr) r(i) = sqrt(x*x+y*y)  !index to radial units, r is array of radial distances
        enddo !i   !radial distance
      enddo   !j   !36 angles (roughly 0 to roughly 180, 5 degree resolution)
      close(15)                               !close file
      print *,'last read position fluid 1:'
      print *, x, y, r(is),r(ie) !r(i) is the actual radius and i is just number
      print *,'read-in complete'


c  another 2D loop for derived variables
      do j=js,jep1    !loop over angles
        th = ((j-js)*5.0d0-2.5d0)*pi18  !angle in radians
        cs = cos(th)
        sn = sin(th)
        do i=is,ie     !loop over radii
          vr(i,j) =  cs*vx(i,j)+sn*vy(i,j)  !radial velocity
          vt(i,j) = -sn*vx(i,j)+cs*vy(i,j)  !tangential velocity
        enddo !i
      enddo   !j

      print *,'sample velocity outer edge, 42.5 degrees:'
      th=(12-js)*5.0d0-2.5d0
      print *, r(ie), th, r(ie)*cos(th*pi18), r(ie)*sin(th*pi18)
      print *, vr(ie,12), vt(ie,12),vx(ie,12),vy(ie,12)

c ====================================================================
c  at this point, there are available six 2D arrays describing inter-
c  stellar neutral hydrogen:
c  d (=density, cm^-3), vx(x-component velocity, km/s), vy(vel. y component),
c  vr(radial velocity), vt(tangential velocity), t(temperature, K).
c  1D array describing the grid:
c    r(i): r location of a grid point with integer r-index i
c    th=(j-js)*5.0d0-2.5d0: angle theta of a grid point with integer index j
c    x = r(i)*cos(th)
c    y = r(i)*sin(th)
c ====================================================================

c !!!NO LONGER USING THE COMMENTED OUT CODE BELOW BECAUSE BELOW DOESN'T TAKE FREQUENCY INTO ACCOUNT!!!
c CALCULATE COLUMN DENSITY WITHOUT CONSIDERING FREQUENCY
c       base_dens = d(ie,js)

c       do j=js,jep1 !loop over all angles
c c       first find the integral of last one
c         begin = .false.
c         total = 0.0
c         do i=ie - 1,is + 1,-1  ! backwards to smallest radii + 1
c           if (begin .eqv. .false.) then

c c           if integral is 5% larger than end, add to the total and set begin=true
c             if (d(i,j) .gt. (baseline * base_dens)) then
c c           find integral
c               integral = (d(i,j) - base_dens) * (r(i) - r(i - 1))
c               total = total + integral
c               begin = .true.
c               print *,"radial distance at ",j," is ",r(i)
c             endif
c           else
c c           find integral and add to total
c             integral = (d(i,j) - base_dens) * (r(i) - r(i - 1))
c             total = total + integral
c           endif
c         enddo
c c       add total to array of totals
c         coldens(j) = total
c       enddo

c ====================================================================
c BELOW CALCULATES THE FREQUENCY-DEPENDENT TAUS (how much light
c goes through the section) AS WELL AS THE OVERALL TAU RADIALLY
c (every 5 degrees)
c ====================================================================

c Calculate the frequency-dependent absorption (a_f) for each degree and radial distance
      
c     section out frequency graph centered on f0
      do f=1, fn
        xf(f) = (f0 - (df*fn/2)) + df*f 
      enddo
      
      do j=js, jep1
        do i=is, ie
c         get doppler shifted frequency
          fd = (1-(vr(i,j) * 1000/c))*f0  !vr is in km, so convert that to meters
c         calculate thermal velocity
          vth = sqrt(2*kb*t(i,j)/massH)    !get thermal velocity in m/s
c         calculate doppler width
          dopp_w = (vth/c)*f0
          do f = 1, fn  !x is just an index num (i), fn is just an arbitrary number splitting the frequency bell curve
c           calculate phi (real part) by calling WOFZ.f
            call WOFZ((xf(f) - fd)/dopp_w, gh/(4 * pi * dopp_w), 
     &      phi_real, phi_im, err)
            phi_r(i,j,f) = phi_real
            if (phi_real < 0) then
              print *,phi_real," phi_real negative at ",i,j,f !seems to be negative starting with i=240, j=17 (angle), f=1 all the way to i=434 j=40, f=25. radial distance negative starting at diff points each angle
            endif

c           calculate delta tau (aka frequency-dependent absorption)
            dtau(f,i,j) = (0.4162D6 * 1.496D11) *1.602D-19*   !fosc - conversion of density cm^3 to m^3 --> 10^6, conversion of AU to meters
     &      1.602D-19*d(i,j)*phi_real/ !phi() function
     &      (4.0D0*c*9.11D-31*8.85D-12*dopp_w*1.7725)
            if (dtau(f,i,j) < 0) then  !dtau is always negative when phi_real is negative
              print *,"d_tau negative at ",i,j,f, " phi_real of ",
     &        phi_real
            endif


          enddo
        enddo
      enddo


c Getting summation of frequency-dependent absorption for each absorption section of curve

      do j=js, jep1  !go through each angle
        do f=1, fn
c         first find the integral of last one
          total = 0.0
          do i=is, ie - 1 ! inner radii to outer radii
c           find integral and add to total
            integral = dtau(f,i,j) * (r(i + 1) - r(i))
            total = total + integral
          enddo
          sum_tau(j,f) = total
        enddo
      enddo


c Calculate frequency tau for ISM
      do j=js,jep1    !loop over angles
        th = ((j-js)*5.0d0-2.5d0)*pi18  !angle in radians
        do i=is,ie     !loop over radii
          vr(i,j) =  -cos(th)* v_ism !radial velocity, vy is 0 because along x axis, negative because coming from right to left
        enddo !i
      enddo   !j

      do j=js, jep1
        do i=is, ie
c         get doppler shifted frequency
          fd = (1-(vr(i,j) * 1000/c))*f0 !vr is in km, convert to m
c         calculate thermal velocity in m/s
          vth = sqrt(2*kb*t_ism/massH)
c         calculate doppler width
          dopp_w = (vth/c)*f0 
          do f = 1, fn  !x is just an index num (i), fn is just an arbitrary number splitting the frequency bell curve
c           calculate phi (real part) by calling WOFZ.f
            call WOFZ((xf(f) - fd)/dopp_w, gh/(4 * pi * dopp_w), 
     &      phi_real, phi_im, err)

c         calculate delta tau (aka frequency-dependent absorption)
          dtau_ism(f,i,j) = (0.4162D6 * 1.496D11) *1.602D-19*   !fosc - conversion of density cm^3 to m^3 --> 10^6, conversion of AU to meters
     &    1.602D-19*d_ism*phi_real/ !phi() function
     &    (4.0D0*c*9.11D-31*8.85D-12*dopp_w*1.7725)
          enddo
        enddo
      enddo

      do j=js, jep1  !go through each angle
        do f=1, fn
c         first find the integral of last one
          total = 0.0
          do i=is, ie - 1 ! inner radii to outer radii
c           find integral and add to total
            integral = dtau_ism(f,i,j) * (r(i + 1) - r(i))
            total = total + integral
          enddo
          sum_tau_ism(j,f) = total
        enddo
      enddo











c ====================================================================
c Lastly, write out the results

c      open(21,file='coldens.dat',status='unknown')
c     write(21,  * ) 'VARIABLES = theta coldens'
c      write(21, 701) jep1-js+1
c      do j=js,jep1
c        th=(j-js)*5.0d0-2.5d0   ! double float
c        
c        write(21,*) th,coldens(j)
c      enddo
c      close(21)

c     print out d_tau for simulation data
      
      open(22, file=dtaufile, status='unknown')    !dtau file
      write(22, * ) 'VARIABLES = frequency,radius,angle,dtau,'//
     &'dtau_ism, diff_tau'
      write(22, 702) fn,ie-is+1,jep1-js+1
      do j=js,jep1
        th=(j-js)*5.0d0-2.5d0   ! double float
        do i=is,ie
          do f=1,fn
            write(22, 713) xf(f),r(i),th,dtau(f,i,j),dtau_ism(f,i,j),
     &      dtau(f,i,j) - dtau_ism(f,i,j)
          enddo
        enddo
      enddo
      close (22)

      open(22, file=absorbfile, status='unknown')  !absorb file
      write(22, * ) 'VARIABLES = frequency,angle,tau, tau_ism, diff'
      write(22, 704) fn, jep1-js+1
      do j=js,jep1
        th=(j-js)*5.0d0-2.5d0   ! double float
          do f=1,fn
            write(22, 714) xf(f),th, sum_tau(j, f), sum_tau_ism(j, f),
     &      sum_tau(j, f) - sum_tau_ism(j, f)
          enddo
      enddo



701   format('ZONE I=',i5,', F=POINT')   !integer saving 5 spaces
702   format('ZONE I=',i4,', J=',i4,', K=',i4)
703   format('ZONE I=',i4,', J=',i4,', K=',i4,', F=',i4)
704   format('ZONE I=',i4,', J=',i4)
712   format(3f10.3,1e11.4,2e11.4,4e12.4,1e11.4)  !add e12.4 for 3D
713   format(6(1pe12.4))
714   format(1pe14.6,4(1pe12.4))
715   format(i5,1pe12.4)
716   format(1pe14.6,6(1pe12.4))

      end
c end of main program (readdata)

c**********************************************************************
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
