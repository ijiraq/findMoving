c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c 
c This program compute the position of objects at different times.
c Should be used once basic detectabilty is determined, ie this uses the 
c output of NHDetectability.
c
c This code add 0.7 to the lorri_mag values as Detectability didn't 
c covert to V from H_r. 
c
c Usage: Obj_at_Date model trajectory obs-date
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      
      include 'SurveySubs.f'

      program EP

      implicit none

      integer*4 n_obj_max, screen, keybd, verbose, 
      integer*4 code, lun_trajectory, line_number
      integer*4 nw_max, nw, ierr, nchar

      parameter
     $     (n_obj_max = 10000, screen = 6, keybd = 5, verbose = 9,
     $     nw_max = 30)
      character line*200, buffer*800,  word(nw_max)*80
      integer*4 lw(nw_max)

      real*8 Pi, TwoPi, drad

      parameter (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi,
     $  drad = Pi/180.0d0)


c     Barycentric Elements
      real*8 a, e, inc, node, peri, M, epoch, 
c     Physical Parameters 
      real*8 h, gb
C     circumstances
      real*8 obs_jday, obspos(3), ros, pos
c     Lorri Observer Computed values
      real*8 ra, dec, r, delta, ph, mag, pos(3)

      parameter (lun_trajectory = 12, code = 500, gb = 0.15)

C     inputs
      character orbit_file*200, trajectory_file*200
     $     word(nw_max)*80

      CALL getarg(1, orbit_file)
      orbit_file = trim(orbit_file)
      CALL getarg(2, trajectory_file)
      trajectory_file = trim(trajectory_file)

C     Orbit file should be formated. 
C     a e inc node peri M epoch
      open (unit=lun_h, file=orbit_file, status='OLD')
      read (lun_h, *, err=9999)  a, e, inc, node, peri, M, epoch
      close (lun_h)

      open (unit=lun_trajectory, file=trajectory_file, status='OLD')
c     Skip First line (should be header line)
      read(lun_trajectory, *) buffer
      if (ierr .ne. 0) then
         write *, 'Error reading trajectory '//trajectory_file
         exit(255)
      end if

c     First possible date for observation
      obs_jday = -1.0

      do while ( rierr < 1 )
c     Get position of NH at obs_jday
         call ObsPos(-lun_trajectory, obs_jday, 
     $        obspos, tmp, ros, ierr)

c     Get the RA/DEC at date
         mt = M
     $        + (TwoPi/(a**1.5d0*365.25d0))*(obs_jday-epoch)
         mt = mt - int(mt/TwoPi)*TwoPi
         call pos_cart(a, e, inc, node, peri, mt, pos(1),
     $        pos(2), pos(3))
         call RADECeclXV(pos, obspos, delta, ra, dec)
c     Compute the magnitude
         r = sqrt(pos(1)**2 + pos(2)**2 + pos(3)**2)
         call AppMag(r, delta, ros, h, gb, alpha, 
     $        mag, ierr) 

         lorri_mag = lorri_mag + 0.7

         write (lun_h, '(6(f16.5,1x))', advance='no') obs_jday, 
     $        ra, dec, delta, mag,  alpha

 8000 end do
      close (lun_h)

 9999 CONTINUE
      end program EP


