c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c 
c This program computes the observational quantities for a KBO observer
c from platform described by a trajectory file.
c
c Usage: EncounterPredict orbit_file.txt trajectory_file.txt
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      
      include 'SurveySubs.f'

      program EP

      implicit none

      integer*4 screen, keybd, verbose
      integer*4 code, lun_trajectory, lun_orbit
      integer*4 nw_max, ierr, stderr
      real*8 mt

      parameter
     $     (screen = 6, keybd = 5, verbose = 9,
     $     nw_max = 30, stderr = 0)
      character buffer*800

      real*8 Pi, TwoPi, drad

      parameter (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi,
     $  drad = Pi/180.0d0)


c     Barycentric Elements
      real*8 a, e, inc, node, peri, M, epoch 
c     Physical Parameters 
      real*8 h, gb
C     circumstances
      real*8 obs_jday, nh_pos(3), nh_vel(3), ros
c     Lorri Observer Computed values
      real*8 ra, dec, r, delta, alpha, mag, pos(3)

      parameter (lun_trajectory = 12, lun_orbit = 11, code = 500, 
     $  gb = 0.15)
C     inputs
      character orbit_file*200, trajectory_file*200

      CALL getarg(1, orbit_file)
      orbit_file = trim(orbit_file)
      CALL getarg(2, trajectory_file)
      trajectory_file = trim(trajectory_file)

C     Orbit file should be formated. 
C     a e inc node peri M epoch H
      open (unit=lun_orbit, file=orbit_file, status='OLD')
      read (lun_orbit, *, err=9999)
     $ a, e, inc, node, peri, M, epoch, h
      inc = inc * drad
      node = node * drad
      peri = peri * drad
      M = M * drad
      close (lun_orbit)

      open (unit=lun_trajectory, file=trajectory_file, status='OLD')
c     Skip First line (should be header line)
      read(lun_trajectory, *) buffer

c     First possible date for observation
      obs_jday = -1.0

      ierr = 0
      do while ( ierr < 1 )
C        Get position of NH at obs_jday
         call ObsPos(-lun_trajectory, obs_jday, 
     $        nh_pos, nh_vel, ros, ierr)
      if (ierr .ne. 0) then
         write (stderr, *) 'Error reading trajectory '//trajectory_file
         call exit(255)
      end if

C        Get the RA/DEC at date
         mt = M
     $        + (TwoPi/(a**1.5d0*365.25d0))*(obs_jday-epoch)
         mt = mt - int(mt/TwoPi)*TwoPi
         call pos_cart(a, e, inc, node, peri, mt, pos(1),
     $        pos(2), pos(3))
         call RADECeclXV(pos, nh_pos, delta, ra, dec)
C        Compute the magnitude
         r = sqrt(pos(1)**2 + pos(2)**2 + pos(3)**2)
         call AppMag(r, delta, ros, h, gb, alpha, 
     $        mag, ierr) 
C        Spacecraft observers in open filter
         mag = mag + 0.7
         if ( mag < 21 ) then
             write (screen, '(6(f16.5,1x))') obs_jday,
     $            ra/drad, dec/drad, delta, mag,  alpha/drad
         end if
      end do

 9999 CONTINUE
      end program EP


