c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      
      include 'SurveySubs.f'

      program XX

      implicit none

      integer*4 n_obj_max, screen, keybd, verbose, lun_h, lun_t
      integer*4 code

      parameter
     $  (n_obj_max = 10000, screen = 6, keybd = 5, verbose = 9)

      real*8 Pi, TwoPi, drad

      parameter (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi,
     $  drad = Pi/180.0d0)

c color array NEEDS to be length 10 or more!
      real*8 a, e, inc, node, peri, M, h, epoch, m_int, d_ra, d_dec,
     $  r, delta, ra, dec, random, mt, extra(10), gb, ph, period, amp,
     $  jday, m_rand, eff, rn_iter, eff_lim, h_rand, pos(3), obpos(3),
     $  ros, tmp(3), obs_jday

      integer*4 n_hits, n_track, ierr, seed, flag, isur, ic, n_iter,
     $  n_track_max, nchar, rierr, co, j, ext_index

      character distri_file*80, survey_dir*100, line*200,
     $  trk_outfile*80, det_outfile*80,
     $  comments*100, surna*10

      logical keep_going

      lun_h = 10
      lun_t = 11
      keep_going = .true.
      code = 500
      CALL getarg(1,distri_file)
!     distri_file='Classical_full_model_used.txt'
!     distri_file='Resonant_full_model_used.txt'
      ext_index = index(distri_file,'.') - 1
      det_outfile = distri_file(1:ext_index)//'_RADEC.txt'
      write(*,*) 'Reading ',distri_file,' and writing output to ', 
     $ det_outfile
      open (unit=lun_h, file=det_outfile,
     $     status='REPLACE')
      open (unit=lun_t, file=distri_file, status='old')
      read (lun_t, '(a)') line
      print *, trim(line)
      write (lun_h, '(a)', advance='no') trim(line)
      write (lun_h, *) ' RA DEC Delta X Y Z'

      rierr=0
      do while ( rierr < 1 ) 
      call read_obj (distri_file, lun_t, a, e, inc, M,
     $     peri, node, h, jday, extra, co, rierr)

      obs_jday = 2459026.5
      obs_jday = 2458727.5
      mt = M
     $     + (TwoPi/(a**1.5d0*365.25d0))*(obs_jday-jday)
      mt = mt - int(mt/TwoPi)*TwoPi

      call ObsPos (code, obs_jday, obpos, tmp, ros, ierr)

      call pos_cart (a, e, inc, node, peri, M, pos(1),
     $     pos(2), pos(3))

      call RADECeclXV (pos, obpos, delta, ra, dec)

C     Only write out objects that might be in the NH Search
C     270:310,-35:-10
      if (       ( ra/drad .LT. 310)
     $     .and. ( ra/drad .GT. 270)
     $     .and. (dec/drad .LT. -10) 
     $     .and. (dec/drad .GT. -30) ) then 
      
      write (lun_h, 9000, advance='no') a, e, inc/drad, node/drad, 
     $     peri/drad, M/drad, H, jday
      do j = 1,co
         write (lun_h, 9001, advance='no') extra(j)
      end do
      write (lun_h, 9002) ra/drad, dec/drad, delta, 
     $     pos(1), pos(2), pos(3)

      end if 
      end do
      
      close (lun_h)
      close (lun_t)


 9000 format (f8.3,1x,f6.3,1x,4(f8.3,1x),f8.1,1x,f16.1,1x)
 9001 format (f8.3,1x)
 9002 format (6(f16.5,1x))
      end program xx
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

      subroutine read_obj (filen, lun_in, a, e, i, capm,
     $  om, capom, h, jday, extra, co, ierr)

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine opens and reads in the object element file.
c Angles are returned in radian.
c Potentially use a common to return the H distribution parameters.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : February 2004
c Version 2 : For L7 data release, June 2010
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     filen : object element file name
c     lun_in: File unit
c
c OUTPUT
c     a     : Semi-major axis (R8)
c     e     : Eccentricity of orbit (R8)
c     i     : Inclination (R8)
c     capm  : Mean anomaly (R8)
c     om    : Argument of pericenter (R8)
c     capom : Longitude of node (R8)
c     h     : Absolute magnitude (R8)
c     jday  : Time of elements (R8)
c     extra : Array of values (10*R8)
c     co    : number of values in 'extra'
c     ierr  : Error code
c                0 : nominal run
c               10 : unable to open filen
c               20 : error reading record
c               30 : end of file reached
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


      real*8
     $  a, e, i, capm, om, capom, h, Pi, drad, jday, extra(10), jd

      integer*4
     $  nw_max

      parameter
     $  (Pi = 3.141592653589793238d0, drad = Pi/180.0D0, nw_max = 20)

      integer
     $  lun_in, ierr, j, nw, lw(nw_max), co

      character
     $  line*200, filen*(*), word(nw_max)*80

      logical
     $  opened

      data opened /.false./

      save opened, jd

      ierr = 0
      if (.not. opened) then
         open (unit=lun_in, file=filen, status='old', err=1000)
         opened = .true.
         jd = -1.d0
      end if

 1500 continue
      do j = 1, len(line)
         line(j:j) = ' '
      end do
      read (lun_in, '(a)', err=2000, end=3000) line
      if (line(1:1) .eq. '#') then
         goto 1500
      end if
      jday = jd
      call parse (line, nw_max, nw, word, lw)
      if (nw .lt. 7) goto 2000
      read (word(1), *) a
      read (word(2), *) e
      read (word(3), *) i
      i = i*drad
      read (word(4), *) capom
      read (word(5), *) om
      read (word(6), *) capm
      read (word(7), *) h
      read (word(8), *) jday
      capom = capom*drad
      om = om*drad
      capm = capm*drad
      if (nw .gt. 19) goto 2000
      do j = 9, nw
         read(word(j), *) extra(j-9+1)
      end do
      co = nw - 9 + 1
      return

 1000 continue
      ierr = 10
      return

 2000 continue
      ierr = 20
      return

 3000 continue
      ierr = 30
      close (lun_in)
      opened = .false.
      return


      end subroutine read_obj
