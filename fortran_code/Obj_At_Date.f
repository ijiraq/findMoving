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

      program XX

      implicit none

      integer*4 n_obj_max, screen, keybd, verbose, lun_h, lun_t
      integer*4 code, lun_trajectory, line_number
      integer*4 nw_max, nw

      parameter
     $     (n_obj_max = 10000, screen = 6, keybd = 5, verbose = 9,
     $     nw_max = 30)

      integer*4 lw(nw_max)

      real*8 Pi, TwoPi, drad

      parameter (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi,
     $  drad = Pi/180.0d0)

c color array NEEDS to be length 10 or more!
      real*8 a, e, inc, node, peri, M, h, epoch, m_int, d_ra, d_dec,
     $     r, delta, ra, dec, random, mt, extra(nw_max), gb, ph, period, 
     $     amp, lorri_mag,
     $     jday, m_rand, eff, rn_iter, eff_lim, h_rand, pos(3), 
     $     obpos(3),
     $     ros, tmp(3), obs_jday, jday_start, jday_step, jday_end, 
     $     mag, mag_best, delta_best, delta_date, llon, llat, trj_jday,
     $     ground_mag, ground_delta, alpha, eos, epos(3)

      integer*4 n_hits, n_track, ierr, seed, flag, isur, ic, n_iter,
     $  n_track_max, nchar, rierr, co, j, ext_index, start_idx, end_idx

      character distri_file*80, survey_dir*100, line*200,
     $  trk_outfile*80, det_outfile*80, buffer*400,
     $  comments*100, surna*10, trajectory_file*200,
     $     start_idx_str*80, end_idx_str*80, jday_str*80,
     $     word(nw_max)*80

      integer*4 ra_col, dec_col, earth_delta_col, earth_mag_col, 
     $     X_col, Y_col, Z_col,
     $     lorri_lon_col, lorri_lat_col, lorri_delta_col,
     $     lorri_best_mag_col, lorri_best_delta_col,
     $     lorri_best_delta_date_col, lorri_mag_col,
     $     lorri_angle_col

      logical keep_going, detectable


c 1 2 3 4     5     6 7 8       9     10 11  12          13        4 5 6 17             18               19              20        21        22
c a e i Omega omega M H epoch_M epoch RA DEC earth_delta earth_mag X Y Z lorri_mag_best lorri_delta_best best_delta_date lorri_lon lorri_lat lorri_delta

      ra_col = 10 - 8
      dec_col = 11 - 8
      earth_delta_col = 12 - 8
      earth_mag_col = 13 - 8
      X_col = 14 - 8
      Y_col = 15 - 8
      Z_col = 16 - 8
      lorri_best_mag_col = 17 - 8
      lorri_best_delta_col = 18 - 8
      lorri_best_delta_date_col = 19 - 8 
      lorri_lon_col = 20 - 8
      lorri_lat_col = 21 - 8
      lorri_delta_col = 22 - 8
      lorri_mag_col = lorri_delta_col + 1
      lorri_angle_col  = lorri_mag_col + 1

      lun_h = 10
      lun_t = 11
      lun_trajectory = 12
      keep_going = .true.
      detectable = .false.
      code = 500
      gb = 0.15

c     Date range to compute encounter condition for
      jday_start = 2459123.5 + 365
      jday_end = 2460584.5 + 365
      jday_step = 1

      CALL getarg(1, distri_file)
      distri_file = trim(distri_file)
      CALL getarg(2, trajectory_file)
      trajectory_file = trim(trajectory_file)
      CALL getarg(3, jday_str)
      read(jday_str, *) obs_jday
      ext_index = index(distri_file,'.') - 1
      det_outfile = distri_file(1:ext_index)//'_'//jday_str//'.txt'
      write(*,*) 'Reading ',distri_file,' and writing output to ', 
     $ det_outfile
      open (unit=lun_trajectory, file=trajectory_file, status='old')
      open (unit=lun_h, file=det_outfile,
     $     status='REPLACE')
      open (unit=lun_t, file=distri_file, status='old')
      read (lun_t, '(a)') line
      call parse (line, nw_max, nw, word, lw)
      do j = 1, nw
         write (lun_h, 9004, advance='no') word(j)
      end do
      write (lun_h, 9004, advance='no') 'lorri_mag'
      write (lun_h, 9004) 'lorri_angle'
      
      rierr=0

c     Get position of NH at time of TCM
      CALL FSEEK(lun_trajectory, 0, 0, ierr)
c     Skip First line (should be header line)
      read(lun_trajectory, *) buffer
      if (ierr .ne. 0) then
         return
      end if
c     First possible date for course change.
      trj_jday = obs_jday + 30*6
      call ObsPos(-lun_trajectory, trj_jday, 
     $     obpos, tmp, ros, ierr)
      call ObsPos(code, obs_jday, 
     $     epos, tmp, eos, ierr)
      do while ( rierr < 1 )
         call read_obj (distri_file, lun_t, a, e, inc, M,
     $        peri, node, h, jday, extra, co, start_idx,
     $        end_idx, rierr)
         if ( rierr .ne. 0) then
            exit
         end if
c         if ( extra(lorri_best_delta_date_col) < trj_jday) then
c            goto 8000
c         end if

c     Get the RA/DEC at survey date
         mt = M
     $        + (TwoPi/(a**1.5d0*365.25d0))*(obs_jday-jday)
         mt = mt - int(mt/TwoPi)*TwoPi
         call pos_cart(a, e, inc, node, peri, mt, pos(1),
     $        pos(2), pos(3))
         call RADECeclXV(pos, epos, ground_delta, ra, dec)
c     Compute the Ground Based magnitude
         r = sqrt(pos(1)**2 + pos(2)**2 + pos(3)**2)
         call AppMag(r, ground_delta, eos, h, gb, alpha, 
     $        ground_mag, ierr) 

c     Change the M and M_epoch to match the RA/DEC values
         M = mt
         jday = obs_jday

         extra(ra_col) = ra/drad
         extra(dec_col) = dec/drad
         extra(earth_delta_col) = ground_delta
         extra(earth_mag_col) = ground_mag
         extra(lorri_best_mag_col) = extra(lorri_best_mag_col) + 0.7
         extra(X_col) = pos(1)
         extra(Y_col) = pos(2)
         extra(Z_col) = pos(3)

c     Get the lon/lat at time of course change
         mt = M
     $        + (TwoPi/(a**1.5d0*365.25d0))*(trj_jday-jday)
         mt = mt - int(mt/TwoPi)*TwoPi
         call pos_cart(a, e, inc, node, peri, mt, pos(1),
     $        pos(2), pos(3))
         call RADECeclXV(pos, obpos, delta, llon, llat)
c     Compute the Ground Based magnitude
         r = sqrt(pos(1)**2 + pos(2)**2 + pos(3)**2)
         call AppMag(r, delta, ros, h, gb, alpha, 
     $        lorri_mag, ierr) 
         lorri_mag = lorri_mag + 0.7
         extra(lorri_lon_col) = llon/drad
         extra(lorri_lat_col) = llat/drad
         extra(lorri_delta_col) = delta
         co = co + 1
         extra(co) = lorri_mag
         co = co + 1
         extra(co) = (extra(lorri_best_delta_col)/
     $        extra(lorri_delta_col))/drad

         write (lun_h, '(8(f16.5,1x))', advance='no') a, e, inc/drad, 
     $        node/drad, peri/drad, M/drad, H, jday
         do j = 1, co
            write (lun_h, 9001, advance='no') extra(j)
         end do
         write(lun_h, *)
 8000 end do
      close (lun_h)
      close (lun_t)


 9000 format (f8.3,1x,f6.3,1x,4(f8.3,1x),f8.1,1x,f16.1,1x)
 9001 format (f16.5,1x)
 9004 format (a16,1x)
 9002 format (6(f16.5,1x))
 9003 format (13(f16.5,1x))
      end program xx
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

      subroutine read_obj (filen, lun_in, a, e, i, capm,
     $  om, capom, h, jday, extra, co, start_idx, end_idx, ierr)

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
c     start_idx : first line to read
c     end_idx : last line to read
c     ierr  : Error code
c                0 : nominal run
c               10 : unable to open filen
c               20 : error reading record
c               30 : end of file reached
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


      real*8
     $  a, e, i, capm, om, capom, h, Pi, drad, jday, extra(20), jd

      integer*4
     $  nw_max, start_idx, end_idx

      parameter
     $  (Pi = 3.141592653589793238d0, drad = Pi/180.0D0, nw_max = 29)

      integer
     $  lun_in, ierr, j, nw, lw(nw_max), co, current_line

      character
     $  line*2000, filen*(*), word(nw_max)*80

      logical
     $  opened

      data opened /.false./

      save opened, jd, current_line

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
      if (nw .gt. 29) goto 2000
      do j = 9, nw
         read(word(j), *) extra(j-8)
      end do
      co = nw - 8
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
