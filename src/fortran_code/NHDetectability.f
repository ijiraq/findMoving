c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      
      include 'SurveySubs.f'

      program XX

      implicit none

      integer*4 n_obj_max, screen, keybd, verbose, lun_h, lun_t
      integer*4 code, lun_trajectory, line_number

      parameter
     $  (n_obj_max = 10000, screen = 6, keybd = 5, verbose = 9)

      real*8 Pi, TwoPi, drad

      parameter (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi,
     $  drad = Pi/180.0d0)

c color array NEEDS to be length 10 or more!
      real*8 a, e, inc, node, peri, M, h, epoch, m_int, d_ra, d_dec,
     $     r, delta, ra, dec, random, mt, extra(10), gb, ph, period, 
     $     amp,
     $     jday, m_rand, eff, rn_iter, eff_lim, h_rand, pos(3), 
     $     obpos(3),
     $     ros, tmp(3), obs_jday, jday_start, jday_step, jday_end, 
     $     mag, mag_best, delta_best, delta_date, llon, llat, trj_jday,
     $     ground_mag, ground_delta, alpha

      integer*4 n_hits, n_track, ierr, seed, flag, isur, ic, n_iter,
     $  n_track_max, nchar, rierr, co, j, ext_index, start_idx

      character distri_file*80, survey_dir*100, line*200,
     $  trk_outfile*80, det_outfile*80, buffer*400,
     $  comments*100, surna*10, trajectory_file*200,
     $     start_idx_str*80

      logical keep_going, detectable

      lun_h = 10
      lun_t = 11
      lun_trajectory = 12
      keep_going = .true.
      detectable = .false.
      code = 500
      gb = 0.15

c     Date range to compute encounter condition for
      jday_start = 2459123.5
      jday_end = 2460584.5
      jday_step = 1

      distri_file = 'model.txt' 
c      CALL getarg(1, distri_file)
c      distri_file = trim(distri_file)
      trajectory_file = 'nhtraj_2019_2024.txt'
c     CALL getarg(2, trajectory_file)
c      trajectory_file = trim(trajectory_file)
      CALL getarg(1, start_idx_str)
      start_idx_str = trim(start_idx_str)
      READ(start_idx_str, *) start_idx
      

      ext_index = index(distri_file,'.') - 1
      write(det_outfile,'(a,a,i0.8,a)') 'nh_observable',
     $ '_', start_idx, '.txt'
      write(*,*) 'Reading ',distri_file,' and writing output to ', 
     $ det_outfile
      open (unit=lun_trajectory, file=trajectory_file, status='old')
      open (unit=lun_h, file=det_outfile,
     $     status='REPLACE')
      open (unit=lun_t, file=distri_file, status='old')
      read (lun_t, '(a)') line
      write (lun_h, '(a)', advance='no') trim(line)
      write (lun_h, '(a)', advance='no') ' RA DEC earth_delta earth_mag'
      write (lun_h, '(a)', advance='no') ' X Y Z '
      write (lun_h, '(a)', advance='no') 'lorri_mag_best '
      write (lun_h, '(a)', advance='no') 'lorri_delta_best '
      write (lun_h, '(a)', advance='no') 'best_delta_date lorri_lon '
      write (lun_h, '(a)') 'lorri_lat lorri_delta'
      
      rierr=0
      do while ( rierr < 1 )
c     2020-06-20T00:00:00  - date the survey will be done.
      obs_jday = 2459020.5
c     time.Time("2036-10-01 00:00:00").jd
c      obs_jday = 2464967.5
c     First possible date for course change.
         trj_jday = obs_jday + 30*4
         call read_obj (distri_file, lun_t, a, e, inc, M,
     $        peri, node, h, jday, extra, co,
     $        rierr)
         
         mt = M
     $        + (TwoPi/(a**1.5d0*365.25d0))*(obs_jday-jday)
         mt = mt - int(mt/TwoPi)*TwoPi         
         call ObsPos (code, obs_jday, obpos, tmp, ros, ierr)
c     Compute position at mt (Mean Anomally at Time)
         call pos_cart (a, e, inc, node, peri, mt, pos(1),
     $        pos(2), pos(3))
         r = sqrt(pos(1)**2 + pos(2)**2 + pos(3)**2)
c     Compute RA/DEC at observation, viewed from Earth
         call RADECeclXV (pos, obpos, ground_delta, ra, dec)

C
C     Only check  objects that might be in the NH Search
C     270:310,-35:-10
C    
        if (( ra/drad .LT. 310)
     $        .and. ( ra/drad .GT. 270)
     $        .and. (dec/drad .LT. -10) 
     $        .and. (dec/drad .GT. -30) ) then 
           
c     Compute the Ground Based magnitude
           call AppMag(r, ground_delta, ros, h, gb, alpha, 
     $          ground_mag, ierr) 
           if ( ground_mag .lt. 28 ) then 
C     Determine best time to observe this source
              obs_jday = jday_start
              detectable = .False.
              mag_best = 1E9
              delta_best = 1E9
              delta_date = 0
              CALL FSEEK(lun_trajectory, 0, 0, ierr)
c     Skip First line (should be header line)
              read(lun_trajectory, *) buffer
              if (ierr .ne. 0) then
                 return
              end if
              
              do while (obs_jday < jday_end)
                 mt = M
     $                + (TwoPi/(a**1.5d0*365.25d0))*(obs_jday-jday)
                 mt = mt - int(mt/TwoPi)*TwoPi
                 call ObsPos(-lun_trajectory, obs_jday, 
     $                obpos, tmp, ros, ierr)
                 if ( ierr .ne. 0 ) then
                    write(screen, *) 'obs pos failed'
                    exit
                 end if
                 call pos_cart(a, e, inc, node, peri, mt, pos(1),
     $                pos(2), pos(3))
                 r = sqrt(pos(1)**2 + pos(2)**2 + pos(3)**2)
                 call RADECeclXV(pos, obpos, delta, llon, llat)
                 call AppMag(r, delta, ros, h, gb, alpha, mag, ierr) 
                 if ( ierr .ne. 0 ) then
                    write(screen, *) 'AppMag Failed'
                    exit
                 end if
                 if ( mag < 23 ) then 
                    detectable = .True.
                    if (mag < mag_best) then
                       mag_best = mag
                    end if
                    if (delta < delta_best) then
                       delta_best = delta
                       delta_date = obs_jday
                    end if
                 end if
                 obs_jday = obs_jday + jday_step
              end do
              
              if ( detectable ) then 
                 
c     Get the lon/lat at time of course change
                 mt = M
     $                + (TwoPi/(a**1.5d0*365.25d0))*(trj_jday-jday)
                 mt = mt - int(mt/TwoPi)*TwoPi
                 call ObsPos(-lun_trajectory, trj_jday, 
     $                obpos, tmp, ros, ierr)
                 call pos_cart(a, e, inc, node, peri, mt, pos(1),
     $                pos(2), pos(3))
                 call RADECeclXV(pos, obpos, delta, llon, llat)
                 
                 write (lun_h, 9000, advance='no') a, e, inc/drad, 
     $                node/drad, peri/drad, M/drad, H, jday
                 do j = 1,co
                    write (lun_h, 9001, advance='no') extra(j)
                 end do
                 write (lun_h, 9003) ra/drad, dec/drad, 
     $                ground_delta, ground_mag,
     $                pos(1), pos(2), pos(3),
     $                mag_best, delta_best, delta_date, 
     $                llon/drad, llat/drad, delta               
              end if
           end if
        end if
      end do
      close (lun_h)
      close (lun_t)
      

 9000 format (f8.3,1x,f6.3,1x,4(f8.3,1x),f8.1,1x,f16.1,1x)
 9001 format (f8.3,1x)
 9002 format (6(f16.5,1x))
 9003 format (13(f16.5,1x))
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
c     drop out if we've read the desired number of lines
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
