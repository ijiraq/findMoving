c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

      include 'SurveySubs.f'

      program XX

      implicit none

      integer*4 obs_code, ierr, niter, nobjects
      integer*4 output_lun, max_objects, max_iters
      integer*4 field_number, i, j, version

      real*8 Pi, TwoPi, drad, mu, theta, phi
      real*8 gb, alpha, mag_limit, H_limit, MA

      parameter (Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi,
     $     drad = Pi/180.0d0, mu = TwoPi**2)
      parameter (gb = 0.15, mag_limit=28, H_limit=14, version=1)
      
      real*8 obs_jday, fake_jday
      real*8 a, e, inc, node, peri, M, H, pos(6), obs_pos(3)
      real*8 epoch_M
      real*8 tmp, ros, fake_obs_pos(3), fake_ros
      real*8 ra, dec, delta, radius, mag, ddec, dra
      real*8 fake_ra, fake_dec, fake_delta, fake_mag, fake_M
      real*8 fake_ra2, fake_dec2, fake_delta2, fake_mag2, fake_M2
      real*8 fake_obs_pos2(3), fake_ros2
      real*8 obs_ra, obs_dec, obs_radius, cos_obs_dec
      real*8 obs_ra_f1, obs_ra_f2, obs_dec_f1, obs_dec_f2
      real*8 r
      character in_str*80
      character output_filename*120
      max_iters = 1E8
      obs_code = 500
C     SatEastY2 20:11:21.16 -20:22:36.0, 59029.57375
C     SatWestY2 20:05:35.75 -20:21:09.0, 59028.54766

      obs_ra_f1 = 15*(20+11/60.0+21.16/3600.0)
      obs_dec_f1 = -1*(20+22/60.0+36/3600.0)
      obs_ra_f2 = 15*(20+05/60.0+35.75/3600.0)
      obs_dec_f2 = -1*(20+21/60.0+09/3600.0)
      obs_jday = 2459030.0
      epoch_M = obs_jday
      call srand(int(epoch_M))

      CALL parse_args(field_number, fake_jday, max_objects)

      if ( field_number .eq. 1 ) then
         obs_ra = obs_ra_f1*drad
         obs_dec = obs_dec_f1*drad
      else if (field_number .eq. 2) then
         obs_ra = obs_ra_f2*drad
         obs_dec = obs_dec_f2*drad
      else 
         CALL usage("Invalid field number, must be 1 or 2", -2)
      end if

      cos_obs_dec = cos(obs_dec)
      obs_radius = 1.2*drad
      write(in_str, *) "NHF",field_number,"_",
     $     fake_jday,"_fk",version,".txt"
      j = 1
      do i = 1, len(in_str) 
         if (in_str(i:i) .ne. ' ') then
            output_filename(j:j) = in_str(i:i)
            j = j+1
         endif
      end do

      write(6,*) "Writing to "//output_filename

      output_lun = 11
      open (unit=output_lun, file=output_filename, status='new')

C     Get observatory location for primary/master location
      call ObsPos (obs_code, obs_jday, obs_pos, 
     $     tmp, ros, ierr)

C     Get the observer location at this fake date.
      call ObsPos (obs_code, fake_jday, fake_obs_pos,
     $     tmp, fake_ros, ierr)

C     Get the observer location at fake date + 1 day (for rates)
      call ObsPos (obs_code, fake_jday+1, fake_obs_pos2,
     $     tmp, fake_ros2, ierr)


C     Loop over making fake objects until we have max_objects in frame
      niter = 1
      nobjects = 1

      write(output_lun, 9002) 'obj_id', 'RA', 'DEC', 'Delta', 
     $     'mag', 'a', 
     $     'e', 'inc',
     $     'node', 'peri', 'M', 'H', 'dra', 'ddec'

      do while ((nobjects .le. max_objects) .and. 
     $     (niter .lt. max_iters))
         niter = niter + 1
         a = 30 + rand()*70
         e = rand()*0.5
         inc = Pi*rand()/2.0
         node = TwoPi*rand()
         peri = TwoPi*rand()
         M = TwoPi*rand()


         call pos_cart(a, e, inc, node, peri, M, pos(1),
     $        pos(2), pos(3))
         call RADECeclXV (pos, obs_pos, delta, ra, dec)

C     Given this random orbit determine if its on the image.
         radius = sqrt(((ra-obs_ra)*cos_obs_dec)**2 + (dec-obs_dec)**2)

         if ( radius .lt. obs_radius ) then
C     Compute the Ground Based magnitude
            r = sqrt(pos(1)**2 + pos(2)**2 + pos(3)**2)
            mag = mag_limit+1

C     Distribute H randomly between 5 and 11
C     Until we get one bright enough
            do while ( mag .gt. mag_limit )
               H = 5+rand()*9
               call AppMag(r, delta, ros, H, gb, alpha, 
     $              mag, ierr) 
            end do

C     This orbit / object is inside our sample. 
C     Compute circumstances at current observation.

            fake_M = MA(M, epoch_M, a, fake_jday)
            call pos_cart(a, e, inc, node, peri, fake_M, pos(1),
     $           pos(2), pos(3))
            call RADECeclXV (pos, fake_obs_pos, fake_delta, 
     $           fake_ra, fake_dec)
            r = sqrt(pos(1)**2 + pos(2)**2 + pos(3)**2)
            call AppMag(r, delta, fake_ros, H, gb, alpha, 
     $           mag, ierr) 

C     Compute position 1 day later to get sky motion rate
            fake_M = MA(M, epoch_M, a, fake_jday+1)

            call pos_cart(a, e, inc, node, peri, fake_M, pos(1),
     $           pos(2), pos(3))
            call RADECeclXV (pos, fake_obs_pos2, fake_delta2, 
     $           fake_ra2, fake_dec2)

            dra = (fake_ra2-fake_ra)*cos(fake_dec2)/24.0
            ddec = (fake_dec2-fake_dec)/24.0

            write (output_lun, 9001) nobjects, 
     $           fake_ra/drad, fake_dec/drad, fake_delta, mag,
     $           a, e, inc/drad, node/drad, peri/drad,
     $           fake_M/drad, H, 3600*dra/drad, 3600.0*ddec/drad

            nobjects = nobjects +1

         end if
      end do

      STOP

 9001 format (i8,13(f16.6,1x))
 9002 format (a8,13(a16,1x))

      end program xx


      real function MA(M, epoch_M, a, jday)

      real*8 M
      real*8 jday
      real*8 MA
      real*8 nu, TwoPi, epoch_M, a
      parameter (TwoPi=3.141592653589793238d0*2d0, nu=TwoPi/365.25d0)

C     Source is on the primary image and detectable by the primary.
      MA = M + (jday-epoch_M)*((nu/a**1.5))
      MA = MA - int(MA/TwoPi)*TwoPi
      return 
      end function MA
      


      subroutine usage(message, ierr)
      character*80 message
      integer*4 ierr

      write(0,*) trim(message)
      write(0,*) ""

      write(0,*) "Usage: MakeFake [1|2] JD NITER"
      write(0,*) " 1|2   --"
      write(0,*) "          1 -- Artificial source for NHF1"
      write(0,*) "          2 -- Artificial source for NHF2"
      write(0,*) " JD    -- Julian Date for image to plant into"
      write(0,*) " NITER -- Number of KBOs to generate"
      write(0,*) ""

      call exit(ierr)
      end subroutine usage


      subroutine parse_args(field_number, fake_jday, max_objects)

C     Get the JD we will generate FAKE objects positions for.
C     Usage;  MakeFake obs_ra obs_dec obs_jday fake_jday number
      character*80 in_str, message, command
      integer*4 length
      real*8 fake_jday
      integer*4 field_number
      integer*4 max_objects
      integer*4 nargs
      integer*4 argn


      nargs = COMMAND_ARGUMENT_COUNT()
      message = 'Wrong number of arguments'
      if (nargs .ne. 3) CALL usage(message, -1)

      argn = 0
      length = 80
      CALL GET_COMMAND_ARGUMENT(argn, command, length, ierr)
      if ( ierr .ne. 0 ) GOTO 1000

      argn = 1
      length = 80
      CALL GET_COMMAND_ARGUMENT(argn, in_str, length, ierr)
      if ( ierr .ne. 0 ) GOTO 1000
      read (in_str, *, err=1000) field_number
      write(6,*) field_number

      argn = 2
      CALL GET_COMMAND_ARGUMENT(argn, in_str, length, ierr)
      if (ierr .ne. 0) GOTO 1000
      read (in_str, *, err=1000) fake_jday

C     Get the number of objects to generate
      argn = 3
      CALL GET_COMMAND_ARGUMENT(argn, in_str, length, ierr)
      if (ierr .ne. 0) GOTO 1000
      read (in_str, *, err=1000) max_objects

      write(0,*) "Executing command: ",command
      write(0,*) "Field number: ", field_number
      write(0,*) "For JD:", fake_jday
      write(0,*) "Number of objects:", max_objects
      RETURN

 1000 CONTINUE

      write(message, 1001) "Error parsing argument", argn
      message = trim(message)
      CALL usage(message, -1)
 1001 FORMAT('A20,I10')
      end subroutine parse_args

