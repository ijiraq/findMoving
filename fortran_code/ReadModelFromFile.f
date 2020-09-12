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

      implicit none

      real*8
     $  a, e, i, capm, om, capom, h, Pi, drad, jday, extra(*), jd

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
         read(word(j), *) extra(j)
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

      end
