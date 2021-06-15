      include 'SurveySubs.f'


      PROGRAM XX


      real*8 ai, ei, ii, noi, pei, mi
      real*8 ao, eo, io, noo, peo, mo
      real*8 deg2rad
      integer*4 ieqec, ierr, lun
      character*20 name

      deg2rad = 3.1415/180.0

      lun = 1
      write(*,*) "targetname a e incl Omega w M"
      open (unit=lun, file='orbits.txt', status='old')
 100  read(lun,*, end=999) name, ai, ei, ii, noi, pei, mi

      ii = ii*deg2rad
      noi = noi*deg2rad
      pei = pei*deg2rad
      mi = mi*deg2rad
      
      ieqec = -1


      call invar_ecl_osc(ieqec, ai, ei, ii, noi, pei, mi,
     $     ao, eo, io, noo, peo, mo, ierr)

      write(*,*) name, ao, eo, io/deg2rad, noo/deg2rad, 
     $     peo/deg2rad, mo/deg2rad

      goto 100

 999  CONTINUE
      END PROGRAM XX
