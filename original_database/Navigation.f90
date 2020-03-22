Program Navigation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Navigating through the whole database starting from the first file 'file100.bi' in folder 'WC2D', through folder 'WC3D', and ending at the last file 'file429.bi' in the folder 'WE3D'.
!! WC2D folder contains 120 binary '.bi' files, each file contains nC =  84,840 SBHB configurations with circular orbits where accretion disk wind occurs only on the two mini-disks.
!! WC3D folder contains  90 binary '.bi' files, each file contains nC =  84,840 SBHB configurations with circular orbits where accretion disk wind occurs on all three disks.
!! WE3D folder contains 240 binary '.bi' files, each file contains nE = 102,275 SBHB configurations with eccentric orbits where accretion disk wind occurs on all three disks.
!! Each SBHB configuration contains the contribution to the total BEL profile from 3 sources: the BLR around the 1) primary mini-disk, 2) secondary mini-disk, 3) circumbinary disk.
!! The total BEL profile of the whole SBHB system can be found by simply adding up the 3 contributions listed above.
!! Each profile contains Nnu = 600 frequency bins. The profile flux values are stored as single precision float numbers.
!! Hence, each '.bi' file in the WC2D, WC3D folder contains 'nC * Nnu * 3' single precision float values.
!! And,   each '.bi' file in the WE3D folder       contains 'nE * Nnu * 3' single precision float values.
!! Since each single precision float value uses a memory of exactly 4 bytes, the total size of the whole database can be calculated as:
!! 120 * 84,840 * 600 * 3 * 4 (bytes) + 90 * 84,840 * 600 * 3 * 4 (bytes) + 2400 * 102,275 * 600 * 3 * 4 (bytes) =  305,009,280,000 (bytes)

!! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!! The script loops through all physical parameters of the SBHB configuration in nested loop fashion.
!! For each SBHB configuration, the program reads and adds up the 3 contributions from the 3 disks and then stores the BEL profile flux into the array Fnu(Nnu) where Nnu=600 is the number of frequency bins.
!! The wavelength of each frequency bin can be found from the array lambda(Nnu). Note that the default emission line is H_beta with the rest wavelength 'lambdanot = 4860.09' Angstrom.
!! For other emission lines, the rest wavelength, lambdanot, needs to be updated.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer:: job
   integer, parameter:: nC = 84840, nE = 102275, njobC = 210, njobE = 240, Nnu = 600

   real, parameter :: pi = 3.1415926536, lambdanot = 4860.09, m = 0.93, n = 1.07, ilower = 80.0*pi/180, iupper = 100.0*pi/180
   integer:: nu
   real, parameter:: dX1 = (n - m)/(Nnu - 1)

   real, dimension(4):: Atheta = (/5.0, 155.0, 55.0, 105.0/)*pi/180

   real, dimension(1):: Aphi = (/0.0/)*pi/180
   real, dimension(5):: AphiE = (/36.0, 108.0, 180.0, 242.0, 324.0/)*pi/180

   real, dimension(6):: Atheta1 = (/0.0, 165.0, 30.0, 135.0, 60.0, 105.0/)*pi/180
   real, dimension(6):: Atheta2 = (/0.0, 165.0, 30.0, 135.0, 60.0, 105.0/)*pi/180

   real, dimension(6):: Aphi1 = (/0.0, 25.0, 60.0, 185.0, 210.0, 235.0/)*pi/180
   real, dimension(6):: Aphi2 = (/0.0, 25.0, 60.0, 185.0, 210.0, 235.0/)*pi/180

   real, dimension(5):: Aa = (/5e3, 1e4, 5e4, 1e5, 1e6/)
   real, dimension(6):: Aq = (/1.0, 9.0/11, 2.0/3, 3.0/7, 1.0/3, 1.0/10/)

   real, dimension(1):: Ae = (/0.0/)
   real, dimension(1):: AeE = (/0.5/)

   real, dimension(2):: ARin1 = (/500.0, 1000.0/)
   real, dimension(2):: ARin2 = (/500.0, 1000.0/)

   real, dimension(4):: AOptDepthC = (/0.0001, 0.1, 1.0, 100.0/)
   real, dimension(3):: AOptDepthE = (/0.0001, 1.0, 1.0/)

   real, dimension(1):: AEta = (/1.0/)
   real, dimension(1):: AWindAngle = (/10.0/)*pi/180

   real, dimension(5):: Aophase = (/0.0, 72.0, 144.0, 216.0, 288.0/)*pi/180

   real:: Rin1, Rin2, Rin3, Rout1, Rout2, Rout3, inc1, inc2
   real:: a, e, q
   real:: ophase, theta, phi, theta1, phi1, theta2, phi2

   real, dimension(Nnu) :: lambda = 0.0, Fnu = 0.0

   real :: A1 = 0.0, A2 = 0.0, A3 = 0.0

   integer:: ia, ie, iRin1, iRin2, itheta1, itheta2, iphi1, iphi2, iq, itheta, iphi, iophase, taujob

   character(10):: filename

   !Wavelength values of the frequency bins are stored in the array lambda(1,Nnu).
   !The plot of the BEL profile flux in term of wavelength can be done by pairing the flux array Fnu(1,Nnu) with the lambda(1,Nnu) array.
   do nu = 1, Nnu
      lambda(nu) = m + (nu - 1)*dX1
      lambda(nu) = lambdanot/lambda(nu)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Working on the WC2D folder with a total of 120 files (jobs are selected continuously from job=100,219), 4 tau values 0.0001, 0.1, 1.0, 100.0 !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do taujob = 1, 4
      do job = 100, 219
         if (1 + mod(int((job - 100)/(size(Aq)*size(Aa))), size(AOptDepthC)) == taujob) then
            write (filename, 300) 'file', job, '.bi'
300         format(A, I3, A)
            open (unit=job*1, file='WC2D'//filename, form='unformatted', access='stream', status='unknown')
! Orientation variables: theta, phi, theta1, phi1, theta2, phi2 - affect inclination angles inc1, inc2.
            do itheta = 1, size(Atheta)
               theta = Atheta(itheta)
               do iphi = 1, size(Aphi)
                  phi = Aphi(iphi)
                  do itheta1 = 1, size(Atheta1)
                     theta1 = Atheta1(itheta1)
                     do iphi1 = 1, size(Aphi1)
                        phi1 = Aphi1(iphi1)
                        do itheta2 = 1, size(Atheta2)
                           theta2 = Atheta2(itheta2)
                           do iphi2 = 1, size(Aphi2)
                              phi2 = Aphi2(iphi2)
                              inc1 = acos(sin(theta1)*cos(phi1)*sin(theta)*cos(phi) + &
                                          sin(theta1)*sin(phi1)*sin(theta)*sin(phi) + cos(theta1)*cos(theta))
                              inc2 = acos(sin(theta2)*cos(phi2)*sin(theta)*cos(phi) + &
                                          sin(theta2)*sin(phi2)*sin(theta)*sin(phi) + cos(theta2)*cos(theta))
                              if ((ilower > inc1 .or. inc1 > iupper) .and. (ilower > inc2 .or. inc2 > iupper)) then
                                 ! Orbital phase or true anomaly: omega - affect the change of frame.
                                 do iophase = 1, size(Aophase)
                                    ophase = Aophase(iophase)
! Circular or Elliptical orbits: e
                                    do ie = 1, size(Ae)
                                       e = Ae(ie)
! Semi-major axis a, and mass ratio q, are determined by the job number. The Semi-major axis a is in terms of total mass: a = a(M).
                                       q = Aq(1 + mod((job - 100), Size(Aq)))
                                       a = Aa(1 + mod(int((job - 100)/size(Aq)), size(Aa)))
                                       !The circumbinary disk size depends on the semi-major axis a.
                                       Rin3 = 2.0*a
                                       Rout3 = 3.0*a
                                       Rout2 = Paczynski(a, e, q)*(1.0 + q)/q
                                       Rout1 = Paczynski(a, e, 1.0/q)*(1.0 + q)
                                       OptDepth = AOptDepthC(1 + mod(int((job - 100)/(size(Aq)*size(Aa))), size(AOptDepthC)))
                                       do iEta = 1, size(AEta)
                                          Eta = AEta(iEta)
                                          do iWindAngle = 1, size(AWindAngle)
                                             WindAngle = AWindAngle(iWindAngle)
                                             do iRin1 = 1, size(ARin1)
                                                Rin1 = ARin1(iRin1)
                                                do iRin2 = 1, size(ARin2)
                                                   Rin2 = ARin2(iRin2)
                                                   !Reading from the database the BEL profile flux in terms of frequency of the binary system.
                                                   do nu = 1, Nnu
                                                      read (job*1) A1, A2, A3
                                                      Fnu(nu) = A1 + A2 + A3
                                                   enddo
                                                enddo ! Rin2
                                             enddo ! Rin1
                                          enddo !WindAngle
                                       enddo !Eta
                                    enddo ! e
                                 enddo ! ophase
                              endif
                           enddo ! phi2
                        enddo ! theta2
                     enddo ! phi1
                  enddo ! theta1
               enddo ! phi
            enddo !theta
            close (job*1)
         endif
      enddo !job
   enddo !taujob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Working on the WC3D folder with a total of 90 files (90 jobs are selected NOT continuously from job=100,219), 3 tau values 0.1,1.0, 100.0 !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do taujob = 2, 4
      do job = 100, 219
         if (1 + mod(int((job - 100)/(size(Aq)*size(Aa))), size(AOptDepthC)) == taujob) then
            write (filename, 300) 'file', job, '.bi'
            open (unit=job*3, file='WC3D'//filename, form='unformatted', access='stream', status='unknown')
! Orientation variables: theta, phi, theta1, phi1, theta2, phi2 - affect inclination angles inc1, inc2.
            do itheta = 1, size(Atheta)
               theta = Atheta(itheta)
               do iphi = 1, size(Aphi)
                  phi = Aphi(iphi)
                  do itheta1 = 1, size(Atheta1)
                     theta1 = Atheta1(itheta1)
                     do iphi1 = 1, size(Aphi1)
                        phi1 = Aphi1(iphi1)
                        do itheta2 = 1, size(Atheta2)
                           theta2 = Atheta2(itheta2)
                           do iphi2 = 1, size(Aphi2)
                              phi2 = Aphi2(iphi2)
                              inc1 = acos(sin(theta1)*cos(phi1)*sin(theta)*cos(phi) + &
                                          sin(theta1)*sin(phi1)*sin(theta)*sin(phi) + cos(theta1)*cos(theta))
                              inc2 = acos(sin(theta2)*cos(phi2)*sin(theta)*cos(phi) + &
                                          sin(theta2)*sin(phi2)*sin(theta)*sin(phi) + cos(theta2)*cos(theta))
                              if ((ilower > inc1 .or. inc1 > iupper) .and. (ilower > inc2 .or. inc2 > iupper)) then
                                 ! Orbital phase or true anomaly: ophase - affect the change of frame.
                                 do iophase = 1, size(Aophase)
                                    ophase = Aophase(iophase)
! Circular or Elliptical orbits: e
                                    do ie = 1, size(Ae)
                                       e = Ae(ie)
                                       q = Aq(1 + mod((job - 100), Size(Aq)))
                                       a = Aa(1 + mod(int((job - 100)/size(Aq)), size(Aa)))
                                       Rin3 = 2.0*a
                                       Rout3 = 3.0*a
                                       Rout2 = Paczynski(a, e, q)*(1.0 + q)/q
                                       Rout1 = Paczynski(a, e, 1.0/q)*(1.0 + q)
                                       OptDepth = AOptDepthC(1 + mod(int((job - 100)/(size(Aq)*size(Aa))), size(AOptDepthC)) + 1) ! offset by 1 since no DW case is not included in the 3DW model.
                                       do iEta = 1, size(AEta)
                                          Eta = AEta(iEta)
                                          do iWindAngle = 1, size(AWindAngle)
                                             WindAngle = AWindAngle(iWindAngle)
                                             do iRin1 = 1, size(ARin1)
                                                Rin1 = ARin1(iRin1)
                                                do iRin2 = 1, size(ARin2)
                                                   Rin2 = ARin2(iRin2)
                                                   !Reading from the database the BEL profile flux in terms of frequency of the binary system.
                                                   do nu = 1, Nnu
                                                      read (job*3) A1, A2, A3
                                                      Fnu(nu) = A1 + A2 + A3
                                                   enddo
                                                enddo ! Rin2
                                             enddo ! Rin1
                                          enddo !WindAngle
                                       enddo !Eta
                                    enddo ! e
                                 enddo ! ophase
                              endif
                           enddo ! phi2
                        enddo ! theta2
                     enddo ! phi1
                  enddo ! theta1
               enddo ! phi
            enddo !theta
            close (job*3)
         endif
      enddo !job
   enddo !taujob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Working on the WE3D folder with a total of 240 files (240 jobs are selected NOT continuously from job=100,459), 2 tau values 0.0001, 1.0 !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do taujob = 1, 2
      do job = 100, 459
         if (1 + mod(int((job - 100)/(size(Aq)*size(Aa))), size(AOptDepthE)) == taujob) then
            write (filename, 300) 'file', job, '.bi'
            open (unit=job*5, file='WE3D'//filename, form='unformatted', access='stream', status='unknown')
! Orientation variables: theta, phi, theta1, phi1, theta2, phi2 - affect inclination angles inc1, inc2.
            do itheta = 1, size(Atheta)
               theta = Atheta(itheta)
               do iphi = 1, size(AphiE)
                  phi = AphiE(iphi)
                  do itheta1 = 1, size(Atheta1)
                     theta1 = Atheta1(itheta1)
                     do iphi1 = 1, size(Aphi1)
                        phi1 = Aphi1(iphi1)
                        do itheta2 = 1, size(Atheta2)
                           theta2 = Atheta2(itheta2)
                           do iphi2 = 1, size(Aphi2)
                              phi2 = Aphi2(iphi2)
                              inc1 = acos(sin(theta1)*cos(phi1)*sin(theta)*cos(phi) + &
                                          sin(theta1)*sin(phi1)*sin(theta)*sin(phi) + cos(theta1)*cos(theta))
                              inc2 = acos(sin(theta2)*cos(phi2)*sin(theta)*cos(phi) + &
                                          sin(theta2)*sin(phi2)*sin(theta)*sin(phi) + cos(theta2)*cos(theta))
                              if ((ilower > inc1 .or. inc1 > iupper) .and. (ilower > inc2 .or. inc2 > iupper)) then
                                 ! Orbital phase or true anomaly: ophase - affect the change of frame.
                                 do iophase = 1, size(Aophase)
                                    ophase = Aophase(iophase)
! Circular or Elliptical orbits: e
                                    do ie = 1, size(AeE)
                                       e = AeE(ie)
                                       q = Aq(1 + mod((job - 100), Size(Aq)))
                                       a = Aa(1 + mod(int((job - 100)/size(Aq)), size(Aa)))
                                       Rin3 = 2.0*a
                                       Rout3 = 3.0*a
                                       Rout2 = Paczynski(a, e, q)*(1.0 + q)/q
                                       Rout1 = Paczynski(a, e, 1.0/q)*(1.0 + q)
                                       if (taujob == 1) OptDepth = 0.0001
                                       if (taujob == 2) OptDepth = 1.0
                                       do iEta = 1, size(AEta)
                                          Eta = AEta(iEta)
                                          do iWindAngle = 1, size(AWindAngle)
                                             WindAngle = AWindAngle(iWindAngle)
                                             Rin1 = ARin1(1 + mod(int((job - 100)/(size(Aq)*size(Aa) &
                                                                                   *size(AOptDepthE))), size(ARin1)))
                                             Rin2 = ARin2(1 + mod(int((job - 100)/(size(Aq)*size(Aa) &
                                                                                   *size(AOptDepthE)*size(ARin1))), size(ARin2)))
                                             !Reading from the database the BEL profile flux in terms of frequency of the binary system.
                                             do nu = 1, Nnu
                                                read (job*5) A1, A2, A3
                                                Fnu(nu) = A1 + A2 + A3
                                             enddo
                                          enddo !WindAngle
                                       enddo !Eta
                                    enddo ! e
                                 enddo ! ophase
                              endif
                           enddo ! phi2
                        enddo ! theta2
                     enddo ! phi1
                  enddo ! theta1
               enddo ! phi
            enddo !theta
            close (job*5)
         endif
      enddo !job
   enddo !taujob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS
!-----------------------------------------------------------------------
   Function Paczynski(a, e, q)
      real:: a, e, q
      If (q == 10.0) then
         Paczynski = 0.4897977*a*(1.0 - e)
      else
         if (q == 9.0) then
            Paczynski = 0.482667*a*(1.0 - e)
         else
            if (q == 3.0) then
               Paczynski = 0.398667*a*(1.0 - e)
            else
               if (q == 7.0/3) then
                  Paczynski = 0.377000*a*(1.0 - e)
               else
                  if (q == 3.0/2) then
                     Paczynski = 0.339000*a*(1.0 - e)
                  else
                     if (q == 11.0/9) then
                        Paczynski = 0.300667*a*(1.0 - e)
                     else
                        if (q == 1.0) then
                           Paczynski = 0.276667*a*(1.0 - e)
                        else
                           if (q == 9.0/11) then
                              Paczynski = 0.255333*a*(1.0 - e)
                           else
                              if (q == 2.0/3) then
                                 Paczynski = 0.237333*a*(1.0 - e)
                              else
                                 if (q == 3.0/7) then
                                    Paczynski = 0.202333*a*(1.0 - e)
                                 else
                                    if (q == 1.0/3) then
                                       Paczynski = 0.185667*a*(1.0 - e)
                                    else
                                       if (q == 1.0/9) then
                                          Paczynski = 0.129000*a*(1.0 - e)
                                       else
                                          if (q == 1.0/10) then
                                             Paczynski = 0.125422*a*(1.0 - e)
                                          else
                                             Print *, "Error mass ratio."
                                          endif
                                       endif
                                    endif
                                 endif
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      endif
   End Function Paczynski
!----------------------------------------------------------------------
End Program Navigation
