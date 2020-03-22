Program Reconstruction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Navigating through the whole database and read in the first 20 principal components of the BEL profiles (data stored in the files in the 3 folders 'WC2D', 'WC3D', 'WE3D'.
!! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!! The program loops through all physical parameters of the SBHB configuration in nested loop fashion.
!! For each SBHB configuration, the program reconstructs the BEL profile flux from the 20 PCs and stores the reconstructed profile into the vector Fnu(Nnu) where Nnu=600 is the number of frequency bins.
!! The wavelength of each frequency bin can be found from the vector lambda(1,Nnu). Note that the default emission line is H_beta with the rest wavelength 'lambdanot = 4860.09' Angstrom.
!! For other emission lines, the rest wavelength, lambdanot, needs to be updated.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer:: job
   integer, parameter:: nC = 84840, nE = 102275, njobC = 210, njobE = 240, Nnu = 600, npc = 20

   real, parameter :: pi = 3.1415926536, lambdanot = 4860.09, m = 0.93, n = 1.07, ilower = 80.0*pi/180, iupper = 100.0*pi/180
   integer:: nu, i, j
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
   real:: a, e, q, BoloF1, BoloF2, F2overF1
   real:: ophase, theta, phi, theta1, phi1, theta2, phi2

   real, dimension(Nnu) :: lambda = 0.0, Fnu, Tmeandata
   real, dimension(npc):: Tscores
   real, dimension(npc, Nnu)::Tloads

   integer:: ia, ie, iRin1, iRin2, itheta1, itheta2, iphi1, iphi2, iq, itheta, iphi, iophase, taujob

   character(*), parameter:: dfolder1 = 'WC2D/', dfolder3 = 'WC3D/', dfolder5 = 'WE3D/'
   character(12):: filename

   !Wavelength values of the frequency bins are stored in the array lambda(1,Nnu).
   !The plot of the BEL profile flux in term of wavelength can be done by pairing the flux array Fnu(1,Nnu) with the lambda(1,Nnu) array.
   do nu = 1, Nnu
      lambda(nu) = m + (nu - 1)*dX1
      lambda(nu) = lambdanot/lambda(nu)
   enddo

   ! Loading the first 20 principal eigenprofiles from the text file Tloads20.txt
   open (unit=9, file='Tloads20.txt')
   read (9, *) Tloads
   close (9)

   ! Loading the average profile of the whole database from the text file Tmeandata.txt
   open (unit=10, file='Tmeandata.txt')
   read (10, *) Tmeandata
   close (10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Working on the WC2D folder with a total of 120 jobs, 4 tau values 0.0001, 0.1, 1.0, 100.0 are stored in 4 binary files 'WC2D20pc1.bi' to 'WC2D20pc4.bi' !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do taujob = 1, 4
      write (filename, 300) 'WC2D20pc', taujob, '.bi'
300   format(A, I1, A)
      open (unit=1, file=dfolder1//filename, form='unformatted', access='stream', status='unknown')
      do job = 100, 219 !100, 219
         if (1 + mod(int((job - 100)/(size(Aq)*size(Aa))), size(AOptDepthC)) == taujob) then
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
                                                   !Reading from the database the first 20 PC of the profile flux.
                                                   read (1) Tscores
                                                   !Reconstructed BEL profile is stored in the array Fnu(1,Nnu)
                                                   Fnu = 0.0
                                                   do i = 1, Nnu
                                                      do j = 1, npc
                                                         Fnu(i) = Fnu(i) + Tscores(j)*Tloads(j, i)
                                                      enddo !j
                                                      Fnu(i) = Fnu(i) + Tmeandata(i)
                                                   enddo !i
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
         endif
      enddo !job
      close (1)
   enddo !taujob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Working on the WC3D folder with a total of 90 jobs, 3 tau values 0.1,1.0, 100.0 are stored in 3 binary files 'WC3D20pc2.bi' to 'WC3D20pc4.bi' !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do taujob = 2, 4
      write (filename, 300) 'WC3D20pc', taujob, '.bi'
      open (unit=1, file=dfolder3//filename, form='unformatted', access='stream', status='unknown')
      do job = 100, 219
         if (1 + mod(int((job - 100)/(size(Aq)*size(Aa))), size(AOptDepthC)) == taujob) then
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
                                                   !Reading from the database the first 20 PC of the profile flux.
                                                   read (1) Tscores
                                                   !Reconstructed BEL profile is stored in the array Fnu(1,Nnu)
                                                   Fnu = 0.0
                                                   do i = 1, Nnu
                                                      do j = 1, npc
                                                         Fnu(i) = Fnu(i) + Tscores(j)*Tloads(j, i)
                                                      enddo !j
                                                      Fnu(i) = Fnu(i) + Tmeandata(i)
                                                   enddo !i
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
         endif
      enddo !job
      close (1)
   enddo !taujob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Working on the WE3D folder with a total of 240 jobs, 2 tau values 0.0001, 1.0 are stored in 2 binary files 'WE3D20pc1.bi' and 'WE3D20pc2.bi' !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do taujob = 1, 2
      write (filename, 300) 'WE3D20pc', taujob, '.bi'
      open (unit=1, file=dfolder5//filename, form='unformatted', access='stream', status='unknown')
      do job = 100, 459
         if (1 + mod(int((job - 100)/(size(Aq)*size(Aa))), size(AOptDepthE)) == taujob) then
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
                                           Rin1 = ARin1(1 + mod(int((job - 100)/(size(Aq)*size(Aa)*size(AOptDepthE))), size(ARin1)))
                               Rin2 = ARin2(1 + mod(int((job - 100)/(size(Aq)*size(Aa)*size(AOptDepthE)*size(ARin1))), size(ARin2)))
                                             !Reading from the database the BEL profile flux in terms of frequency of the binary system.
                                             !Reading from the database the first 20 PC of the profile flux.
                                             read (1) Tscores
                                             !Reconstructed BEL profile is stored in the array Fnu(1,Nnu)
                                             Fnu = 0.0
                                             do i = 1, Nnu
                                                do j = 1, npc
                                                   Fnu(i) = Fnu(i) + Tscores(j)*Tloads(j, i)
                                                enddo !j
                                                Fnu(i) = Fnu(i) + Tmeandata(i)
                                             enddo !i
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
         endif
      enddo !job
      close (1)

   enddo !taujob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
End Program Reconstruction
