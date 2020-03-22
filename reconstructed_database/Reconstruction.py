import numpy as np
import pandas as pd
import struct


# Constants:
nC = 84840
nE = 102275
njobC = 210
njobE = 240
Nnu = 600
npc = 20
lambdanot = 4860.09
m = 0.93
n = 1.07
#ilower = np.around(80.0*np.pi/180,decimals=9)
#iupper = np.around(100.0*np.pi/180,decimals=9)
ilower = 1.3962634 # ~ 80 degree
iupper = 1.7453293 # ~100 degree
dX1 = (n - m)/(Nnu - 1)

# Arrays of Parameter Choices:
Atheta = np.array([5.0, 155.0, 55.0, 105.0], dtype='f')*np.pi/180.0
AphiC = np.array([0.0], dtype='f')*np.pi/180
AphiE = np.array([36.0, 108.0, 180.0, 242.0, 324.0], dtype='f')*np.pi/180.0
Atheta1 = np.array([0.0, 165.0, 30.0, 135.0, 60.0, 105.0], dtype='f')*np.pi/180.0
Atheta2 = np.array([0.0, 165.0, 30.0, 135.0, 60.0, 105.0], dtype='f')*np.pi/180.0
Aphi1 = np.array([0.0, 25.0, 60.0, 185.0, 210.0, 235.0], dtype='f')*np.pi/180.0
Aphi2 = np.array([0.0, 25.0, 60.0, 185.0, 210.0, 235.0], dtype='f')*np.pi/180.0
Aa = np.array([5e3, 1e4, 5e4, 1e5, 1e6])
Aq = np.array([1.0, 9.0/11, 2.0/3, 3.0/7, 1.0/3, 1.0/10])
AeC = np.array([0.0])
AeE = np.array([0.5])
ARin1 = np.array([500.0, 1000.0])
ARin2 = np.array([500.0, 1000.0])
AOptDepthC = np.array([0.0001, 0.1, 1.0, 100.0])
AOptDepthE = np.array([0.0001, 1.0, 1.0])
AEta = np.array([1.0])
AWindAngle = np.array([10.0])*np.pi/180
Aophase = np.array([0.0, 72.0, 144.0, 216.0, 288.0], dtype='f')*np.pi/180.0

# A few array
wl = np.empty(Nnu)
Fnu = np.empty(Nnu)
Tscores=np.empty(npc)
Tloads=np.empty((Nnu,npc))
Tmeandata=np.empty(Nnu)

# Convert frequency bins (nu) to wavelength (wl).

for nu in range(0, Nnu):
 wl[nu] = m + (nu - 1)*dX1
 wl[nu] = lambdanot/wl[nu]

# Loading the first 20 principal eigenprofiles from the text file Tloads20.txt
Tloads=pd.read_csv('Tloads20.txt', delim_whitespace=True, header=None).astype(float).values

#Loading the average profile of the whole database from the text file Tmeandata.txt
Tmeandata=pd.read_csv('Tmeandata.txt', delim_whitespace=True, header=None).astype(float).values

# Some functions:
def Paczynski(a, e, q):
    if np.isclose(q, 10.0):
        return 0.4897977*a*(1.0 - e)
    elif np.isclose(q, 9.0):
        return 0.482667*a*(1.0 - e)
    elif np.isclose(q, 3.0):
        return 0.398667*a*(1.0 - e)
    elif np.isclose(q, 7.0/3):
        return 0.377000*a*(1.0 - e)
    elif np.isclose(q, 3.0/2):
        return 0.339000*a*(1.0 - e)
    elif np.isclose(q, 11.0/9):
        return 0.300667*a*(1.0 - e)
    elif np.isclose(q, 1.0):
        return 0.276667*a*(1.0 - e)
    elif np.isclose(q, 9.0/11):
        return 0.255333*a*(1.0 - e)
    elif np.isclose(q, 2.0/3):
        return 0.237333*a*(1.0 - e)
    elif np.isclose(q, 3.0/7):
        return 0.202333*a*(1.0 - e)
    elif np.isclose(q, 1.0/3):
        return 0.185667*a*(1.0 - e)
    elif np.isclose(q, 1.0/9):
        return 0.129000*a*(1.0 - e)
    elif np.isclose(q, 1.0/10):
        return 0.125422*a*(1.0 - e)
    else:
        print('Error mass ratio.')





###### Working on the WC2D folder with a total of 120 jobs, 4 tau values 0.0001, 0.1, 1.0, 100.0 are stored in 4 binary files 'WC2D20pc1.bi' to 'WC2D20pc4.bi' ########

counting=0

for taujob in range(0, 4): #for taujob in range(0, 4):
 with open('WC2D/WC2D20pc'+str(taujob+1)+'.bi','rb') as f:		
  for job in range(100, 220): #for job in range(100, 220):
   if int((job - 100)/(len(Aq)*len(Aa))) % len(AOptDepthC) == taujob:           
    # Orientation variables: theta, phi, theta1, phi1, theta2, phi2 - affect inclination angles inc1, inc2.
    for theta in Atheta:
     for phi in AphiC:
      for theta1 in Atheta1:
       for phi1 in Aphi1:
        for theta2 in Atheta2:
         for phi2 in Aphi2:
          inc1 = np.arccos(np.sin(theta1)*np.cos(phi1)*np.sin(theta)*np.cos(phi) + np.sin(
                 theta1)*np.sin(phi1)*np.sin(theta)*np.sin(phi) + np.cos(theta1)*np.cos(theta))
          inc2 = np.arccos(np.sin(theta2)*np.cos(phi2)*np.sin(theta)*np.cos(phi) + np.sin(
                 theta2)*np.sin(phi2)*np.sin(theta)*np.sin(phi) + np.cos(theta2)*np.cos(theta))
          if (ilower > inc1 or inc1 > iupper) and (ilower > inc2 or inc2 > iupper):
           # Orbital phase or true anomaly: omega - affect the change of frame.
           for ophase in Aophase:
            # Circular or Elliptical orbits: e
            for e in AeC:
             # Semi-major axis a, and mass ratio q, are determined by the job number. The Semi-major axis a is in terms of total mass: a = a(M).
             q = Aq[(job-100) % len(Aq)]
             a = Aa[int((job-100)/len(Aq)) % len(Aa)]
             # The circumbinary disk size depends on the semi-major axis a.
             Rin3 = 2.0*a
             Rout3 = 3.0*a
             Rout2 = Paczynski(a, e, q)*(1.0 + q)/q
             Rout1 = Paczynski(a, e, 1.0/q)*(1.0 + q)
             OptDepth = AOptDepthC[int((job-100)/(len(Aq)*len(Aa))) % len(AOptDepthC)]
             for Eta in AEta:
              for WindAngle in AWindAngle:
               for Rin1 in ARin1:
                for Rin2 in ARin2:
                 # Reading from the database the 20 Principle component and reconstruct the BEL profile flux.
                 counting=counting+1
                 (Tscores) = struct.unpack('20f',f.read(80))							
                 Fnu=np.matmul(Tscores,Tloads.transpose())+Tmeandata

if counting==120*nC:
 print 'Finished WC2D.'
else:
 print 'Error WC2D.'


####### Working on the WC3D folder with 90 jobs, 3 tau values 0.1,1.0, 100.0 are stored in 3 binary files 'WC3D20pc2.bi' to 'WC3D20pc4.bi'#############

counting=0

for taujob in range(1, 4):
 with open('WC3D/WC3D20pc'+str(taujob+1)+'.bi','rb') as f:
  for job in range(100, 220):
   if int((job - 100)/(len(Aq)*len(Aa))) % len(AOptDepthC) == taujob:
    # Orientation variables: theta, phi, theta1, phi1, theta2, phi2 - affect inclination angles inc1, inc2.
    for theta in Atheta:
     for phi in AphiC:
      for theta1 in Atheta1:
       for phi1 in Aphi1:
        for theta2 in Atheta2:
         for phi2 in Aphi2:
          inc1 = np.arccos(np.sin(theta1)*np.cos(phi1)*np.sin(theta)*np.cos(phi) + np.sin(
                 theta1)*np.sin(phi1)*np.sin(theta)*np.sin(phi) + np.cos(theta1)*np.cos(theta))
          inc2 = np.arccos(np.sin(theta2)*np.cos(phi2)*np.sin(theta)*np.cos(phi) + np.sin(
                 theta2)*np.sin(phi2)*np.sin(theta)*np.sin(phi) + np.cos(theta2)*np.cos(theta))
          if (ilower > inc1 or inc1 > iupper) and (ilower > inc2 or inc2 > iupper):
           # Orbital phase or true anomaly: omega - affect the change of frame.
           for ophase in Aophase:
            # Circular or Elliptical orbits: e
            for e in AeC:
             # Semi-major axis a, and mass ratio q, are determined by the job number. The Semi-major axis a is in terms of total mass: a = a(M).
             q = Aq[(job-100) % len(Aq)]
             a = Aa[int((job-100)/len(Aq)) % len(Aa)]
             # The circumbinary disk size depends on the semi-major axis a.
             Rin3 = 2.0*a
             Rout3 = 3.0*a
             Rout2 = Paczynski(a, e, q)*(1.0 + q)/q
             Rout1 = Paczynski(a, e, 1.0/q)*(1.0 + q)
             OptDepth = AOptDepthC[int((job-100)/(len(Aq)*len(Aa))) % len(AOptDepthC)]
             for Eta in AEta:
              for WindAngle in AWindAngle:
               for Rin1 in ARin1:
                for Rin2 in ARin2:
                 # Reading from the database the 20 Principle component and reconstruct the BEL profile flux.
                 counting=counting+1
                 (Tscores) = struct.unpack('20f',f.read(80))							
                 Fnu=np.matmul(Tscores,Tloads.transpose())+Tmeandata

if counting==90*nC:
 print 'Finished WC3D.'
else:
 print 'Error WC3D.'


####### Working on the WE3D folder with 240 jobs, 2 tau values 0.0001, 1.0 are stored in 2 binary files 'WE3D20pc1.bi' and 'WE3D20pc2.bi' #########

counting=0

for taujob in range(0, 2):
 with open('WE3D/WE3D20pc'+str(taujob+1)+'.bi','rb') as f:
  for job in range(100, 460):
   if int((job - 100)/(len(Aq)*len(Aa))) % len(AOptDepthE) == taujob:
    # Orientation variables: theta, phi, theta1, phi1, theta2, phi2 - affect inclination angles inc1, inc2.
    for theta in Atheta:
     for phi in AphiE:
      for theta1 in Atheta1:
       for phi1 in Aphi1:
        for theta2 in Atheta2:
         for phi2 in Aphi2:
          inc1 = np.arccos(np.sin(theta1)*np.cos(phi1)*np.sin(theta)*np.cos(phi) + np.sin(
                 theta1)*np.sin(phi1)*np.sin(theta)*np.sin(phi) + np.cos(theta1)*np.cos(theta))
          inc2 = np.arccos(np.sin(theta2)*np.cos(phi2)*np.sin(theta)*np.cos(phi) + np.sin(
                 theta2)*np.sin(phi2)*np.sin(theta)*np.sin(phi) + np.cos(theta2)*np.cos(theta))
          if (ilower > inc1 or inc1 > iupper) and (ilower > inc2 or inc2 > iupper):
           # Orbital phase or true anomaly: omega - affect the change of frame.
           for ophase in Aophase:
            # Circular or Elliptical orbits: e
            for e in AeE:
             # Semi-major axis a, and mass ratio q, are determined by the job number. The Semi-major axis a is in terms of total mass: a = a(M).
             q = Aq[(job-100) % len(Aq)]
             a = Aa[int((job-100)/len(Aq)) % len(Aa)]
             # The circumbinary disk size depends on the semi-major axis a.
             Rin3 = 2.0*a
             Rout3 = 3.0*a
             Rout2 = Paczynski(a, e, q)*(1.0 + q)/q
             Rout1 = Paczynski(a, e, 1.0/q)*(1.0 + q)
             if taujob == 0:
              OptDepth = 0.0001
             elif taujob == 1:
              OptDepth = 1.0
             for Eta in AEta:
              for WindAngle in AWindAngle:
               Rin1 = ARin1[int((job-100)/(len(Aq)*len(Aa)*len(AOptDepthE))) % len(ARin1)]
               Rin2 = ARin2[int((job-100)/(len(Aq)*len(Aa)*len(AOptDepthE)*len(ARin1))) % len(ARin2)]
               # Reading from the database the 20 Principle component and reconstruct the BEL profile flux.
               counting=counting+1
               (Tscores) = struct.unpack('20f',f.read(80))
               Fnu=np.matmul(Tscores,Tloads.transpose())+Tmeandata

if counting==240*nE:
 print 'Finished WE3D.'
else:
 print 'Error WE3D.'

###################################################################################################################################################################
