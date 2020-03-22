Documentation for the database of broad emission-line (BEL) profiles calculated for sub-parsec binary supermassive black holes.

================================================================
 
 Database: Broad Emission Line (BEL) Profiles 

 Code names: Navigation.f90, Navigation.py
 
 Code author: Khai Nguyen
 
 Date: 14 August 2019
 
 Languages: Fortran 90, Python 2
 
================================================================

I. Introduction:

This synthetic database and scripts are resources accompanying the following manuscripts: 

[1] Khai Nguyen, Tamara Bogdanovic 2016 "Emission Signatures from Sub-parsec Binary Supermassive Black Holes. I. Diagnostic Power of Broad Emission Lines", ApJ, 828, 68 (arXiv:1605.09389)

[2] Khai Nguyen, Tamara Bogdanovic, Jessie C. Runnoe, Michael Eracleous, Steinn Sigurdsson, and Todd Boroson 2019 "Emission Signatures from Sub-parsec Binary Supermassive Black Holes. II. Effect of Accretion Disk Wind", ApJ, 870, 16 (arXiv:1807.09782)

[3] Khai Nguyen, Tamara Bogdanovic, Jessie C. Runnoe, Michael Eracleous, Steinn Sigurdsson, and Todd Boroson 2020 "Emission Signatures from Sub-parsec Binary Supermassive Black Holes. III. Comparison of Models with Observations", ApJ, submitted (arXiv:1908.01799)

The resources are provided free of charge and with no technical support. We ask you to please cite the original manuscripts if you use the provided codes or data.


II. The database of BEL profiles:

+ The entire database contains BEL profiles of 42,362,400 different SBHB configurations as described in paper [2]. They are stored as binary files with extension '.bi' in 3 different folders:

 1) 'WC2D' contains 120 '.bi' files, each file contains  nC = 84,840 SBHB configurations with circular orbits where accretion disk wind occurs only on the two mini-disks.
 2) 'WC3D' contains  90 '.bi' files, each file contains  nC = 84,840 SBHB configurations with circular orbits where accretion disk wind occurs on all three disks.
 3) 'WE3D' contains 240 '.bi' files, each file contains nE = 102,275 SBHB configurations with eccentric orbits where accretion disk wind occurs on all three disks.

+ For each SBHB configuration, we provide 3 profiles that correspond to the BEL profiles from the broad line region (BLR) around the primary black hole, the secondary black hole, and the circumbinary disk. They can be added up to obtain the total BEL profile for the whole system.

+ Each BEL profile is defined over Nnu = 600 frequency bins. The profile flux values are stored as single precision float numbers. Hence, each '.bi' file in the WC2D, WC3D folder contains (nC * Nnu * 3) single precision float values.  And    each '.bi' file in the WE3D folder contains (nE * Nnu * 3) single precision float values. Since each single precision float value uses a memory of exactly 4 bytes, the total size of the whole database can be calculated as:
120 * 84,840 * 600 * 3 * 4 (bytes) + 90 * 84,840 * 600 * 3 * 4 (bytes) + 240 * 102,275 * 600 * 3 * 4 (bytes) =  305,009,280,000 (bytes) ~ 305 Gb . This corresponds to the size of the entire database.


III. The Navigation.f90 script: 

The script loops through all physical parameters of all the SBHB configurations in nested loop fashion.
For each SBHB configuration, the program reads and adds up contributions from the 3 BLRs, and then stores the BEL profile flux into the 1D array Fnu(Nnu) where Nnu=600 is the number of frequency bins. The wavelength of each frequency bin can be found from the 1D array lambda(Nnu). Note that the default emission line is H_beta with the rest wavelength 'lambdanot = 4860.09' Angstrom. For other emission lines, the rest wavelength, lambdanot, needs to be updated.


IV. The Navigation.py script:

This script similar to the above script but written in Python 2 language. 
