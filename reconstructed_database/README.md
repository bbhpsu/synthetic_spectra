================================================================
 
 Database: Reconstructed database using the first 20 principal components
 
 Link to the database: https://www.dropbox.com/sh/ww9dsy28i20m6al/AABaR8ztRBGoFvEzPEpclpxma?dl=0

 Database size: ~3.4Gb

 Code name: Reconstruction.f90, Reconstruction.py
 
 Code author: Khai Nguyen
 
 Date: 14 August 2019
 
 Language: Fortran 90, Python 2
 
================================================================

I. Introduction:

This database is a reconstruction of the BEL profile database using the Principal Component Analysis (PCA). Each BEL profile in the original database is reconstructed using the first 20 principal components (PCs). Each reconstructed profile is a close replica of the original profile but is not identical to it (as it is represented by a final series of PCs). The main advantage of the reconstructed database is that it can be written in a more compact form and is a factor of 100 smaller than the original. Note that, unlike the original database of modeled profiles, the PCA reconstructed database does NOT contain information about the contribution from individual  BELs associated with the primary, secondary mini-disk, or the circumbinary disk, but only contains information about the total, composite BEL profile of the whole system. The database and scripts here are resources accompanying the following manuscripts: 

[1] Khai Nguyen, Tamara Bogdanovic 2016 "Emission Signatures from Sub-parsec Binary Supermassive Black Holes. I. Diagnostic Power of Broad Emission Lines", ApJ, 828, 68 (arXiv:1605.09389)

[2] Khai Nguyen, Tamara Bogdanovic, Jessie C. Runnoe, Michael Eracleous, Steinn Sigurdsson, and Todd Boroson 2019 "Emission Signatures from Sub-parsec Binary Supermassive Black Holes. II. Effect of Accretion Disk Wind", ApJ, 870, 16 (arXiv:1807.09782)

[3] Khai Nguyen, Tamara Bogdanovic, Jessie C. Runnoe, Michael Eracleous, Steinn Sigurdsson, and Todd Boroson 2020 "Emission Signatures from Sub-parsec Binary Supermassive Black Holes. III. Comparison of Models with Observations", ApJ, submitted (arXiv:1908.01799)

The resources are provided free of charge and with no technical support. We ask you to please cite the original manuscripts if you use the provided codes or data.


II. The PCA reconstructed database:

+ The reconstructed database contains BEL profiles of 42,362,400 different SBHB configurations as described in paper [2]. They are stored as binary files with an extension '.bi' in 3 different folders:

 1) 'WC2D' contains 4 '.bi' files, each file contains  30 *  84,840 SBHB configurations with circular orbits where accretion disk wind occurs only on the two mini-disks.
 2) 'WC3D' contains 3 '.bi' files, each file contains  30 *  84,840 SBHB configurations with circular orbits where accretion disk wind occurs on all three disks.
 3) 'WE3D' contains 2 '.bi' files, each file contains 120 * 102,275 SBHB configurations with eccentric orbits where accretion disk wind occurs on all three disks.

+ For each SBHB configuration, only the total BEL profile is provided.

+ Each BEL profile contains npc = 20 single precision components.

Hence, each '.bi' file in the WC2D or WC3D folder contains '30  *  84,840 * 20' single precision float values. 
And each '.bi' file in the WE3D folder contains '120 * 102,275 * 20' single precision float values. Since each single precision float value uses a memory of exactly 4 bytes, the total size of the whole PCA reconstructed database can be calculated as:
4 * 30 * 84,840 * 20 * 4 (bytes) + 3 * 30 * 84,840 * 20 * 4 (bytes)  + 2 * 120 * 102,275 * 20 * 4 (bytes)  =  3,388,992,000 (bytes) ~ 3.4Gb . This corresponds to the size of the entire reconstructed database.



III. The included ASCII files: 

Two text files are included as they are required to reconstruct the profiles from the 20 PCs.

 1) Tmeandata.txt 
This text file contains the average profile of all BEL profiles in the main database. Hence it's a list of 600 numbers that correspond to the flux of the average profile in 600 frequency bins.

 2) Tloads20.txt
This text file contains the first 20 principal eigenprofiles. Hence it's an array of 600 rows x 20 columns.


IV. The Reconstruction.f90 script:

+ The script reads the average profile from the file 'Tmeandata.txt' and stores it in a 1D array Tmeandata(600)
+ The script reads the first 20 eigenprofiles from the file 'Tloads20.txt' and stores it in a 2D array Tloads(20x600)
+ The script loops through all physical parameters of all the SBHB configurations in nested loop fashion, similar to the Navigation.f90 script.
+ For each SBHB configuration, the script reads the first 20 PCs of each BEL profile from the PCA reduction database, and stores them in the 1D array Tscores(20).
+ The script then reconstructs the BEL profile and stores it in the 1D array Fnu(600) by using the following relation: Fnu = (Tscores * Tloads) + Tmeandata

The wavelength of each frequency bin can be found from the array lambda(Nnu). Note that the default emission line is H_beta with the rest wavelength 'lambdanot = 4860.09' Angstrom. For other emission lines, the rest wavelength, lambdanot, needs to be updated.


