# synthetic_spectra

The two databases and associated scripts are resources accompanying the following manuscripts: 

[1] Khai Nguyen, Tamara Bogdanovic 2016 "Emission Signatures from Sub-parsec Binary Supermassive Black Holes. I. Diagnostic Power of Broad Emission Lines", ApJ, 828, 68 (arXiv:1605.09389)

[2] Khai Nguyen, Tamara Bogdanovic, Jessie C. Runnoe, Michael Eracleous, Steinn Sigurdsson, and Todd Boroson 2019 "Emission Signatures from Sub-parsec Binary Supermassive Black Holes. II. Effect of Accretion Disk Wind", ApJ, 870, 16 (arXiv:1807.09782)

[3] Khai Nguyen, Tamara Bogdanovic, Jessie C. Runnoe, Michael Eracleous, Steinn Sigurdsson, and Todd Boroson 2020 "Emission Signatures from Sub-parsec Binary Supermassive Black Holes. III. Comparison of Models with Observations", ApJ, submitted (arXiv:1908.01799)

The resources are provided free of charge and with no technical support. We ask you to please cite the original manuscripts if you use the provided codes or data.

+ The folder "original_database" contains the entire database of synthetic broad emission-line (BEL) profiles calculated for 42,362,400 different SBHB configurations as described in paper [2]. (Size ~305 Gb)

+ The folder "reconstructed_database" contains reconstruction of the original BEL profile database using the Principal Component Analysis (PCA). Each BEL profile in the original database is reconstructed using the first 20 principal components. Each reconstructed profile is a close replica of the original profile but is not identical to it (because it is represented by a final series of principal components). The main advantage of the reconstructed database is that it can be written in a more compact form and is a factor of 100 smaller than the original database. (Size ~3.4 Gb)
