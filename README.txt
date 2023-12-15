The current folder contains data and code to perform the analyses described in the accompanying manuscript and supplementary information. 

This text file describes each of the subfolders. 

- data: contains raw and cleaned incidence data. Files whose name begins with 'h' or 'c' refer to hospitalisations and reported cases, respectively. This identifier is followed by the name of the Brazilian state and the period it refers to. Other files include serotype information and population counts. Files 'mca_c_fit.csv' and 'mca_h_fit.csv' contain the mean age of dengue cases (c) and hospitalisations (h) by year and state in Brazil; these files are necessary for gaussian process analysis. 'uf_locations.csv' contains spatial coordinates of representative points for each state in Brazil. The folder 3vFSA contains weekly incidence of DENV and ZIKV in Feira de Santana, Bahia (Fig 1 in the main manuscript).

- shapefiles: contains shapefiles of Brazil for plotting purposes.

- r_scripts: contains an R script (glm_age_shift.R) that performs a Bayesian GLM analysis of age shifts of DENV incidence across Brazilian states. Relevant data is contained in the file data_fit.csv in the same folder. Required R packages are: rstanarm, bayesplot and ggplot2.

- pyflavi: contains C++ code to simulate DENV and ZIKV transmission. The subfolder pyflavi/pyflavi contains code to create a python module (pyflavi) wrapping C++ code. pyflavi must be compiled before it can be imported and requires pybind11. Pybind11 is hosted at https://github.com/pybind/pybind11 and additional details on installation can be found at https://pybind11.readthedocs.io/en/stable/installing.html. To compile pyflavi, open the terminal and navigate to the inner pyflavi folder (where a file named 'makefile' is located). Once there, type 'make' in the command line to compile pyflavi. Compilation requires a working C++11 compiler, like g++ for Ubuntu or clang for macOS.
We were able to compile pyflavi on a macOS and an Ubuntu machine. No test has been carried out on Windows machines.

- python_scripts: contains four Jupyter notebooks:
    - plot_data.ipynb: contains code to generate all figures that show data. Some intermediate files are also found in the same folder.
    - gaussian_process_trend_analysis.ipynb: contains code to analyse trends of mean age of dengue cases using gaussian process regression.
    - probability_of_infection_children.ipynb: contains code to calculate the probability of a child being infected with dengue virus in 2018 or 2019 using estimates of force of infection in the Northeast in Brazil.  
    - running_simulations.ipynb: contains minimal examples to run simulations of ZIKV and DENV using different initialisation schemes for each pathogen (e.g. direct simulation vs serology-based initialisation of DENV dynamics). 
Shown examples also demonstrate how to manipulate simulated incidence to obtain reported incidence, measure mean case age and age shifts. Dependencies include numpy, scipy, pandas, geopandas and pyflavi (which has to be compiled before it can be imported, as described above).