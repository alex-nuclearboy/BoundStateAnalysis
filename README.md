# Search for eta-mesic Helium with the WASA-at-COSY

Meson–nucleus bound systems are considered to be very interesting objects in the modern nuclear and hadronic physics. 
The existence of the eta-nuclear bound states was first predicted by Haider and Liu [Q.Haider, L.C.Liu, Phys. Lett. B 172, 257 (1986)]. 
Experimental searches have been performed by several experiments, although, so far, there is no direct experimental confirmation of the existence of mesic nuclei.

The high statistics experiment devoted to the search for eta-mesic Hepium in the pd -> dp pi0 reaction was carried out with the WASA (Wide Angle Shower Apparatus) detection setup installed at the COSY accelerator in the Forschungszentrum Juelich. 
The WASA detector consisted of two main parts: the Forward Detector (FD) and Central Detector (CD) optimized for tagging the recoil particles and registering the meson decay products, respectively.

The measurement was performed changing the proton beam momentum very slowly and continuously around the eta-meson production threshold in each acceleration cycle from 1.426 to 1.635 GeV/c,  corresponding to the 3He-eta excess energy range Q from -70 to 30 MeV. 
The application of this so-called ramped beam technique allowed us to reduce the systematic uncertainties with respect to separate runs at fixed beam energies. 

Possible resonance-like structure below the eta-meson production threshold associated with the 3He-eta bound state was searched for via  measurement of the excitation function for the pd -> dp pi0 reaction. 

The events corresponding to formation of 3He-eta bound states were selected with appropriate conditions based on the Monte Carlo simulation of the pd -> (3He-eta)_(bound) -> dp pi0 reaction. 
The proton deuteron collision leads to the formation of a 3He nucleus bound with the eta-meson via strong interactions. Then, the eta-meson can be absorbed by one of the nucleons inside the helium exciting it to the N* (1535) nucleon resonance until the resonance decays into a proton pi0 pair, with the pion subsequently decaying into two photons.

Events selection for this process started with particles identification in the Central Detector. 
Protons were identified based on the energy deposited in the Scintillator Electromagnetic Calorimeter (SEC) combined with the energy loss in the Plastic Scintillator Barrel (PSB). 
The neutral pions pi0 were identified on the basis of the invariant mass of two photons originating from their decays and measured in the SEC.
Deuterons which were not directly registered in the experiment were identified via the missing mass technique.
The events corresponding to eta-mesic bound states were selected by applying cuts in the pi0-proton opening angle in the CM frame, in the missing mass as well as in the deuteron momentum distributions.

All files are distributed under the terms of the GNU General Public Licence Version 3.

## Software

### ROOT 
ROOT is an object-oriented framework for large scale data analyses based on C++.
It was developed for experiments at CERN, in which an impressive amount of data has to be processed. It was originally designed for particle physics data analysis.
ROOT has an extensive library with modules in physics, mathematics, statistics, histograms and graphics. It contains statistics tools, which can be used for data analysis.

### RootSorter
RootSorter is based on the ROOT data analysis framework. The framework is organized in different parts to handle different tasks like, decoding of data, calibration of the individual detector, track and energy reconstruction, particle identification etc.
The experimental data which come directly from the electronics are stored in a Hit Bank Raw. The simulated data are stored in Hit Bank MC. The experimental data goes through the calibration whereas the simulated data passes through filters before the reconstruction starts.
Then the track reconstruction, energy reconstruction and particle identification are done taking calibrated data or the filtered data. The reconstructed kinematic information are written to a ROOT formatted file for the physics analysis.

## Needed environment variables:

path where ROOT is installed

    ROOTSYS

path where RootSorter is instaled

    ROOTSORTERSYS

path where data from experiment are located

    RUNS_DATA

path where the simulation results are stored

    WMC_DATA

path to store analysis results for raw data

    OUTPUT_DATA

path to store analysis results for MC data

    OUTPUT_MC

## Running analysis

    ./runAnalysis-mc.sh <reaction>

for analysis of simulation results or

    ./runAnalysis-data.sh

for the analysis of raw data.
