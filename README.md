# Bessel-droplet 

Codes were developed in MATLAB R2018a with Windows 10 installed on a Dell workstation (Xeon CPU E5-2667, RAM 128 GB). The main codes to generate phase patterns on SLM for Bessel-droplet foci are “Bessel_droplet_SIM.m” and “Bessel_droplet_SIM_dPhase.m”. The codes to calculate the aberrated Bessel and Bessel-droplet foci are “AstiBesselPSF.m” and “AstiBesselDropletPSF.m”. Related functions are included in the accompanying folder ‘Functions’. The purposes and usages of the Supplementary Codes are listed below with the code procedures describing the computation flow.

Bessel_droplet_SIM.m
Purpose: This code is used to (1) Find the optimized illumination annulus diameters for Bessel-droplet foci; (2) Generate the phase mask for the optimized Bessel-droplet foci of different numerical apertures (e.g., NA=0.4, 0.5, 0.6, 0.7), and (3) Generate optional 3D PSF for Bessel-droplet foci.
Usage: Open the "Bessel_droplet_SIM.m" file in the folder to run the main function. Specify the parameter used for simulation: 
flag_ScanR2: Set flag_ScanR2=1 to vary the R2 value from 0.2*R1 to R1, where R1 defines the NA. The program will calculate the side-ring ratio between the peak of the most prominent side ring and the peak of the central peak of the lateral PSF, and save the results as .tif images in the “Results” folder. The optimized R2 is found when the side-ring ratio is minimum. Then set the flag_ScanR2=0 and flag_saveTemResults=1 to use the optimized R2 for simulation.
flag_3D: Set flag_3D=1 to calculate and save 3D PSF (a very time-consuming calculation). Set flag_3D=0 when scanning R2, set flag_3D=1 when simulating the optimized R2.
flag_mask: Set flag_mask ='A','B','DR' to simulate the cases for outer-annulus-only, inner-annulus-only, or double-annulus illumination at the objective back pupil plane.
Depending on the workstation performance and the simulation settings, this computation takes about 5~10 minutes.
Code procedures:
Step 1. Start the main program “Bessel_droplet_SIM.m”
Specify the system parameters (SLM pixel size, pixel number, wavelength, magnification, et al)
Specify the mask parameters used for Bessel-droplet foci (inner diameter, outer diameter, annulus thickness)
Specify the simulation purpose: scanning R2 or only computing at the optimized R2. Then based on this setting, start iterative calculations or a single calculation.
Step 2. Initialize an optical field (field_Pupil_ring) for the Bessel SLM located at the pupil plane
Assign Gaussian amplitude profile.
Assign the phase values (e.g, 0 for both the outer ring and inner ring).
Construct the desired optical complex field at the objective pupil plane.
Step 3. Fourier transform field_Pupil_ring to the SLM conjugated to the objective focal plane
Do a one-dimensional Fourier transform with “Fourier_CircularLens.m”
Display the amplitude and phase patterns on the SLM which is conjugated to the objective focal plane.
Step 4. Validation: Propagate the phase patterns on the SLM to the mask (the mask is conjugated to the objective pupil plane) 
Propagate the phase patterns on the SLM with a Gaussian amplitude profile to the mask.
Plot the resulting amplitude and phase profile at the mask, with red lines indicate the transmissive area of the mask.
Step 5. (optional) Calculate the 3D PSF using Richard Wolf PSF model for high NA objective.
Further propagate the complex optical field at the mask to the back pupil plane of the objective and calculate the PSF for Bessel-droplet foci with a high NA objective.

Bessel_droplet_SIM_dPhase.m
Purpose: This code is used to (1) Generate the phase masks to interferometrically scan the Bessel-droplet foci in Z axis (e.g., introduce a phase difference of π between the two annular illuminations at the pupil plane). (2) Generate optional 3D PSF for Bessel-droplet foci interferometrically scanned in Z axis.
Usage: Open the " Bessel_droplet_SIM_dPhase.m" file in the folder to run the main function. Specify the parameter used for simulation: 
Used the optimized R2 from the “Bessel_droplet_SIM.m”for the corresponding R1.
dPhase: Set dPhase=0 to generate two illumination annuli with the same phase value; set dPhase=π to introduce a phase difference of π between the two annuli;
flag_3D: Set flag_3D=1 to calculate and save the 3D PSF. Set flag_3D=0 without calculating the 3D PSF..
flag_mask: Set flag_mask ='DR'.
Depending on the workstation performance and the simulation settings, this computation takes about 5~10 minutes.
Code procedures:
Step 1. Start the main program “Bessel_droplet_SIM_dPhase.m”
Specify the system parameters (SLM pixel size, pixel number, wavelength, magnification, et al)
Specify the mask parameters used for Bessel-droplet foci (inner diameter, outer diameter, annulus thickness). Used the optimized R2 from the “Bessel_droplet_SIM.m”.
Set the phase offset value between the two annuli (e.g., dPhase = pi).

Step 2. Initialize an optical field (field_Pupil_ring) for the Bessel SLM located at the pupil plane
Assign Gaussian amplitude profile.
Assign the phase values (e.g, 0 for the outer ring and π for the inner ring).
Construct the desired optical complex field at the objective pupil plane.
Step 3. Fourier transform field_Pupil_ring to the SLM conjugated to the objective focal plane
Do a one-dimensional Fourier transform with “Fourier_CircularLens.m”
Display the amplitude and phase patterns on the SLM which is conjugated to the objective focal plane.
Step 4. Validation: Propagate the phase patterns on the SLM to the mask (the mask is conjugated to the objective pupil plane) 
Propagate the phase patterns on the SLM with a Gaussian amplitude profile to the mask.
Plot the resulting amplitude and phase profile at the mask, with red lines indicate the transmissive area of the mask.
Step 5. (optional) Calculate the 3D PSF using Richard Wolf PSF model for high NA objective.
Further propagate the complex optical field at the mask to the back pupil plane of the objective and calculate the PSF for Bessel-droplet foci with a high NA objective.

AstiBesselPSF.m and AstiBesselDropletPSF.m
Purpose: These codes calculate aberrated Bessel and Bessel-droplet foci (e.g., NA=0.7, aberrated by astigmatism) at the focal plane of an objective (Olympus 25x, NA 1.05). 
Usage: Open the " AstiBesselPSF.m " or "AstiBesselDropletPSF.m" files in the folder. Run the program and select the aberration phase pattern you want to apply at the back pupil of the objective. One example aberration phase pattern as astigmatism can be loaded from the folder. The lateral profiles of the aberrated foci are save in the result folders generated by the program.
