nustar_stray_light
==================

- Lightweight version of the stray light checking code for NuSTAR. Heavily indebted to the "nuplan.pro" fully-featured code written by Roman Krivonos and Kaya Mori (among others who contibuted to earlier versions of this code) which allows the user to generate pixel masks to ingore the stray light as well as many other things. This repo is intended to allow for a quick check of the dependence of the stray light source on the PA angle.
- Updated 2014/08/15 by Brian Grefenstette (bwgref@srl.caltech.edu)

- This code runs in two modes: "Single PA" mode where IDL will display the stray light pattern (SLP) and "PA scan" mode where the PA is scanned from 0 to 360 degrees in 5 degree steps. A summary PDF is produced (scan.pdf) that shows the SLP for each PA angle and a summary data file (pa_scan.dat) file gives some statistcs of the stray light coverage for FPMA/B and Det 0 on FPMA/B for each stray light step as well as the list of possible sources of stray light for the target.

- To run:
- It is assumed that you have the nustar-lib git repo installed or you've got the AstroLib and Coyote librarie installed. If the latter is true, do NOT use the setup_stray_light.sh initialization script. If you ARE using the nustar-lib it is recommended that you edit the setup_stray_light.sh where indicated and use the custom IDL startup.pro script provided with this repo.

- Start IDL:

- For single PA mode, you need to provide the RA and Dec in degrees, as well as the PA angle. The coordinates of M83 are provided as en example in the repo, but you could also set these coordinates yourself for the desired target. If these values are not provided then the script will prompt you for RA, DEC, and PA.
- IDL> ra = 204.254 & dec = -29.8654 & pa = 10
- Now call nustar_stray light:
- IDL> nustar_stray_light, ra, dec, pa =pa
- Some diagnosic information should be displayed in the terminal window, while IDL will display the stray light patterns in the FoV on the current graphics device (e.g. an "X" window).

- For PA scan mode, set the RA and Dec as above, but enable the /do_scan keyword:
- IDL> nustar_stray_light, ra, dec, /do_scan

- Note that the PA scan can take a little bit of time as it depends on the number of stray light sources to be rendered. Also note that the output files (scan.pdf and pa_scan.dat) will be overwritten every time the script is run, so please make a copy of a particular run if you want to save the results.

- There are numberous keywords that can affec the code; please look into the nustar_stray_light.pro script for details.
