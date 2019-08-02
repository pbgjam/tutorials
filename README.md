# NEON-GJAM
NASA NEON data project

#### Personnel
Dr. Jim Clark (Co-PI)

Dr. Jennifer Swenson (Co-PI)

John Fay (Data analyst)

Amanda Schwantes (Data Analyst, Post doc)

Brad Tomasek (Data analyst)

Christ Kilner (Phd Student)

#### Overview
This NSF funded project integrates key remote sensing variables with continental scale ecological data to provide broadly accessible ecological forecasts to determine which species and communities are most vulnerable to climate change and to forecast responses of entire communities. 

The scripts hosted here consist of the following sections:
 - **NEON**: Scripts to extract and synthesize species observation records from the [NEON repository](http://data.neonscience.org/data-api). 
 - **DAAC**: Scripts to extract remotely sensed data for specific locations from the [NASA LP-DAAC](https://lpdaac.usgs.gov/).
 - **GJAM**: Scripts to use the Generalized Joint Attribute Model ([GJAM](https://cran.r-project.org/web/packages/gjam/vignettes/gjamVignette.html)) package using the extracted data above. 
 - **Interface**: Scripts to assemble the web interface allowing users to visualize impacts under different climate scenarios.
