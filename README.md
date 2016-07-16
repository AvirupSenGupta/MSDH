## MERRA Spatial Downscaling for Hydrology (MSDH)
* Description: an R based tool to downscale globally available Modern Era Retrospective-Analysis for Research and Applications (MERRA) temperature, precipitation, wind speed, relative humidity, incoming shortwave and longwave radiation data to overcome data scarcity for hydrological analysis and modeling.

### Software Information ###

* Version: 1.0
* Developers: Avirup Sen Gupta, and David Tarboton
* Year first available:  2013
* Source Code: https://bitbucket.org/AvirupSenGupta/msdh.usu/src.
* License: GNU General Public License version 3, http://www.gnu.org/licenses/gpl-3.0.html 
* Program language: R
* Hardware:  PC running Microsoft Windows
* Software dependencies:  
  - netCDF Operator (NCO), available at http://nco.sourceforge.net/
  - Climate Data Operators (CDO), available at https://code.zmaw.de/projects/cdo/files.  
  - GTK+, available at http://www.gtk.org/
  - R, available at http://cran.r-project.org/bin/windows/base/



### Using the software ###

Running "MERRA Spatial Downscaling for Hydrology (MSDH)" in windows computer

1. download and install GTK2+ from: http://sourceforge.net/projects/gtk-win/files/latest/download

2. The bin directory of GTK2+ installation must be set as path in environement variable
Please follow the following procedure to do so.
Click on -> Control Panel -> System -> Advanced
Click on Environment Variables, under System Variables, find PATH, and click on it.
In the Edit windows, modify PATH by adding the location of the class to the value for PATH. 
If you do not have the item PATH, you may select to add a new variable and add PATH as the name 
and the location of the class as the value.

3. Download and install R from http://cran.r-project.org/bin/windows/base, if you have not done that already.

4. The bin directory of R installation must be set as Environmental Variable path by similar procedure described above.

5. Download and install NCO from http://nco.sourceforge.net/ and set its installaton directory as Environmental Variable path.

6. Download and install CDO from https://code.zmaw.de/projects/cdo/files and set its installaton directory as Environmental Variable path.

7. Make sure RunMSDH.exe and MSDH_USU.R are under the same folder.

8. Click on RunMSDH.exe.


### Contact Information ###

If you wish to use or incorporate this program (or parts of it) into 
other software that does not meet the GNU General Public License 
conditions contact the author to request permission.

Avirup Sen gupta  
Utah State University  
4110 Old Main Hill  
Logan, UT 84322-4110  
USA   
email:  avirup.sengupta@aggiemail.usu.edu
