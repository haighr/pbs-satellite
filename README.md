## PBSsatellite: Spatio-temporal analysis of satellite data ##

<font color="red">&copy; Fisheries and Oceans Canada (2015-2023)</font>

**PBSsatellite** provides software designed to simplify the extraction and statistical analysis of gridded satellite data. This software extends the R Project for Statistical Computing, and it uses <a href="https://cran.r-project.org/package=PBSmapping">PBSmapping</a>, an existing R package, to aid in spatial analysis and the production of plots. The tools found in this package provide users with the functionality necessary to work with data from a variety of sources. Additionally, users are able to write their own data interpretation algorithms and provide them as arguments to some analysis functions within this package. The first version of <a href="https://doi.org/10.31542/r.gm:1451">PBSsatellite</a> was created by Nicholas Lefebvre and Nicholas Boers at MacEwan University in collaboration with Lyse Godbout and Rowan Haigh of Fisheries and Oceans Canada.

Three formats are widely used for exchanging and storing meteorological data: (a) Extensible Markup Language (XML), (b) Network Common Data Format (NetCDF), and (c) Hierarchical Data Format (HDF). XML is substantially more verbose than the other two formats, which leads to unnecessarily large files. Therefore, XML was not seriously considered for integration into **PBSsatellite**.

NetCDF was selected over HDF primarily due to the availability and quality of R packages for importing these files. At the time of writing, the Comprehensive R Archive Network (CRAN) did not host any packages explicitly for importing HDF files. The package **rgdal** provided bindings for the Geospatial Data Abstraction Library, which can import HDF files when appropriately configured. Unfortunately, this package was not configured to support HDF for versions built on Mac OS X and (reportedly) Windows. In contrast, three available NetCDF packages were hosted on CRAN: **ncdf**, **ncdf4**, and **RNetCDF**.

Of the three NetCDF packages, **ncdf4** was selected for the following advantages:<br>
• supports NetCDF versions 3 & 4;<br>
• supports offsetting into data ﬁles.<br>
While **ncdf4** requires an external package called **netcdf**, the offset functionality can effectively skip irrelevant data to reach desired data, making processing significantly faster.

As with any freely available product, there is no warranty or promise that **PBSsatellite** will perform adequately for all circumstances. Additionally, coding errors are possible, and users should contact the package maintainer if bugs are detected.

Maintainer: <a href="mailto:rowan.haigh@dfo-mpo.gc.ca">Rowan Haigh</a>

<p align="right"><img src="DFOlogo_small.jpg" alt="DFO logo" style="height:30px;"></p> 





