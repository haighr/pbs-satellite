\name{download.ncdfData}
\alias{download.ncdfData}
\title{
  Download NetCDF Satellite Data Files
}
\description{
  Download NetCDF data files from specified URL. 
  Currently, two sources are recognised:\cr
  \href{https://www.ncdc.noaa.gov/data-access/marineocean-data/extended-reconstructed-sea-surface-temperature-ersst-v4}{Extended Reconstructed Sea Surface Temperature (ERSST) v4}\cr
  \href{https://oceancolor.gsfc.nasa.gov/about/}{NASA's Ocean Biology Processing Group}
}
\usage{
download.ncdfData(url, save.dir="./rawdata",
   years, months, days, res, pname, period, sensor, curlbin,
   need.seq=FALSE, overwrite=FALSE, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{url}{\code{character} -- string specifying URL (uniform resource locator) for satellite data; current choices:\cr
    ERSST: \code{"https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v4/netcdf"}
    Ocean Color Data (OCD): \code{"https://oceandata.sci.gsfc.nasa.gov/"} }
  \item{save.dir}{\code{character} -- string specifying a directory on user's computer
    (can be relative to current working directory) to save data files.}
  \item{years}{\code{numeric} -- vector of years for data files (relevant for ERSST, OCD).}
  \item{months}{\code{numeric} -- vector of months for data files (ERSST).}
  \item{days}{\code{numeric} -- vector of Julian days (1 to 365/366) for data files (OCD).}
  \item{res}{\code{character} -- resolution of data, e.g. \code{"4km"}, \code{"9km"} (OCD).}
  \item{pname}{\code{character} -- parameter/product name specifying which data, e.g. \code{"CHL_chlor_a"} (OCD).}
  \item{period}{\code{character} -- periodicity of data, e.g. \code{"DAY"} (OCD).}
  \item{sensor}{\code{character} -- name of sensor/satellite, e.g. \code{"MODIS-Aqua"} (OCD).}
  \item{curlbin}{\code{character} -- if available, path to \code{bin} directory of user's \code{curl} installation.}
  \item{need.seq}{\code{logical} -- if \code{TRUE}, enforce sequential \code{years}, \code{months}, and \code{days} (ERSST, OCD).}
  \item{overwrite}{\code{logical} -- if \code{TRUE}, overwrite existing files in \code{save.dir}.}
  \item{verbose}{\code{logical} -- if \code{TRUE}, display some information about the data source.}
}
\details{
  This function evolved from Michael Malick's \code{'sst_download'} function for ERSST data.
  See Malick's \href{https://github.com/michaelmalick/sst-ersst-v4}{GitHub} site for more details.

  A more generalised set of data exist from NASA's \href{https://oceancolor.gsfc.nasa.gov/about/}{OBPG website},
  specifically their \href{https://oceandata.sci.gsfc.nasa.gov/}{Direct Data Access} portal.
  The ocean color data (OCD) by sensor (or satellite) are offered in five levels of formatting:
  \itemize{
  \item \bold{Level 0} -- data are unprocessed instrument/payload data at full resolution.
        Any artifacts of the communication (e.g. synchronization frames, communication headers)
        of these data from the spacecraft to the ground station have been removed.
        These data are the most raw format available, and are only provided for a few of the missions that OBPG supports.
  \item \bold{Level 1A} -- data are reconstructed, unprocessed instrument data at full resolution,
        time-referenced and annotated with ancillary information including radiometric and geometric
        calibration coefficients and georeferencing parameters (e.g. platform ephemeris data) computed
        and appended but not applied to the Level 0 data.
        It is this last point that makes Level 1A the preferred archival data level.
        If the sensor calibration changes, the data \emph{do not need to be reacquired} -- no mean feat where large data sets are concerned.
        The same is not true for Level 1B data (described below) which must be replaced every time the sensor calibration changes.
  \item \bold{Level 1B} -- data are Level 1A data that have had instrument/radiometric calibrations applied.
  \item \bold{Level 2} -- data consist of derived (by the \code{l2gen} program) geophysical variables 
        at the same resolution as the source Level 1 data. 
        These variables are grouped into a few product suites (e.g. OC [ocean color], SST, and SST4 for MODIS).
  \item \bold{Level 3} -- data are derived geophysical variables that have been aggregated/projected onto a
        well-defined spatial grid over a well-defined time period. OBPG archives two types of Level 3 data:
    \itemize{
    \item \bold{Binned} -- each Level 3 binned data product consists of the accumulated data for all L2 products
          in a product suite, for the specified instrument and resolution, corresponding to a period of time 
          (e.g. daily, 8 days, monthly, etc.) and stored in a global, nearly equal-area, integerized sinusoidal grid.
    \item \bold{Mapped} -- the Level 3 Standard Mapped Image (SMI) products are created from the corresponding Level 3 binned products. 
          Each SMI file contains a Plate Carre\ifelse{latex}{\out{\'{e}}}{\ifelse{html}{\out{&eacute;}}{e}}, 
          pixel-registered grid of floating-point values (or scaled integer representations of the values) for a single geophysical parameter. 
          A color look-up table is also provided in each file that may be used to generate an image from the data. 
    }
  \item \bold{Level 4} -- data are model output or results from analyses of lower level data
        (e.g., variables derived from multiple measurements).
        Ocean Primary Productivity is a good example of a Level 4 product.
  }
  \bold{IMPORTANT} -- This function can only access the L3SMI products.

  The OCD products all include combinations of the following derived parameters (\code{'pname'} in parentheses):
  \itemize{
  \item \bold{chl-a} -- chlorophyll-A concentration (\code{CHL_chlor_a})
  \item \bold{Rrs}   -- remote sensing reflectance (\code{RRS_Rrs_nnn})\cr
        where wavelength 'nnn' is \ifelse{latex}{\out{$\in$}}{\ifelse{html}{\out{&#8712;}}{element of}}
        \{412,443,469,488,531,547,555,645,667,678\}
  \item \bold{AOT}   -- aerosol optical thickness, \ifelse{latex}{\out{$\tau$}}{\ifelse{html}{\out{&#120591;}}{tau}},
        in one NIR or red band (\code{RRS_aot_869})
  \item \bold{A}     -- angstrom coefficient (\code{RRS_angstrom})
  \item \bold{K490}  -- diffuse attenuation coefficient at 490 nm (\code{KD490_Kd_490})
  \item \bold{PIC}   -- calcite concentration or particulate inorganic carbon (\code{PIC_pic})
  \item \bold{POC}   -- particulate organic carbon (\code{POC_poc})
  \item \bold{PAR}   -- photosynthetically available radiation (\code{PAR_par})
  \item \bold{FLH}   -- fluorescence line height (\code{FLH_ipar}, \code{FLH_nflh})
  \item \bold{IOPs}  -- inherent optical properties (\code{IOP}s, numerous), which include
        absorption and backscattering coefficients in the visible bands.
  \item \bold{SST}   -- sea surface temperature (\code{NSST_sst}, \code{SST4_sst}, \code{SST_sst})
  }

  [NASA's \href{https://oceancolor.gsfc.nasa.gov/browse_help/l3/}{Level 3 Browse}]\cr
  Level-3 products are produced for several different compositing periods. 
  In general, longer time periods fill in more of the naturally occurring data gaps
  (caused by, for example, clouds, sun glint, inter-orbits gaps, ice, low light, etc.)
  at the expense of short-lived features which tend to get smoothed out in longer-period composites.

  The following \code{'period'} composites are available from OCD's Direct Data Access:
  \itemize{
  \item \code{DAY} -- daily snapshot: one data day's worth of measurements.
  \item \code{8D}  -- 8-day (weekly): eight consecutive days' worth of data with predetermined start and stop days;
        the first 8-day period of each year always begins with January 1, the second with January 9, the third with January 17, etc.
        The final \dQuote{8-day} composite of each year comprises only five days in non-leap years (27-31 December)
        or six days in leap years (26-31 December).
  \item \code{MO}  -- monthly: one calendar month's worth of data; this will be 28, 29, 30, or 31 days' worth of data depending on the month.
  \item \code{R32} -- rolling 32-day: 32 days' worth of data with each subsequent composite shifted 8 days further along the calendar;
        the start dates of the composites are the same as those of the corresponding 8-day composites;
        those that start in December wrap into the subsequent year. 
        These composites are intended to be used as frames in animations;
        the fact that neighboring frames share 24 of the 32 days that they are composed of makes for smoother animations.
  \item \code{YR}  -- annual: one calendar year's worth of data.
  }

  Sensors with L3SMI data are:
  \itemize{
  \item \bold{Aquarius} -- The Aquarius instrument was the first satellite-based ocean salinity sensor.
        Aquarius combined a passive radiometer with an active pulsed radar scatterometer,
        which share the same optics (feed horns and antenna). 
        The radiometer was built at GSFC, and the scatterometer and overall integration were the responsibility of the 
        NASA Jet Propulsion Laboratory (JPL). 
        \tabular{ll}{
        Satellite:            \tab SAC-D (Argentina)\cr
        Launch Date:          \tab Jun 10, 2011\cr
        Data Collection:      \tab Aug 25, 2011 to Jun 7, 2015\cr
        Orbit Type:           \tab Polar Sun-synchronous\cr
        Altitude (nominal):   \tab 657 km\cr
        Orbit Period:         \tab 98 minutes (maintained 7-day repeat cycle)\cr
        Equator crossing time/direction: \tab 6:00 PM ascending
        }
  \item \bold{CZCS} -- CZCS was the first satellite-based ocean color sensor, and was flown on Nimbus-7 as a proof-of-concept experiment.
        The onboard data storage was limited, and data collection was scheduled over designated locations.
        The instrument was officially declared dead on December 18, 1986.
        \tabular{ll}{
        Satellite:            \tab Nimbus-7 (USA)\cr
        Launch Date:          \tab Oct 25, 1978\cr
        Data Collection:      \tab Oct 30, 1978 to Jun 22, 1986\cr
        Orbit Type:           \tab Polar Sun-synchronous\cr
        Altitude (nominal):   \tab 955 km\cr
        Orbit Period:         \tab ~104 minutes\cr
        Equator crossing time/direction: \tab 12 noon ascending
        }
  \item \bold{MERIS} -- The MERIS instrument was a programmable, medium-spectral resolution, imaging spectrometer operating
        in the solar reflective spectral range. Fifteen spectral bands can be selected by ground command, although the 
        selected band suite remained fixed for the duration of the mission.
        \tabular{ll}{
        Satellite:            \tab Envisat (Europe)\cr
        Launch Date:          \tab March, 2002\cr
        Data Collection:      \tab Apr 29, 2002 to Apr 8, 2012\cr
        Orbit Type:           \tab Polar Sun-synchronous\cr
        Altitude (nominal):   \tab 800 km\cr
        Orbit Period:         \tab 101 minutes\cr
        Equator crossing time/direction: \tab 12 noon ascending
        }
  \item \bold{MODIS-Terra | MODIS-Aqua} -- Terra and Aqua were the first two (of three, the third being Aura)
        large Earth Observing System (EOS) satellites launched and operated by NASA. 
        Both EOS missions are multi-sensor and multi-disciplinary, and continue to operate without significant problems.
        \tabular{ll}{
        Satellite:            \tab Terra | Aqua (USA)\cr
        Launch Date:          \tab Dec 18, 1999 | May 4, 2002\cr
        Data Collection:      \tab Feb 24, 2000 to present | Jul 4, 2002 to present\cr
        Orbit Type:           \tab Polar Sun-synchronous\cr
        Altitude (nominal):   \tab 700 km\cr
        Orbit Period:         \tab 98.9 minutes (maintains 16-day repeat cycle)\cr
        Equator crossing time/direction: \tab 10:30 AM descending | 1:30 PM ascending
        }
  \item \bold{OCTS} -- ADEOS was launched and operated by the Japanese space agency (then NASDA, now JAXA).
        It was designed for a multi-year mission but suffered a fatal spacecraft control failure after nine months of operations (Jun 30, 1997).
        OCTS included bands for both ocean color (visible/NIR) and thermal IR measurements.
        \tabular{ll}{
        Satellite:            \tab ADEOS (Japan)\cr
        Launch Date:          \tab Aug 17, 1996\cr
        Data Collection:      \tab Nov 1, 1996 to Jun 29, 1997\cr
        Orbit Type:           \tab Polar Sun-synchronous\cr
        Altitude (nominal):   \tab 800 km\cr
        Orbit Period:         \tab 101 minutes\cr
        Equator crossing time/direction: \tab 10:30 AM descending
        }
  \item \bold{SeaWiFS} -- SeaWiFS was the first sensor designed to support global ocean color data collection, 
        although it was preceded into orbit by OCTS.
        Global Area Coverage (GAC) data, which were recorded globally, were subsampled at 4x4 from the full-resolution sensor data,
        and also truncated at the start and end of the scan. 
        Local Area Coverage (LAC) data were collected at full resolution, and recorded in limited quantities over selected targets;
        full-resolution were also direct-broadcast to ground stations in High Resolution Picture Transmission (HRPT) format.
        The sensor continued to operate normally with no apparent problem until the end of the mission.
        \tabular{ll}{
        Satellite:            \tab OrbView-2 (USA)\cr
        Launch Date:          \tab Aug 1, 1997\cr
        Data Collection:      \tab Sep 18, 1997 to Dec 11, 2010\cr
        Orbit Type:           \tab Polar Sun-synchronous\cr
        Altitude (nominal):   \tab 705 km at launch (decayed to 690km, raised to 782km)\cr
        Orbit Period:         \tab ~99 minutes\cr
        Equator crossing time/direction: \tab 12 noon descending (at launch)
        }
  \item \bold{VIIRS-SNPP | VIIRS-JPSS1} -- The SNPP satellite was built and launched by NASA as a risk reduction mission for JPSS.
        SNPP carries four other sensors besides VIIRS. The first JPSS satellite, JPSS-1, has an identical sensor suite.
        VIIRS will be flown on every JPSS satellite as well as on SNPP. 
        It is conceptually similar to MODIS, with comparable resolution, spectral coverage and support for multidisciplinary data product generation.
        \tabular{ll}{
        Satellite:            \tab SNPP | JPSS1 (USA)\cr
        Launch Date:          \tab Oct 28, 2011 | Nov 10, 2017\cr
        Data Collection:      \tab Jan 2, 2012 to present\cr
        Orbit Type:           \tab Polar Sun-synchronous\cr
        Altitude (nominal):   \tab 830 km\cr
        Orbit Period:         \tab 101.5 minutes (maintains 16-day repeat cycle)\cr
        Equator crossing time/direction: \tab 1:30 PM ascending
        }
  }
}
\value{
  Downloads NetCDF files (\code{'.nc'}) from a website to the user's computer.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2018-10-31}
}
\note{
  When specifying the Ocean Color Data source, the R code can only handle \code{period="DAY"} at present.
}
\seealso{
  \code{\link[PBSsatellite]{convert.ncdfData}},
  \code{\link[PBSsatellite]{create.ncdfData}},
  \code{\link[PBSsatellite]{merge.ncdfData}},
  \code{\link[PBSsatellite]{plot.ncdfData}},
  \code{\link[PBSsatellite]{print.ncdfData}},
  \code{\link[PBSsatellite]{read.ncdfData}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}% use one of  RShowDoc("KEYWORDS")
\keyword{connection}% __ONLY ONE__ keyword per line

