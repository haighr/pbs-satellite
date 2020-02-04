# Taking cue from Roger Bivand's maptools:
.PBSsatEnv <- new.env(FALSE, parent=globalenv())  # be sure to exportPattern("^\\.PBS") in NAMESPACE

.onAttach <- function(lib, pkg)
{
	pkg_info = utils::sessionInfo( package="PBSsatellite" )$otherPkgs$PBSsatellite
	if( is.character( pkg_info$Packaged ) )
		pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]
	else
		pkg_date  <- date()

	userguide_path <- system.file( "doc/PBSsatellite-UG.pdf", package = "PBSsatellite" )
	year <- substring(date(),nchar(date())-3,nchar(date()))

	packageStartupMessage("
-----------------------------------------------------------
PBS Satellite ", pkg_info$Version, " -- Copyright (C) 2015-",year," Fisheries and Oceans Canada | MacEwan University

A complete user guide 'PBSsatellite-UG.pdf' is located at 
", userguide_path, "

Packaged on ", pkg_date, "
Pacific Biological Station, Nanaimo
-----------------------------------------------------------

")
}
.onUnload <- function(libpath) {
	rm(.PBSsatEnv)
}

# No Visible Bindings
# ===================
if(getRversion() >= "2.15.1") utils::globalVariables(names=c(
	"nc"
	), package="PBSsatellite")

