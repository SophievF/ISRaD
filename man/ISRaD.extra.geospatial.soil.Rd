% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ISRaD.extra.geospatial.soil.R
\name{ISRaD.extra.geospatial.soil}
\alias{ISRaD.extra.geospatial.soil}
\title{ISRaD.extra.geospatial.soil}
\usage{
ISRaD.extra.geospatial.soil(database, geodata_soil_directory)
}
\arguments{
\item{database}{soilcarbon dataset object}

\item{geodata_soil_directory}{directory where geospatial soil data are found}
}
\value{
returns new ISRaD_extra object with extracted bulk density and clay, carbon, and coarse fragment content up to 200 cm soil depth
}
\description{
Extracts modeled soil properties from 250m resolution Soil Grids spatial products
}
\details{
Uses filled geographic coordinates of profiles to extract estimated (observations + machine learning predictions) clay content (kg/kg), organic carbon content (x 5 g/kg), carbon stock (kg/m2), bulk density (kg/m3), and coarse fragments (% volumetric) from rasters at 0, 10, 30, 60, 100, and 200 cm soil depth. For function to work, you must first download the appropriate raster files from SoilGrids (doi: 10.5281/zenodo.2525663, doi.org/10.5281/zenodo.1475457, doi.org/10.5281/zenodo.1475970, doi.org/10.5281/zenodo.2525681, doi.org/10.5281/zenodo.1475453). To convert organic carbon content to %, divide by 2. In added columns, "SG" denotes "Soil Grids". For more information see ISRaD Extra info file at <http://soilradiocarbon.org>
}
