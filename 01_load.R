# Copyright 2018 Province of British Columbia
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.

source("header.R")

# Load Provincial Boundary, Land Form, BTM, GBPU, CE data,
# Livestock Density, Human Density, GBPU

#Rasterize the Province for subsequent masking
ProvRast<-raster(nrows=15744, ncols=17216, xmn=159587.5, xmx=1881187.5,
                 ymn=173787.5, ymx=1748187.5,
                 crs="+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0",
                 res = c(100,100), vals = 0)

BCr_file <- file.path(dataOutDir,"BCr.tif")
if (!file.exists(BCr_file)) {
  BCr <- fasterize(bcmaps::bc_bound_hres(class='sf'),ProvRast)
  writeRaster(BCr, filename=BCr_file, format="GTiff", overwrite=TRUE)
} else {
  BCr <- raster(BCr_file)
}

#BTM - for settlement, rural, agriculture, range, mining
VRI_file <- file.path("tmp/VRI_file")
if (!file.exists(VRI_file)) {
  # Link to VRI file download from BCDC:
  # devtools: (idea) nstall_github("bcgov/bcgovr", ref = "devel")

  #Dowload file manually and put *.zip in this script and place file in the data directory
 VRIZip <- 'VEG_COMP_LYR_R1_POLY.gdb.zip'
  unzip(file.path(DataDir, VRIZip), exdir = file.path(DataDir, "VRI"))

  # List feature classes in the geodatabase
  VRI_gdb <- list.files(file.path(DataDir, "FullVRI"), pattern = ".gdb", full.names = TRUE)[1]
  fc_list <- st_layers(VRI_gdb)

  #I have the same memory issue - need to allocate more memory to work
  VRI <- read_sf(VRI_gdb, layer = "VEG_COMP_LYR_R1_POLY")

  saveRDS(VRI, file = VRI_file)


  #
  #BCLCS_LEVEL_1
  #BCLCS_LEVEL_2
  #BCLCS_LEVEL_3
  #BCLCS_LEVEL_4
  #BCLCS_LEVEL_5
  #
  BCLCS_L1 <- VRI[VRI$BCLCS_LEVEL_1] %>%
    fasterize(ProvRast, background=NA)

#read in shape
  TestVRI <- read_sf(file.path(DataDir,'VRI/VRI.shp'))
#make a LUT and join codes for each BCLCS_LE_1 type
  BCLCS_LE_1_LUT<-data.frame(BCLCS_LE_1=c('L','N','T','W'),
                             BCLCS_LE_1Code=c(1,2,3,4))
  TestVRI_1<-left_join(TestVRI,BCLCS_LE_1_LUT,by='BCLCS_LE_1')
#fasterize the BCLCS_LE_1 code
  BCLCS_LE_1r <- fasterize(TestVRI_1, ProvRast, field = 'BCLCS_LE_1Code', background=0)
