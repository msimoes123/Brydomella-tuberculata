setwd("C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp")

#Organizing occ--------
library(rworldmap)
pts_A <- read.delim("Extant_.txt")
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}
pts_a <- round_df(pts_A, 3) 
dt <- na.omit(unique( pts_a[ , 1:2 ] ))
dt$species <- 'brydo'
dt <- cbind(dt$species, dt$longitude, dt$latitude)
colnames(dt) <- c('species', 'long', 'lat')
write.csv(dt, 'extant.csv', row.names = F )

pts_I <- read.delim("Extinct_.txt")
pts_i <- round_df(pts_I, 3) 
dt <- na.omit(unique( pts_i[ , 1:2 ] ))
dt$species <- 'brydo'
dt <- cbind(dt$species, dt$longitude, dt$latitude)
colnames(dt) <- c('species', 'long', 'lat')
write.csv(dt, 'extinct.csv', row.names = F )


a <- read.csv('extant.csv') 
i<- read.csv('extinct.csv')
b <- rbind(a,i) #all records
write.csv(b, 'total.csv', row.names = F )

#Exploring options for M size----------
library(rgeos)
library(rgdal)
library(dplyr)
library(ggplot2)
t <- SpatialPoints(b[,c("long","lat")])
x <-gDistance(t, t, byid = TRUE)
write.csv(x, 'DistanceOcc.csv', row.names = F)
x <- read.csv('DistanceOcc.csv') #matrix of distances
max(x)
hist( x[x > 0 & x < 150], main = 'Distances between Occurences of Brydomella', breaks=30) #5 and 10km and more frequent than any other


#Visualization---------
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(0, 120), ylim = c(40, 71), asp = 1)
points(b$long, b$lat, col = "blue", cex = .9)
points(b$longitude, b$latitude, col = "red", cex = .9)
r <- raster("C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\Bios_all\\bios_brydo\\bio1.asc")
d <- read.csv('total.csv')
t <- extract(r, d[,2:3])
x <- cbind(d, t)
c<-na.omit(x)
write.csv(c[,1:3], 'total.csv', row.names = F)
b<-read.csv('total.csv')

#Crop variables----
#mask variables Worldclim
require(raster)
path <- ".\\Bios_all\\"
varaibles_list <-list.files(path = path, pattern = ".bil", full.names = TRUE) #variables 
variables <- stack(varaibles_list) #create a stack
require(rgdal)
shape <- readOGR(dsn = "C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp",layer = "brydo_cont") #SHAPEFILE YOU CREATED FOR PROJECTION AREA
var_mask <- mask(crop(variables, shape), shape)

## names for layers
rnames <- paste0("C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\Bios_all\\bios_brydo\\", names(variables), ".asc") # users select the format

## saving layers in new folder
sav <- lapply(1:nlayers(var_mask), function(x) {
  writeRaster(var_mask[[x]], filename = rnames[x], format = "ascii") # change format accordingly
})

# PCA------
library(kuenm)
var_folder <- "C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\Bios_all\\bios_brydo" # name of folder with variables to be combined in distinct sets
out_folder <- "C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\pc_all" # name of folder that will contain the sets
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
npcs <- 8 # number of pcs you want as rasters, if not defined all pcs are returned as rasters

# PCA of variables for models
kuenm_rpca(variables  = var_folder, in.format = in_format, var.scale = TRUE, write.result = TRUE, 
           out.format = "ascii", out.dir = out_folder, n.pcs = npcs)

#M mask-------
#Based on the distance analysis we will create 2 buffers:5 and 10km 
#head(b)
#variables <- raster::stack(list.files("C:/Users/mari1/Desktop/Lara_ExtinctSp/pc_temp/Initial/", 
 #                                     pattern = "pc", full.names = TRUE))
#b_5 <- buffer_area(b, longitude = "long", latitude = "lat", buffer_distance = 50)
#b_10 <- buffer_area(b, longitude = "long", latitude = "lat", buffer_distance = 10)
#b_100 <- buffer_area(b, longitude = "long", latitude = "lat", buffer_distance = 100)
#plot(b_100)
#var_mask <- mask(crop(variables, b_100), b_100)
#rnames <- paste0("C:/Users/mari1/Desktop/Lara_ExtinctSp/Ms_variables/100_all/", names(variables), ".asc") # 
## saving layers in new folder
sav <- lapply(1:nlayers(var_mask), function(x) {
  writeRaster(var_mask[[x]], filename = rnames[x], format = "ascii", overwrite=T) # change format accordingly
})


#Ellispsoids-------------

if(!require(devtools)){
  install.packages("devtools")
}
if(!require(ellipsenm)){
  devtools::install_github("marlonecobos/ellipsenm")
}
library(ellipsenm)
setwd("C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp")
# raster layers of environmental data (this ones are masked to the accessible area)
# users must prepare their layers accordingly if using other data
vars <- raster::stack(list.files("C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\pc_all\\Initial", pattern = "pc", full.names = TRUE))
crs(vars) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 '


# preparing training and testing data
data_split <- split_data(b, method = "random", longitude = "long", 
                         latitude = "lat", train_proportion = 0.75, 
                         save = TRUE, name = "occ")

# sets of variables (example)
setwd('C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\pc_all\\Initial')
sets <- list(set_1 = c('pc_1', 'pc_2', 'pc_3', 
                       'pc_4', 'pc_5', 'pc_6')) # change as needed
variable_sets <- prepare_sets(vars, sets)

# methods to create ellipsoids
methods <- c("covmat", "mve1")

# model calibration process
calib <- ellipsoid_calibration(data_split, species = "species", longitude = "long", 
                               latitude = "lat", variables= variable_sets, methods = methods, level = 99, 
                               selection_criteria = "S_OR_P",
                               error = 5, iterations = 500, percentage = 50,
                               output_directory = "calibration_results2")

#replicated was used default = "bootstrap"; prediction mahala and suitability; return_numeric - TRUE to get Mahal v
ell_model <- ellipsoid_model(data = b, species = "species",
                             longitude = "long", latitude = "lat",
                             replicate_type = "bootstrap", bootstrap_percentage = 75,
                             raster_layers = vars,  
                             method = "covmat",  level = 99,
                             replicates = 10, prediction = "both", 
                             return_numeric = TRUE, format = "GTiff", color_palette = viridis::magma, 
                             overwrite = TRUE, output_directory = "ellipsenm_model_pcall5")

#We will use the mean, so no need to calculate the median- takes too long
#rast <- raster::stack(list.files("C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\ellipsenm_model_pcall5\\mahalanobis", pattern = ".tif", full.names = TRUE))
#crs(vars) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 '
#r_median <- calc(rast, median)
#writeRaster(r_median, 'C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\ellipsenm_model_pcall5\\mahalanobis\\median', format="GTiff")
#extract values to all points
a <- read.csv('extant.csv')
r <- raster("C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\Bios_all\\bios_brydo\\bio1.asc")
g<- extract(r, a[,2:3])
x<-cbind(a,g)
a <-na.omit(x)
i$status <- 'extinct'
i <- read.csv('extinct.csv')
a <- a[,1:3]
a$status <- 'extant'
head(i)
t <- rbind(a,i)
write.csv(h, 'MahalanobisDistance.csv', row.names = F)
head(h)
#Pcs6
require(raster)
r <- raster("C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\ellipsenm_model_pcall6\\mean_mahalanobis_calibration_brydo.tif.tif")
pc6 <- extract(r, t[,2:3])
r <- raster("C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\ellipsenm_model_pcall5\\mean_mahalanobis_calibration_brydo.tif.tif")
pc5 <- extract(r, t[,2:3])
r <- raster("C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\ellipsenm_model_pctemp\\mean_mahalanobis_calibration_brydo.tif.tif")
pcT <- extract(r, t[,2:3])
r <- raster("C:\\Users\\mari1\\Desktop\\Lara_ExtinctSp\\ellipsenm_model_pcprec\\mean_mahalanobis_calibration_brydo.tif.tif")
pcP <- extract(r, t[,2:3])
h<- cbind(t, pc6,pc5,pcT,pcP) #all distances in a matrix
require(ggplot2)
i<- h[ which(h$status=='extinct'), ]
a<- h[ which(h$status=='extant'), ]
x<- i$pc6
ggplot(a, aes(x=pc6)) +
  geom_histogram(fill="white", color="black")+
  labs(title="Mahalanobis Distance",x="Distance (D2) ", y = "Number of Occ records")+
  abline(v=x, col="red", lty=2, lwd=3)+
  theme_classic()


head(h)
p<-ggplot(h, aes(x=pcT, color=status, fill=status)) +
  geom_histogram(alpha=0.8, position="identity", binwidth = 0.5)+
  labs(title="Principal components of Temperature Variables",x="Mahalanobis Distance (D2) ", y = "Number of Occ records")+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0,0))+
  coord_cartesian(ylim = c(0, 40)) +   theme_minimal() 
p

  