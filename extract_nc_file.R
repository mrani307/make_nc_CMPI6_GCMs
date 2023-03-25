#Aniruddha Saha
#aniruddha_s@hy.iitr.ac.in
#CODE RUN at line 147

#library ---------------
library(raster)
library(rgdal)
library(ncdf4)
library(rgdal)
  

#functions----------------------

prepare_ncFile<-function(prec,y,var,save.dest) {
  
  #converts a raster stack into a nc file
  #INPUTS : prec - raster stack, 
  #            y - year for which the precipitation/temperature values are stored in the raster stack
  #    save.dest - destination where the .nc file it to be stored
  
  #selecting the variable name and units for the .nc file
  {
    if(var=="Precip")
    {
      varname<-"RAINFALL"
      units<-"mm"
    }
    if(var=="Tmin")
    {
      varname<-"MINIMUM TEMPERATURE"
      units<-"deg C"
    }
    if(var=="Tmax")
    {
      varname<-"MAXIMUM TEMPERATURE"
      units<-"deg C"
    }
  }
  
  # Change NA values to -999
  prec[is.na(prec)] <- (-999)
  
  # output filename
  filename <- paste0(save.dest,"/",var,"_",y,".nc")
  
  # Longitude and Latitude data
  xvals <- unique(values(init(prec, "x")))
  yvals <- unique(values(init(prec, "y")))
  nx <- length(xvals)
  ny <- length(yvals)
  lon <- ncdim_def("LONGITUDE", "degrees_east", xvals)
  lat <- ncdim_def("LATITUDE", "degrees_north", yvals)
  
  # Missing value to use
  mv <- (-999)
  
  # Time component
  time.val.start<-as.Date(paste0(y,"-01-01"))-as.Date("1900-12-31")
  time.val.end<-as.Date(paste0(y+1,"-01-01"))-as.Date("1900-12-31")-1
  time <- ncdim_def(name = "TIME", units = "days since 1900-12-31", 
                    vals = time.val.start:time.val.end, unlim = TRUE)
  
  # Define the variables
  var_prec <- ncvar_def(name = varname,units = units,
                        dim = list(lon, lat, time),missval = mv,compression = NA)
  
  # Add the variables to the file
  ncout <- nc_create(filename, list(var_prec), force_v4 = TRUE)
  
  # add some global attributes
  ncatt_put(ncout, 0, "Title", "Bias corrected climate projections from CMIP6 models for Indian sub-continental river basins")
  ncatt_put(ncout, 0, "Source", "https://zenodo.org/record/3874046#.ZBiWxnZBxEY")
  ncatt_put(ncout, 0, "References", "https://doi.org/10.1038/s41597-020-00681-1")
  
  # Placing the precip/tmax/tmin values in the file need to loop through the layers to get them to match to correct time index
  for (i in 1:nlayers(prec)) { 
    ncvar_put(nc = ncout, varid = var_prec, vals = values(prec[[i]]), 
              start = c(1, 1, i), count = c(-1, -1, 1))
    
  }
  
  # Close the netcdf file when finished adding variables
  nc_close(ncout)
}


makeRaster<-function(r){
  
  #convert XYZ data into raster
  #INPUT : XYZ data
  #OUTPUT : raster with grographic projection
  
    r<-as.data.frame(r)
    sp::coordinates(r)<- ~ long+lat
    sp::gridded(r)<- TRUE
    raster::crs(r)<-"+proj=longlat +datum=WGS84"
    r<-raster::raster(r)
    return(r)
}

convert_to_nc<-function(file.dest,start.year,end.year,var, save.dest, study.region){
  
  #'1.Reads the .txt precip/tmin/tmax file,
  #'2.Extracts year wise data for the whole region, 
  #'3.Converts the year wise data into a raster
  #'4.Masks the raster to your study location  (if required)
  #'5.Converts the raster to .nc file and saves to the given location
  #'6.Prepares .nc files for the mentioned years (if required)
  
  #INPUT  : file.dest    - path of the PrecipData/TMinData/TmaxData file
  #         start.year   - starting year (starts from 2015)
  #         end.year     - ending year   (ends at 2100)
  #         var          - name of variable ("Precip"/"Tmin"/"Tmax")
  #         save.dest    - path of the directory where the yearly .nc files are to be saved
  #         study.region - (OPTIONAL) shapefile for which the .nc files shall be clipped
  #OUTPUT : yearly .nc files shall be saved in save.dest
  
  #read the file
    d<-read.table(file.dest, header = FALSE, sep = " ", dec = ".")
                
  #extracting coordinate data
    coord<-d[1:2,4:ncol(d)]
    rownames(coord)<-c("long","lat")
  
  #extracting date info
    date.df<-d[(3:nrow(d)),1:3]
    colnames(date.df)<-c("year","month","date")
    date<-as.Date(paste0(date.df$year,"-",date.df$month,"-",date.df$date),format="%Y-%m-%d")
    years<-unique(date.df$year)
  
  #removing the date and coordinates from the dataframe
    d<-d[-(1:2),-(1:3)]
    rownames(d)<-c(1:nrow(d))
    
  #checking if the study region data exists
    flag<-1
    if(is.null(study.region)==F)
      flag<-0
  
  #extracting data for each date and storing it year wise
    prepare_raster_and_nc_yearWise<-function(y)
    {
      pos<-which(date.df$year==y)
      r.list<-list()
      for(i in 1:length(pos))
      {
        data<-t(rbind(coord,d[pos[i],])) #prepering xyz data for a particular date for the whole dataset
        r<-makeRaster(data) #converting the xyz data into a raster
         if(flag==1)#checking if a shapefile is provided as argument...
          r<-mask(crop(r,study.region)) #... then masking the raster to the shapefile
        r.list[[i]]<-r
      }
      names(r.list)<-date[pos]
      r.stack<-stack(r.list)   
      prepare_ncFile(prec = r.stack,y= y, var= var,save.dest=save.dest)
    }
    lapply(c(start.year:end.year), prepare_raster_and_nc_yearWise)
}

#RUN-------------------------------

  #reading a shapefile (only if I want the netcdf files to be cropped to my study region)
    shp<-readOGR(dsn="MahiBasin/Datasets/Catchment_Boundary",
                 layer="Mahi_WtBoundary")
  
  #converts the Tmin data for the year 2030 to 2040 and crop it to the study region shapefile

    #INPUT  : file.dest    - path of the "Tmin" file
    #         start.year   - starting year (starts from 2015)
    #         end.year     - ending year   (ends at 2100)
    #         var          - name of variable ("Precip"/"Tmin"/"Tmax")
    #         save.dest    - path of the directory where the yearly .nc files are to be saved
    #         study.region - (OPTIONAL) shapefile for which the .nc files shall be clipped
    #OUTPUT : yearly .nc files shall be saved in save.dest
    convert_to_nc(file.dest = "C:/Users/merul/OneDrive - iitr.ac.in/PhD_Research/GCM/Mahi/ACCESS-CM2/ssp126/TMinData", #path of PrecipData/TMinData/TMaxData file
                  start.year = 2018,
                  end.year = 2020, 
                  var = "Tmin", #either of "Precip","Tmax" or "Tmin"
                  save.dest = "C:/Users/merul/OneDrive - iitr.ac.in/PhD_Research/GCM/Mahi/ACCESS-CM2/ssp126", #path of folder where the yearly .nc files are to be stored
                  #study.region = shp
                  )

