run_field_drainage_analysis <- function() {
  
  # --- 1. DASHBOARD ---
  # Where your field boundary files (KML, SHP, GeoJSON) are stored
  input_folder      <- "input_fields"
  
  # Where the final maps, CSVs, and audit reports will be saved
  output_folder     <- "Field_Zoning_Results"
  
  # Resolution level (13 is ~15-30m pixels). Higher (14) = more detail but more noise.
  # Lower (12) = smoother terrain but loses small draws.
  dem_zoom          <- 13
  
  # Gaussian smoothing intensity. Higher numbers (3.0+) "melt" digital artifacts
  # and neon-tube effects. Lower (1.0) keeps sharp edges but risks noise.
  smoothing_sigma   <- 3.0    
  
  # Post-processing "clumpiness." Higher (3+) creates larger, simpler blocks
  # for machinery. Lower (1) allows for smaller, more complex zone pockets.
  cleanup_level     <- 2      
  
  # The number of management zones to create (Dry, Med, Wet). 
  # Use the Elbow Plot tool to see if you need more than 3.
  num_zones         <- 3
  
  # Vertical distance between labeled contour lines on the map (in meters).
  # 1.0 is standard; use 0.25 or 0.5 for very flat fields.
  contour_interval  <- 1.0     
  
  # If TRUE, deletes the large intermediate TIFF files to save hard drive space.
  # Set to FALSE if you want to keep the raw TWI and Slope rasters for GIS.
  clean_temp        <- TRUE
  
  # --- 2. ENVIRONMENT SETUP ---
  pacman::p_load(terra, sf, elevatr, whitebox, tidyverse, tidyterra, ggspatial, patchwork)
  wbt_init()
  out_base_abs <- file.path(getwd(), output_folder)
  
  all_files <- list.files(path = input_folder, pattern = "\\.(kml|gpkg|shp|geojson)$", full.names = TRUE)
  
  for (i in seq_along(all_files)) {
    file_path  <- all_files[i]; field_name <- tools::file_path_sans_ext(basename(file_path))
    
    tryCatch({
      message(paste0("\n--- Zoning Field: ", field_name, " ---"))
      field_folder <- file.path(out_base_abs, paste0(field_name, "_zones"))
      if(!dir.exists(field_folder)) dir.create(field_folder, recursive = TRUE)
      
      # Step A: Data Prep
      field_sf <- st_read(file_path, quiet=T) %>% st_zm() %>% st_make_valid()
      field_utm <- st_transform(field_sf, as.integer(paste0("326", floor((st_coordinates(st_centroid(st_union(field_sf)))[1] + 180) / 6) + 1)))
      
      # Step B: TWI Generation
      p_dem <- file.path(field_folder, "dem.tif"); t_raw <- file.path(field_folder, "twi_raw.tif")
      dem <- rast(get_elev_raster(field_utm, z=dem_zoom, clip="bbox", progress=F))
      writeRaster(dem, p_dem, overwrite=T)
      
      # Whitebox Pipeline
      s_tif <- file.path(field_folder, "s.tif"); f_tif <- file.path(field_folder, "f.tif")
      sl_tif <- file.path(field_folder, "sl.tif"); ac_tif <- file.path(field_folder, "ac.tif")
      
      wbt_gaussian_filter(p_dem, s_tif, sigma=smoothing_sigma)
      wbt_fill_depressions(s_tif, f_tif)
      wbt_slope(f_tif, sl_tif)
      wbt_d8_flow_accumulation(f_tif, ac_tif)
      wbt_wetness_index(ac_tif, sl_tif, t_raw)
      
      # Explicitly load and mask
      twi <- rast(t_raw) %>% crop(vect(field_utm)) %>% mask(vect(field_utm))
      elev <- rast(f_tif) %>% crop(vect(field_utm)) %>% mask(vect(field_utm))
      names(twi) <- "TWI"
      
      # --- Technique 1: K-Means ---
      twi_df <- as.data.frame(twi, na.rm=TRUE)
      if(nrow(twi_df) == 0) stop("No valid TWI data found within field boundary.")
      
      set.seed(42)
      km <- kmeans(scale(twi_df$TWI), centers=num_zones)
      twi_df$KM_Zone <- as.factor(match(km$cluster, order(km$centers)))
      
      z_km <- twi; values(z_km)[!is.na(values(z_km))] <- as.numeric(twi_df$KM_Zone)
      z_km_clean <- focal(z_km, w=(cleanup_level*2)+1, fun="modal", na.rm=T)
      
      # --- Technique 2: Quantiles ---
      cuts <- quantile(twi_df$TWI, probs = seq(0, 1, length.out = num_zones + 1), na.rm = TRUE)
      twi_df$Quant_Zone <- cut(twi_df$TWI, breaks = cuts, labels = 1:num_zones, include.lowest = TRUE)
      
      z_quant <- twi; values(z_quant)[!is.na(values(z_quant))] <- as.numeric(twi_df$Quant_Zone)
      z_quant_clean <- focal(z_quant, w=(cleanup_level*2)+1, fun="modal", na.rm=T)
      
      # --- Step C: PLOTTING ---
      map_base <- list(
        theme_minimal(),
        annotation_scale(location = "bl"),
        annotation_north_arrow(location = "tr", style = north_arrow_minimal()),
        geom_spatraster_contour(data = elev, binwidth = contour_interval, color = "black", linewidth = 0.2),
        geom_spatraster_contour_text(data = elev, binwidth = contour_interval, size = 2.5, color = "black", check_overlap = TRUE),
        scale_fill_manual(values = c("#d7191c", "#ffffbf", "#2c7bb6"), 
                          name="Zone", labels=c("Dry", "Med", "Wet"), na.translate=F)
      )
      
      p1 <- ggplot() + geom_spatraster(data = as.factor(z_km_clean)) + map_base + labs(title="K-Means (Natural)")
      p2 <- ggplot() + geom_spatraster(data = as.factor(z_quant_clean)) + map_base + labs(title="Quantiles (Equal Area)")
      
      # Audit Plot
      p3 <- ggplot(twi_df) +
        geom_jitter(aes(x = TWI, y = 1, color = KM_Zone), alpha = 0.2, size = 0.4) +
        geom_jitter(aes(x = TWI, y = 2, color = Quant_Zone), alpha = 0.2, size = 0.4) +
        scale_color_manual(values = c("#d7191c", "#ffffbf", "#2c7bb6"), guide = "none") +
        scale_y_continuous(breaks = c(1, 2), labels = c("K-Means", "Quantile")) +
        theme_minimal() + labs(title = "TWI Distribution Audit", x = "TWI Value", y = "")
      
      final_report <- (p1 + p2) / p3 + plot_layout(heights = c(4, 1))
      ggsave(file.path(field_folder, "Zoning_Final.png"), final_report, width=14, height=10, bg="white")
      
      if(clean_temp) unlink(list.files(field_folder, pattern="\\.tif$", full.names=T))
      message("Success!")
      
    }, error = function(e) { message("Error: ", e$message) })
  }
}

run_field_drainage_analysis()

#Slope analysis
run_slope_analysis <- function() {
  
  # --- 1. DASHBOARD ---
  input_folder      <- "input_fields"
  output_folder     <- "Field_Slope_Results"
  dem_zoom          <- 13
  smoothing_sigma   <- 3.0    # Essential to avoid "stepped" slope artifacts
  contour_interval  <- 1.0    
  clean_temp        <- TRUE
  
  # --- 2. ENVIRONMENT SETUP ---
  pacman::p_load(terra, sf, elevatr, whitebox, tidyverse, tidyterra, ggspatial)
  wbt_init()
  out_base_abs <- file.path(getwd(), output_folder)
  
  all_files <- list.files(path = input_folder, pattern = "\\.(kml|gpkg|shp|geojson)$", full.names = TRUE)
  
  for (i in seq_along(all_files)) {
    file_path  <- all_files[i]; field_name <- tools::file_path_sans_ext(basename(file_path))
    
    tryCatch({
      message(paste0("\n--- Mapping Slope: ", field_name, " ---"))
      field_folder <- file.path(out_base_abs, paste0(field_name, "_slope"))
      if(!dir.exists(field_folder)) dir.create(field_folder, recursive = TRUE)
      
      # Step A: Data Prep
      field_sf <- st_read(file_path, quiet=T) %>% st_zm() %>% st_make_valid()
      field_utm <- st_transform(field_sf, as.integer(paste0("326", floor((st_coordinates(st_centroid(st_union(field_sf)))[1] + 180) / 6) + 1)))
      
      # Step B: DEM & Slope Calculation
      p_dem <- file.path(field_folder, "dem.tif"); s_tif <- file.path(field_folder, "smoothed_dem.tif")
      dem <- rast(get_elev_raster(field_utm, z=dem_zoom, clip="bbox", progress=F))
      writeRaster(dem, p_dem, overwrite=T)
      
      # Use Whitebox to smooth the DEM first (raw DEMs produce very "noisy" slope)
      wbt_gaussian_filter(p_dem, s_tif, sigma=smoothing_sigma)
      
      # Calculate Slope (in degrees)
      elev <- rast(s_tif) %>% crop(vect(field_utm)) %>% mask(vect(field_utm))
      slope_rast <- terrain(elev, v="slope", unit="degrees")
      names(slope_rast) <- "Slope"
      
      # Step C: Plotting
      slope_plot <- ggplot() +
        geom_spatraster(data = slope_rast) +
        # "Inferno" or "Magma" scales are great for slope—dark is flat, bright is steep
        scale_fill_viridis_c(option = "inferno", name = "Slope (°)", na.value = "transparent") +
        geom_spatraster_contour(data = elev, binwidth = contour_interval, color = "white", alpha = 0.3, linewidth = 0.2) +
        geom_spatraster_contour_text(data = elev, binwidth = contour_interval, size = 2, color = "white", alpha = 0.6) +
        annotation_scale(location = "bl") +
        annotation_north_arrow(location = "tr", style = north_arrow_minimal()) +
        theme_minimal() +
        labs(title = paste("Slope Analysis:", field_name),
             subtitle = paste("Elevation Contours:", contour_interval, "m | Data: AWS Terrain Tiles"),
             caption = paste("Projection:", st_crs(field_utm)$input))
      
      ggsave(file.path(field_folder, "Slope_Map.png"), slope_plot, width=10, height=10, bg="white")
      
      if(clean_temp) unlink(list.files(field_folder, pattern="\\.tif$", full.names=T))
      message("Success!")
      
    }, error = function(e) { message("Error: ", e$message) })
  }
}

# Run the slope analysis
run_slope_analysis()
