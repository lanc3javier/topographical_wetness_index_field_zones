# topographical_wetness_index_field_zones
Some code that uses .kml field boundaries, grabs DEM data, and splits a field into zones using k-means analysis and a median-based approach. Used to compare the two methods for zone assignment of pixels in the raster.

Ag-Terrain-Zoner: Topographic Wetness & Slope Analysis
This tool uses R to turn field boundaries into drainage and management maps. It calculates where water flows and pools based on the shape of the land.

What it does
The script takes a field boundary and automatically downloads elevation data. It then calculates the Topographic Wetness Index (TWI). This index identifies ridges (dry) and draws (wet).

The tool provides two different ways to divide the field:

Natural Groups (K-Means): This groups the soil by its actual characteristics. If most of the field is similar, the middle zone will be the largest. This is best for deciding where to take soil samples.

Equal Areas (Quantiles): This forces the field into equal parts (33% each).

How to use it
Put your field files (KML, Shapefile, or GeoJSON) into a folder named input_fields.

Run the script. It will create a results folder for each field.

Check the Zoning_Report_Final.png to see the maps side-by-side with elevation contours.

Settings you can change
Smoothing: Adjusts how "noisy" the map is. Higher smoothing makes for cleaner zones that are easier for a tractor to follow.

Contours: You can change how often elevation labels appear (for example, every 1 meter or every 0.5 meters).

Number of Zones: Usually set to 3 (Dry, Medium, Wet), but you can increase this if the field is very complex.

Required R Packages
You will need to install these libraries to run the code:

terra (for map math)

whitebox (for hydrology calculations)

sf (for handling boundaries)

elevatr (for downloading elevation data)

tidyverse & tidyterra (for making the maps)

Disclaimer
These maps are mathematical estimates based on satellite elevation data. Always verify the results in the field with a shovel or soil probe before installing tile drainage or changing your management practices.
