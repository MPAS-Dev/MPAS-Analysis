### make single file from yearly Altiberg merged files
### downloaded from 
### ftp://ftp.ifremer.fr/ifremer/cersat/projects/altiberg/v2/data/merged/grid/geographic/

# make time record dimension
ncks --mk_rec_dmn time prod_latlon_merged_1991.nc prod_latlon_merged_1991_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_1992.nc prod_latlon_merged_1992_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_1993.nc prod_latlon_merged_1993_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_1994.nc prod_latlon_merged_1994_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_1995.nc prod_latlon_merged_1995_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_1996.nc prod_latlon_merged_1996_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_1997.nc prod_latlon_merged_1997_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_1998.nc prod_latlon_merged_1998_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_1999.nc prod_latlon_merged_1999_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2000.nc prod_latlon_merged_2000_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2001.nc prod_latlon_merged_2001_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2002.nc prod_latlon_merged_2002_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2003.nc prod_latlon_merged_2003_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2004.nc prod_latlon_merged_2004_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2005.nc prod_latlon_merged_2005_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2006.nc prod_latlon_merged_2006_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2007.nc prod_latlon_merged_2007_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2008.nc prod_latlon_merged_2008_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2009.nc prod_latlon_merged_2009_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2010.nc prod_latlon_merged_2010_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2011.nc prod_latlon_merged_2011_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2012.nc prod_latlon_merged_2012_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2013.nc prod_latlon_merged_2013_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2014.nc prod_latlon_merged_2014_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2015.nc prod_latlon_merged_2015_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2016.nc prod_latlon_merged_2016_timerec.nc
ncks --mk_rec_dmn time prod_latlon_merged_2017.nc prod_latlon_merged_2017_timerec.nc

# concatanate files
ncrcat prod_latlon_merged_1991_timerec.nc \
prod_latlon_merged_1992_timerec.nc \
prod_latlon_merged_1993_timerec.nc \
prod_latlon_merged_1994_timerec.nc \
prod_latlon_merged_1995_timerec.nc \
prod_latlon_merged_1996_timerec.nc \
prod_latlon_merged_1997_timerec.nc \
prod_latlon_merged_1998_timerec.nc \
prod_latlon_merged_1999_timerec.nc \
prod_latlon_merged_2000_timerec.nc \
prod_latlon_merged_2001_timerec.nc \
prod_latlon_merged_2002_timerec.nc \
prod_latlon_merged_2003_timerec.nc \
prod_latlon_merged_2004_timerec.nc \
prod_latlon_merged_2005_timerec.nc \
prod_latlon_merged_2006_timerec.nc \
prod_latlon_merged_2007_timerec.nc \
prod_latlon_merged_2008_timerec.nc \
prod_latlon_merged_2009_timerec.nc \
prod_latlon_merged_2010_timerec.nc \
prod_latlon_merged_2011_timerec.nc \
prod_latlon_merged_2012_timerec.nc \
prod_latlon_merged_2013_timerec.nc \
prod_latlon_merged_2014_timerec.nc \
prod_latlon_merged_2015_timerec.nc \
prod_latlon_merged_2016_timerec.nc \
prod_latlon_merged_2017_timerec.nc \
-o Altiberg_1991-2017.nc
