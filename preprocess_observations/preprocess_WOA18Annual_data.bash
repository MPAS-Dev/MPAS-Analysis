#!/bin/bash

outfile='woa18_decav_04_TS_ann.nc'

infileTemp=./annual/woa18_decav_t00_04.nc
infileSalt=./annual/woa18_decav_s00_04.nc
ncwa -O -a time ${infileTemp} tmp1.nc
ncwa -O -a time ${infileSalt} tmp2.nc

ncks -O -v t_an tmp1.nc ${outfile}
ncks -A -v s_an tmp2.nc ${outfile}

# Clean up
rm -f tmp*nc
