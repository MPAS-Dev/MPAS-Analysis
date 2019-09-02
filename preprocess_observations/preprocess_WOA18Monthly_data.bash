#!/bin/bash

outfile='woa18_decav_04_TS_mon.nc'

infilesTemp=`ls monthly/woa18_decav_t*.nc`
infilesSalt=`ls monthly/woa18_decav_s*.nc`

#months=(1 2 3 4 5 6 7 8 9 10 11 12)
#months=($(seq 1 1 12))

i=1
for infile in ${infilesTemp[@]}; do
  ncwa -O -a time ${infile} tmp${i}.nc
  i=$[${i}+1]
done
# Make sure monthly files are added in the right order:
ncecat -O -u month tmp1.nc tmp2.nc tmp3.nc tmp4.nc tmp5.nc tmp6.nc tmp7.nc \
 tmp8.nc tmp9.nc tmp10.nc tmp11.nc tmp12.nc allmonthsT.nc
ncap2 -O -s 'month=array(1,1,$month)' allmonthsT.nc allmonthsT.nc
ncatted -a standard_name,month,a,c,'climatological month' allmonthsT.nc

i=1
for infile in ${infilesSalt[@]}; do
  ncwa -O -a time ${infile} tmp${i}.nc
  i=$[${i}+1]
done
# Make sure monthly files are added in the right order:
ncecat -O -u month tmp1.nc tmp2.nc tmp3.nc tmp4.nc tmp5.nc tmp6.nc tmp7.nc \
 tmp8.nc tmp9.nc tmp10.nc tmp11.nc tmp12.nc allmonthsS.nc
ncap2 -O -s 'month=array(1,1,$month)' allmonthsS.nc allmonthsS.nc
ncatted -a standard_name,month,a,c,'climatological month' allmonthsS.nc

ncks -O -v t_an allmonthsT.nc ${outfile}
ncks -A -v s_an allmonthsS.nc ${outfile}

# Clean up
rm -f tmp*nc allmonths*nc
