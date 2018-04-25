% Converts txt PIOMAS data to netcdf files
%
% Data source: http://psc.apl.uw.edu/research/projects/arctic-sea-ice-volume-anomaly/
%

clear all;
close all;

workdir = '/lustre/atlas1/cli115/proj-shared/milena/observations/SeaIce/PIOMAS';
infile = 'PIOMASvolume_monthly';
outfile1 = [infile '.nc'];
outfile2 = [infile '_climo.nc'];

cwd = pwd;
eval(sprintf('cd %s;',workdir));

% the following loads year (first column) and ice volume (in 10^3 km^3)
% for each month of that year (subsequent 12 columns):
eval(sprintf('load %s.txt;',infile));

eval(sprintf('years = %s(:,1);',infile));
eval(sprintf('obs = %s(:,2:13);',infile));

obs(find(obs==-1)) = nan;

% Reorganize monthly observations:
ind = 1;
for iy=1:length(years),
  for im=1:12,
    tdate(ind)  = datenum(years(iy),im,15);
    date(ind,:) = datestr(tdate(ind),'yyyy-mm-dd_HH:MM:SS');
    icevol(ind) = obs(iy,im);
    ind = ind+1;
  end
end

% Compute climatological annual cycle:
tdatevec = datevec(tdate);
for im=1:12,
  indmonth = find(tdatevec(:,2)==im);
  icevol_climo(im) = nanmean(icevol(indmonth));
end
tdate_climo = datenum(1,1:12,15);
date_climo  = datestr(tdate_climo,'yyyy-mm-dd_HH:MM:SS');

% First create netcdf file with monthly data:
%
% create netcdf file
ncid = netcdf.create(outfile1,'clobber');
% define dimension(s)
t_dimid = netcdf.defDim(ncid,'Time',netcdf.getConstant('NC_UNLIMITED'));
strLen_dimid = netcdf.defDim(ncid,'StrLen',64);
% define variables and attributes
t_varid = netcdf.defVar(ncid,'xtime','NC_CHAR',[strLen_dimid,t_dimid]);
icevol_varid = netcdf.defVar(ncid,'IceVol','NC_DOUBLE',t_dimid);
netcdf.putAtt(ncid,t_varid,'long_name','calendar date');
netcdf.putAtt(ncid,t_varid,'format','YYYY-MM-DD_HH:MM:SS character string');
netcdf.putAtt(ncid,icevol_varid,'long_name','PIOMAS ice volume');
netcdf.putAtt(ncid,icevol_varid,'units','10^3 km^3');

% leave define mode and enter data mode to write data
netcdf.endDef(ncid);

% add variables
ntot = ind-1;
for n=1:ntot,
  netcdf.putVar(ncid,t_varid,[0,n-1],[size(date,2),1],date(n,:));
end
netcdf.putVar(ncid,icevol_varid,0,ntot,icevol);

% close netcdf file
netcdf.close(ncid);

% Then create netcdf file with climatological annual cycle:
%
% create netcdf file
ncid = netcdf.create(outfile2,'clobber');
% define dimension(s)
t_dimid = netcdf.defDim(ncid,'Time',netcdf.getConstant('NC_UNLIMITED'));
strLen_dimid = netcdf.defDim(ncid,'StrLen',64);
% define variables and attributes
t_varid = netcdf.defVar(ncid,'xtime','NC_CHAR',[strLen_dimid,t_dimid]);
icevol_varid = netcdf.defVar(ncid,'IceVol','NC_DOUBLE',t_dimid);
netcdf.putAtt(ncid,t_varid,'long_name','climatological date');
netcdf.putAtt(ncid,t_varid,'format','YYYY-MM-DD_HH:MM:SS character string');
netcdf.putAtt(ncid,icevol_varid,'long_name','PIOMAS ice volume');
netcdf.putAtt(ncid,icevol_varid,'units','10^3 km^3');

% leave define mode and enter data mode to write data
netcdf.endDef(ncid);

% add variables
for n=1:12,
  netcdf.putVar(ncid,t_varid,[0,n-1],[size(date_climo,2),1],date_climo(n,:));
end
netcdf.putVar(ncid,icevol_varid,0,12,icevol_climo);

% close netcdf file
netcdf.close(ncid);

eval(sprintf('cd %s;',cwd));
