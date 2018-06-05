% Converts txt IceArea files to netcdf files
%
% Data source: http://neptune.gsfc.nasa.gov/csb/index.php?section=59
%

workdir = '/lustre/atlas1/cli115/proj-shared/milena/observations/SeaIce/IceArea_timeseries';
files =  {'iceAreaNH','iceAreaSH'};

cwd = pwd;
eval(sprintf('cd %s;',workdir));

for i=1:length(files),
  varname = char(files(i));
  infile = sprintf('%s_year.txt',varname);
  outfile1 = sprintf('%s_year.nc', varname);
  outfile2 = sprintf('%s_climo.nc', varname);
  varname = [varname '_year'];

  % the following loads t (yyyy.yearfraction),icearea (km^2)
  eval(sprintf('load %s;',infile));

  % First create netcdf file with original data:
  %
  % create netcdf file
  ncid = netcdf.create(outfile1,'clobber');
  % define dimension(s)
  t_dimid = netcdf.defDim(ncid,'Time',netcdf.getConstant('NC_UNLIMITED'));
  strLen_dimid = netcdf.defDim(ncid,'StrLen',64);
  % define variables and attributes
  t_varid = netcdf.defVar(ncid,'xtime','NC_CHAR',[strLen_dimid,t_dimid]);
  %t_varid = netcdf.defVar(ncid,'time','NC_DOUBLE',t_dimid);
  icearea_varid = netcdf.defVar(ncid,'IceArea','NC_DOUBLE',t_dimid);
  netcdf.putAtt(ncid,t_varid,'long_name','calendar date');
  %netcdf.putAtt(ncid,t_varid,'long_name','days since 0001-01-01');
  netcdf.putAtt(ncid,t_varid,'format','YYYY-MM-DD_HH:MM:SS character string');
  netcdf.putAtt(ncid,icearea_varid,'long_name','SSM/I derived ice area');
  netcdf.putAtt(ncid,icearea_varid,'units','km^2');

  % leave define mode and enter data mode to write data
  netcdf.endDef(ncid);

  % add variables
  eval(sprintf('t = %s(:,1);',varname));
  eval(sprintf('var = %s(:,2);',varname));
  % go from year.fraction_of_year to date (string format)
  tdate = yf2num(t);
  date = datestr(tdate,'yyyy-mm-dd_HH:MM:SS');
  ntot = length(t);
  %netcdf.putVar(ncid,t_varid,0,ntot,t);
  for n=1:ntot,
    netcdf.putVar(ncid,t_varid,[0,n-1],[size(date,2),1],date(n,:));
  end
  netcdf.putVar(ncid,icearea_varid,0,ntot,var);

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
  %t_varid = netcdf.defVar(ncid,'time','NC_DOUBLE',t_dimid);
  icearea_varid = netcdf.defVar(ncid,'IceArea','NC_DOUBLE',t_dimid);
  netcdf.putAtt(ncid,t_varid,'long_name','climatological date');
  %netcdf.putAtt(ncid,t_varid,'long_name','days since 0001-01-01');
  netcdf.putAtt(ncid,t_varid,'format','YYYY-MM-DD_HH:MM:SS character string');
  netcdf.putAtt(ncid,icearea_varid,'long_name','SSM/I derived ice area');
  netcdf.putAtt(ncid,icearea_varid,'units','km^2');

  % leave define mode and enter data mode to write data
  netcdf.endDef(ncid);

  % add variables
  tdatevec = datevec(tdate);
  for im=1:12,
    indmonth = find(tdatevec(:,2)==im);
    var_climo(im) = nanmean(var(indmonth));
  end
  tdate = datenum(1,1:12,15);
  date = datestr(tdate,'yyyy-mm-dd_HH:MM:SS');
  %netcdf.putVar(ncid,t_varid,0,ntot,t);
  for n=1:12,
    netcdf.putVar(ncid,t_varid,[0,n-1],[size(date,2),1],date(n,:));
  end
  netcdf.putVar(ncid,icearea_varid,0,12,var_climo);

  % close netcdf file
  netcdf.close(ncid);
end % loop on data files

eval(sprintf('cd %s;',cwd));
