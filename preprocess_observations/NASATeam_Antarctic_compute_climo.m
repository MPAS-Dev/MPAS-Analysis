%
% Combined SSM/I-SSMIS and SMMR sea-ice concentration data for both the 
% Arctic and Antarctic.
% *** NASA Team algorithm ***
%
% Source: http://nsidc.org/data/NSIDC-0051
%
clear all;
close all;

homedir = getenv('HOME');

maindir = [homedir '/ACME/observations/obsdir/SeaIce/SSMI/'...
           'NASATeam_NSIDC0051'];
datadir = [maindir '/south/monthly'];
griddir = [homedir '/ACME/observations/obsdir/SeaIce/ICESat/Antarctic/climo'];
lonlatfile = 'spring_ICESat_gridded_mean_thickness';
months4climo = [12  1  2;  % (boreal) winter months
                 3  4  5;  % (boreal) spring months
                 6  7  8;  % (boreal) summer months
                 9 10 11]; % (boreal) fall   months
seasons = {'djf','mam','jja','son'};

M = 332; % Polar Stereographic grid for Antarctica, # of rows
N = 316; % Polar Stereographic grid for Antarctica, # of columns

[err,cwd] = unix('pwd');
% Read in lon,lat:
eval(sprintf('cd %s;',griddir));
eval(sprintf('load %s.txt',lonlatfile));
eval(sprintf('cd %s;',cwd));
eval(sprintf('lat = %s(:,1);',lonlatfile));
eval(sprintf('lon = %s(:,2);',lonlatfile));
x = reshape(lon,M,N);
y = reshape(lat,M,N);

eval(sprintf('cd %s;',datadir));
for is=1:length(seasons),
  season = char(seasons(is));

  % Get list of files to compute seasonal climatologies from:
  filenames = sprintf('nt_????%02d*',months4climo(is,1));
  for i=2:size(months4climo,2),
    filenames = [filenames sprintf(' nt_????%02d*',months4climo(is,i))];
  end
  eval(sprintf('[err,filelist] = unix(''/bin/ls %s'');',filenames));
  filelist = strsplit(filelist); % convert to cell of strings
  filelist = filelist(2:end-1); % remove redudant first and last space chars

  % Compute seasonal climatology:
  for ifile=1:length(filelist);
    infile = char(filelist(ifile));
    % Read in binary data:
    fid = fopen(infile,'r');
    header = fread(fid,300,'char=>char');
    aice(ifile,:,:) = fread(fid,[N M],'uint8');
    fclose(fid);
    aice(find(aice>250)) = nan;
    aice(find(aice==0))=nan;
  end
  aice = squeeze(nanmean(aice,1));
  aice = aice/250;
  aice(find(isnan(aice)==1)) = -999;
  aice = aice';

  % Write to file:
  outfile = ['SSMI_NASATeam_gridded_concentration_SH_' season];
  fid = fopen(sprintf('%s/%s.txt',maindir,outfile),'w');
  fprintf(fid,'%12.4f %12.4f %12.5f\n',[y(:)'; x(:)'; aice(:)']);
  fclose(fid);

  clear aice
end
eval(sprintf('cd %s;',cwd));
