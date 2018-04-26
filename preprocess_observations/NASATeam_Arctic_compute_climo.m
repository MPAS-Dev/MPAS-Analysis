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
datadir = [maindir '/north/monthly'];
griddir = [homedir '/ACME/observations/obsdir/SeaIce/ICESat/Arctic/' ...
           'NSIDC0393_GLAS_SI_Freeboard_v01/glas_seaice_grids'];
lonfile = [griddir '/PS25km_north_lon.img'];
latfile = [griddir '/PS25km_north_lat.img'];
months4climo = [ 1  2  3;  % (boreal) winter months
                 4  5  6;  % (boreal) spring months
                 7  8  9;  % (boreal) summer months
                10 11 12]; % (boreal) fall   months
seasons = {'jfm','amj','jas','ond'};

M = 304; % Polar Stereographic grid for the Arctic, # of rows
N = 448; % Polar Stereographic grid for the Arctic, # of columns

% Read in lon,lat:
fid = fopen(lonfile,'r','l');
x = fread(fid,[M N],'single');
fclose(fid);
fid = fopen(latfile,'r','l');
y = fread(fid,[M N],'single');
fclose(fid);

[err,cwd] = unix('pwd');
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
    aice(ifile,:,:) = fread(fid,[M N],'uint8');
    fclose(fid);
    aice(find(aice>250)) = nan;
    aice(find(aice==0))=nan;
  end
  aice = squeeze(nanmean(aice,1));
  aice = aice/250;
  aice(find(isnan(aice)==1)) = -999;

  % Write to file:
  outfile = ['SSMI_NASATeam_gridded_concentration_NH_' season];
  fid = fopen(sprintf('%s/%s.txt',maindir,outfile),'w');
  fprintf(fid,'%12.4f %12.4f %12.5f\n',[y(:)'; x(:)'; aice(:)']);
  fclose(fid);

  clear aice
end
eval(sprintf('cd %s;',cwd));
