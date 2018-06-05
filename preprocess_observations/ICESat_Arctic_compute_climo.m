%
% ICESat sea-ice thickness data for the Arctic.
%
% Source: http://nsidc.org/data/NSIDC-0393
%
clear all;
close all;

homedir = getenv('HOME');

datadir = [homedir '/ACME/observations/obsdir/SeaIce/ICESat/Arctic/' ...
           'NSIDC0393_GLAS_SI_Freeboard_v01/glas_seaice_grids'];
% Infile names (austral spring (Oct-Nov), max Antarctic ice extent;
%               austral late summer (Feb-Mar), min Antarctic ice extent):
on_files_th = {'laser2a_thickness_mskd.img','laser3a_thickness_mskd.img',...
               'laser3d_thickness_mskd.img','laser3g_thickness_mskd.img',...
               'laser3i_thickness_mskd.img','laser3k_thickness_mskd.img'};
fm_files_th = {'laser1_thickness_mskd.img' ,'laser2b_thickness_mskd.img',...
               'laser3b_thickness_mskd.img','laser3e_thickness_mskd.img',...
               'laser3h_thickness_mskd.img','laser3j_thickness_mskd.img'};
on_files_fb = {'laser2a_freeboard_mskd.img','laser3a_freeboard_mskd.img',...
               'laser3d_freeboard_mskd.img','laser3g_freeboard_mskd.img',...
               'laser3i_freeboard_mskd.img','laser3k_freeboard_mskd.img'};
fm_files_fb = {'laser1_freeboard_mskd.img' ,'laser2b_freeboard_mskd.img',...
               'laser3b_freeboard_mskd.img','laser3e_freeboard_mskd.img',...
               'laser3h_freeboard_mskd.img','laser3j_freeboard_mskd.img'};
outfiles = {'ICESat_gridded_mean_thickness_NH_on',...
            'ICESat_gridded_mean_thickness_NH_fm'};
lonfile = 'PS25km_north_lon.img';
latfile = 'PS25km_north_lat.img';

M = 304; % Polar Stereographic grid for the Arctic, # of rows
N = 448; % Polar Stereographic grid for the Arctic, # of columns

eval(sprintf('cd %s;',datadir));

% Read in lon,lat:
fid = fopen(lonfile,'r','l');
x = fread(fid,[M N],'single');
fclose(fid);
fid = fopen(latfile,'r','l');
y = fread(fid,[M N],'single');
fclose(fid);

% Compute Oct-Nov climatologies:
for ifile=1:length(on_files_th),
  on_file = char(on_files_th(ifile));
  fid = fopen(on_file,'r','l');
  hi(ifile,:,:) = fread(fid,[M N],'single');
  fclose(fid);
end
hi(find(hi<0)) = nan;
hi = squeeze(nanmean(hi,1));
hi(find(isnan(hi)==1)) = -999;
for ifile=1:length(on_files_fb),
  on_file = char(on_files_fb(ifile));
  fid = fopen(on_file,'r','l');
  fb(ifile,:,:) = fread(fid,[M N],'single');
  fclose(fid);
end
fb(find(fb<0)) = nan;
fb = squeeze(nanmean(fb,1));
fb(find(isnan(fb)==1)) = -999;
% Write to file:
outfile = char(outfiles(1));
fid = fopen(sprintf('%s.txt',outfile),'w');
fprintf(fid,'%12.4f %12.4f %12.5f %12.5f\n',[y(:)'; x(:)'; fb(:)'; hi(:)']);
fclose(fid);

clear hi fb

% Compute Feb-Mar climatologies:
for ifile=1:length(fm_files_th),
  fm_file = char(fm_files_th(ifile));
  fid = fopen(fm_file,'r','l');
  hi(ifile,:,:) = fread(fid,[M N],'single');
  fclose(fid);
end
hi(find(hi<0)) = nan;
hi = squeeze(nanmean(hi,1));
hi(find(isnan(hi)==1)) = -999;
for ifile=1:length(fm_files_fb),
  fm_file = char(fm_files_fb(ifile));
  fid = fopen(fm_file,'r','l');
  fb(ifile,:,:) = fread(fid,[M N],'single');
  fclose(fid);
end
fb(find(fb<0)) = nan;
fb = squeeze(nanmean(fb,1));
fb(find(isnan(fb)==1)) = -999;
% Write to file:
outfile = char(outfiles(2));
fid = fopen(sprintf('%s.txt',outfile),'w');
fprintf(fid,'%12.4f %12.4f %12.5f %12.5f\n',[y(:)'; x(:)'; fb(:)'; hi(:)']);
fclose(fid);
