function [m] = data2niivol(values, blacklistFileName, roisFileName, niiOutFileName, extdataPath)

% path(path, '/m/nbe/scratch/braindata/shared/toolboxes/NIFTI')
% path(path, '/proj/networks/rmkujala/brain_fmri/code/src')
% path(path, '/proj/braindata/eglerean/net_viz')

%might be unnecessary

values = reshape(values, numel(values), 1);

load(roisFileName);
load(blacklistFileName);
psess.rois=rois;        % your rois
psess.datasize=[91 109 91]; % the size of your data

% a list of nodes that are not good because for example
% they are too close to the edges

%num2str(length(blacklist)))
psess.blacklist = find(blacklist < 0.5);
% 0.5 = hack

cfg=[];
cfg.psess=psess;
cfg.vals=values;
cfg.filename=niiOutFileName;
cfg.box=3;  % if == 0, then no smoothing


funpsy_interpvol(cfg);
nii = load_nii(cfg.filename);
sum(nii.img(find(nii.img < 0) ) );


nii=fixOriginator(cfg.filename, extdataPath);
nii.img(find(isnan(nii.img)))=0;

size(nii.img);
max(max(max(nii.img)));
save_nii(nii,cfg.filename);

m = 0; % just to return something
