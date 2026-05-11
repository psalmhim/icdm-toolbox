function subj = icdm_compose_subject(icdm, property_names, property_values, name, dartel_template, outpath, opts)
% ICDM_COMPOSE_SUBJECT  Create and cache a subject struct for the iCDM pipeline.
%
%   subj = icdm_compose_subject(icdm, property_names, property_values,
%                                name, dartel_template, outpath, opts)
%
%   Assembles the subject-level data structure used throughout the iCDM
%   framework.  Demographic covariates are stored as arbitrary key-value
%   pairs in subj.property, the DARTEL deformation field is located
%   automatically, and a cached data file (icdm_subj_<name>.mat) is
%   created containing the WM mask, native indices, and prebuilt
%   trilinear warp operators.
%
%   Inputs
%     icdm             : path to the 4-D ICDM NIfTI (native space).
%     property_names   : cell array of covariate names, e.g. {'age','sex'}.
%     property_values  : numeric vector of matching covariate values.
%     name             : subject identifier string.
%     dartel_template  : path to the DARTEL template NIfTI.
%     outpath          : directory for cached output files (default: same
%                        as icdm file).
%     opts             : optional struct; opts.thresh sets the WM mask
%                        threshold (default 5).
%
%   Output
%     subj : struct with fields id, icdm_4d, template, property,
%            dartel_flow, outpath, datafile.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_BUILD_SUBJECT_WARPS, ICDM_RUN_SUBJECT_VB

if nargin < 7, opts = struct(); end
if nargin < 6, outpath = ''; end
if nargin < 5, dartel_template = ''; end
if nargin < 4, name = ''; end

if ~exist(icdm,'file')
    error('ICDM file not found: %s', icdm);
end

% -------------------------------------------------------------------------
% 0) Create subj struct
% -------------------------------------------------------------------------
subj.id = name;
subj.icdm_4d = icdm;
subj.template = dartel_template;

% property_names must match property_values
if ~iscell(property_names) || ~isvector(property_values)
    error('property_names must be cell array, property_values must be vector/array.');
end
if numel(property_names) ~= numel(property_values)
    error('property_names and property_values must have same length.');
end

% Create subj.property.* dynamically
subj.property = struct();
for i = 1:numel(property_names)
    pname = property_names{i};
    pval  = property_values(i);
    subj.property.(pname) = pval;
end

% Output path
if isempty(outpath)
    subj.outpath = fileparts(icdm);
else
    subj.outpath = outpath;
end

% -------------------------------------------------------------------------
% 1) Attach DARTEL deformation field
% -------------------------------------------------------------------------
flow_y = fullfile(fileparts(icdm),'dartel','y_fs_t1.nii');
if exist(flow_y,'file')
    subj.dartel_flow = flow_y;
else
    warning('Missing DARTEL flow field for %s', subj.id);
    subj.dartel_flow = '';
end

% -------------------------------------------------------------------------
% 2) subj info file: mask_3d, idx_native, dim_native
% -------------------------------------------------------------------------

datafile = fullfile(subj.outpath, sprintf('icdm_subj_%s.mat', subj.id));
subj.datafile = datafile;

if exist(datafile,'file')
    subj.datafile= datafile;
else
    % Read ICDM 4D
    V  = spm_vol(subj.icdm_4d);
    C4 = spm_read_vols(V);
    thr = getfield_default(opts,'thresh',5);
    sumC = sum(C4,4);
    mask_3d = sumC > thr;
    idx_native = find(mask_3d(:));
    icdm2d = reshape(C4, [], size(C4,4));
    icdm2d = icdm2d(idx_native, :);
    dim_native = V(1).dim;
    V=V(1);

    warp.to_native = icdm_build_warp_to_native( ...
            subj.icdm_4d, subj.dartel_flow,subj.template);

    warp.to_mni = icdm_build_warp_to_mni( ...
            subj.icdm_4d, subj.dartel_flow, subj.template);
    save(datafile, ...
        'icdm2d','V','idx_native','dim_native','warp','thr','-v7.3');
end

% % -------------------------------------------------------------------------
% % 3) warp operator
% % -------------------------------------------------------------------------
% warpfile = fullfile(subj.outpath, sprintf('warp_subj_%s.mat', subj.id));
% subj.warpfile = warpfile;

% if exist(warpfile,'file')
%     fprintf('[ICDM] Loading warp file: %s\n', warpfile);

%     W = load(warpfile,'warp_native','warp_mni');
%     subj.warp.to_native = W.warp_native;
%     subj.warp.to_mni    = W.warp_mni;

% else
%     fprintf('[ICDM] Creating warp file: %s\n', warpfile);

%     if ~isfield(opts,'dim_mni')
%         error('opts.dim_mni is required to build warp operator.');
%     end

%     % Build warp operators
%     subj.warp.to_native = build_warp_to_native( ...
%         subj.dartel_flow, subj.template, subj.icdm_4d, 50000);

%     subj.warp.to_mni = build_warp_to_mni( ...
%         subj.dartel_flow, subj.template, opts.dim_mni, 50000);

%     % Save warp operator
%     warp_native = subj.warp.to_native;
%     warp_mni    = subj.warp.to_mni;

%     save(warpfile, 'warp_native','warp_mni', '-v7.3');
% end

end
