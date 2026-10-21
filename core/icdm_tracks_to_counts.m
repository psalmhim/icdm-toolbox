function [count_map, out_file, adj_matrix] = icdm_tracks_to_counts(track_file, atlas_file, out_file, opts)
% ICDM_TRACKS_TO_COUNTS  Voxelwise streamline--target incidence counts from a tractogram and a label atlas.
%
%   [COUNT, OUT, ADJ] = icdm_tracks_to_counts(track_file, atlas_file, out_file, opts)
%
% Inputs
%   track_file : MRtrix .tck (world coordinates, mm) or TrackVis .trk (voxel coordinates)
%   atlas_file : integer label image (NIfTI); label value k = target channel k (1..L), 0 = unlabelled.
%                The count image inherits this image's grid and affine.
%   out_file   : 4-D NIfTI written as uint16 [X Y Z L] ('' = do not write)
%   opts       : .labelidx  restrict to these label values (default: all labels > 0)
%
% Counting rule (see docs/CONTRACT.md)
%   1. Each vertex is mapped to the atlas voxel grid: inv(atlas.mat) * [x;y;z;1], rounded to the nearest
%      voxel; vertices outside the grid are dropped.
%   2. If an end vertex lies in an unlabelled voxel, the streamline is extended along its terminal direction
%      by one and two voxels; a labelled voxel reached this way is appended (snap to the adjacent target).
%   3. Consecutive repeats of the same voxel are collapsed to one; a voxel re-entered after leaving it is
%      visited again.
%   4. The streamline's target set is the set of distinct labels over all its voxels (endpoint regions and
%      any labelled structure on its path). A streamline with no label contributes nothing.
%   5. Every visited voxel receives +1 in each channel of the target set (unit weight per streamline).
%   No streamline weights (e.g. SIFT2) are applied.
%
% Outputs
%   count_map  : [X Y Z L] uint16 counts n_{v,k}
%   out_file   : path of the written 4-D image
%   adj_matrix : [L L] target co-occurrence (streamlines whose target set contains both labels)
%
% Requires SPM (spm_vol, spm_read_vols, spm_write_vol).

if nargin < 3, out_file = ''; end
if nargin < 4, opts = struct; end
if ~isfield(opts,'labelidx'), opts.labelidx = []; end
if ~exist(track_file,'file'), error('icdm_tracks_to_counts: track file not found: %s', track_file); end

V       = spm_vol(atlas_file);
atlas3d = uint16(round(spm_read_vols(V(1))));
if ~isempty(opts.labelidx)
    atlas3d(~ismember(atlas3d, opts.labelidx)) = 0;
end
dim   = V(1).dim(1:3);
Lmax  = double(max(atlas3d(:)));
if Lmax < 1, error('icdm_tracks_to_counts: atlas has no labels'); end
M_vox = inv(V(1).mat);

count_map  = zeros([dim Lmax], 'uint16');
adj_matrix = zeros(Lmax, Lmax, 'uint32');
total_fibres = 0;

[~,~,ext] = fileparts(lower(track_file));
switch ext
    case '.trk', voxel_space = true;  trk_reader(@process_fibre);
    case '.tck', voxel_space = false; tck_reader(@process_fibre);
    otherwise,   error('icdm_tracks_to_counts: unsupported extension "%s"', ext);
end
fprintf('[icdm_tracks_to_counts] %d streamlines processed\n', total_fibres);

if ~isempty(out_file)
    if exist(out_file,'file'), delete(out_file); end
    for l = 1:Lmax
        Vo = V(1); Vo.fname = out_file; Vo.dt = [spm_type('uint16') 0];
        Vo.pinfo = [1;0;0]; Vo.n = [l 1];
        spm_write_vol(Vo, double(count_map(:,:,:,l)));
    end
end

    function tf = in_bounds(ii)
        tf = all(ii >= 1) && ii(1) <= dim(1) && ii(2) <= dim(2) && ii(3) <= dim(3);
    end

    function ijk = snap_end(ijk, at_head)
        % extend an unlabelled end by one and two voxels along the terminal direction
        if at_head, e = 1; nb = 2; nb2 = 3; else, e = size(ijk,2); nb = e-1; nb2 = e-2; end
        if size(ijk,2) < 2 || atlas3d(ijk(1,e),ijk(2,e),ijk(3,e)) ~= 0, return; end
        p1 = ijk(:,e); v = p1 - ijk(:,nb); n = norm(v);
        if n < eps
            if size(ijk,2) < 3, return; end
            v = p1 - ijk(:,nb2); n = norm(v);
            if n < eps, return; end
        end
        u = v / n;
        lab = 0; p2 = round(p1 + u);
        if in_bounds(p2)
            lab = atlas3d(p2(1),p2(2),p2(3));
            if lab > 0, ijk = add(ijk, p2, at_head); end
        end
        plab = lab; p2 = round(p1 + 2*u);
        if in_bounds(p2)
            lab = atlas3d(p2(1),p2(2),p2(3));
            if plab == 0, plab = lab; end
            if lab > 0 && plab == lab, ijk = add(ijk, p2, at_head); end
        end
    end

    function ijk = add(ijk, p, at_head)
        if at_head, ijk = [p ijk]; else, ijk = [ijk p]; end
    end

    function process_fibre(XYZ)
        % mapping is done in single precision, as in the production counter
        if voxel_space
            ijk = round(XYZ);
        else
            ijk = M_vox * [XYZ; ones(1,size(XYZ,2))];
            ijk = round(ijk(1:3,:));
        end
        ijk = ijk(:, all(ijk >= 1 & ijk <= dim(:), 1));
        if isempty(ijk), return; end
        ijk = snap_end(ijk, true);
        ijk = snap_end(ijk, false);
        dup = all(diff(ijk,1,2) == 0, 1);
        ijk = ijk(:, [true, ~dup]);
        lin  = sub2ind(dim, ijk(1,:), ijk(2,:), ijk(3,:));
        labs = atlas3d(lin);
        labs = unique(labs(labs > 0 & labs <= Lmax));
        if isempty(labs), return; end
        for a = 1:numel(labs)
            for b = a+1:numel(labs)
                adj_matrix(labs(a),labs(b)) = adj_matrix(labs(a),labs(b)) + 1;
                adj_matrix(labs(b),labs(a)) = adj_matrix(labs(b),labs(a)) + 1;
            end
        end
        for lab = labs(:)'
            for p = 1:size(ijk,2)
                count_map(ijk(1,p),ijk(2,p),ijk(3,p),lab) = count_map(ijk(1,p),ijk(2,p),ijk(3,p),lab) + 1;
            end
        end
    end

    function trk_reader(cb)
        fid = fopen(track_file,'rb','ieee-le');
        if fid < 0, error('Cannot open %s', track_file); end
        fread(fid, 1000, '*uint8');
        while ~feof(fid)
            np = fread(fid,1,'int32');
            if isempty(np) || np <= 0, break; end
            xyz = fread(fid,[3 np],'float32=>single');
            if numel(xyz) ~= 3*np, break; end
            cb(xyz); total_fibres = total_fibres + 1;
        end
        fclose(fid);
    end

    function tck_reader(cb)
        fid = fopen(track_file,'rb','ieee-le');
        if fid < 0, error('Cannot open %s', track_file); end
        n_scalars = 0; data_file = '.'; data_off = 0;
        while true
            line = fgetl(fid);
            if ~ischar(line), error('Unexpected end of .tck header'); end
            if isempty(line) || strcmp(strtrim(line),'END'), break; end
            if startsWith(line,'scalar_name')
                n_scalars = n_scalars + 1;
            elseif startsWith(line,'property_count')
                n_scalars = max(n_scalars, sscanf(line,'property_count: %d'));
            elseif startsWith(line,'file:')
                C = textscan(line,'file: %s %d'); data_file = strtrim(C{1}{1}); data_off = C{2};
            end
        end
        stride = 3 + n_scalars;
        if ~strcmp(data_file,'.')
            fclose(fid);
            fid = fopen(fullfile(fileparts(track_file), data_file),'rb','ieee-le');
            if fid < 0, error('Cannot open payload file %s', data_file); end
        end
        fseek(fid, data_off, 'bof');
        chunk = 2e6; tail = single([]); done = false;
        while ~done
            buf = fread(fid, chunk, 'float32=>single');
            if isempty(buf), done = true; end
            vec = [tail; buf];
            infpos = find(isinf(vec), 1);
            if ~isempty(infpos), vec = vec(1:infpos-1); done = true; end
            isn = isnan(vec); d = diff([0; isn; 0]);
            st = find(d == 1); en = find(d == -1) - 1;
            sep = st((en - st + 1) >= stride);
            start = 1;
            for s = sep(:)'
                stop = s - 1;
                if stop >= start
                    block = vec(start:stop);
                    if mod(numel(block), stride) == 0
                        xyzs = reshape(block, stride, []);
                        cb(single(xyzs(1:3,:))); total_fibres = total_fibres + 1;
                    else
                        warning('icdm_tracks_to_counts: discarded malformed streamline (%d values)', numel(block));
                    end
                end
                start = s + stride;
            end
            tail = vec(start:end);
        end
        fclose(fid);
    end
end
