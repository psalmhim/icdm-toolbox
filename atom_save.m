% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI", NeuroImage (2026).
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%   3. Neither the name of the copyright holder nor the names of its
%      contributors may be used to endorse or promote products derived
%      from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES.
%
% REQUIRES: SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)

function atom_save(fn, varargin)
% ATOM_SAVE  Atomic, parfor-safe save to MAT file.
%
%   atom_save(fn, 'var1', 'var2', ..., '-v7.3')
%   atom_save(fn, '-struct', 'S', '-v7.3')
%
%   Writes workspace variables to a MAT file via an atomic
%   write-to-temp-then-rename strategy, ensuring that partially written
%   files never appear at the target path.  Safe for use inside parfor
%   loops when each worker targets a unique output path.
%
%   INPUT
%     fn       : target MAT-file path (created if parent directory missing)
%     varargin : variable names to save, OR '-struct','S' for struct mode.
%                Optional save flags (e.g. '-v7.3') may be appended.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also SAVE, ICDM_POPULATION_EB, ICDM_RUN_SUBJECT_VB

    if nargin < 2
        error('atom_save: not enough input arguments.');
    end

    % Ensure parent directory exists
    [pth,~,~] = fileparts(fn);
    if ~isempty(pth) && exist(pth,'dir') ~= 7
        mkdir(pth);
    end

    % Split varargin into mode & flags
    mode_struct = false;
    flags = {};
    varnames = {};
    i = 1;
    while i <= numel(varargin)
        arg = varargin{i};
        if ischar(arg) || isstring(arg)
            arg = char(arg);
            if strcmp(arg, '-struct')
                if i+1 > numel(varargin)
                    error('atom_save: "-struct" must be followed by a struct variable name.');
                end
                mode_struct = true;
                struct_name = char(varargin{i+1});
                i = i + 2;
                continue;
            elseif startsWith(arg, '-') % save flag (e.g., -v7.3)
                flags{end+1} = arg; %#ok<AGROW>
            else
                varnames{end+1} = arg; %#ok<AGROW>
            end
        else
            error('atom_save: unsupported argument type at position %d.', i);
        end
        i = i + 1;
    end

    % Build temp filename in same directory (unique)
    tmp = [tempname(pth) '.part'];

    % Make sure temp is cleaned on error
    c = onCleanup(@() cleanup_temp(tmp));

    % Do the save to temp
    if mode_struct
        S = evalin('caller', struct_name);
        if ~isstruct(S)
            error('atom_save: "%s" is not a struct in caller workspace.', struct_name);
        end
        args = [{tmp, '-struct', 'S'}, flags];
        save(args{:});
    else
        % Collect variables from caller into a struct, then save -struct
        S = struct();
        for k = 1:numel(varnames)
            vn = varnames{k};
            S.(vn) = evalin('caller', vn);
        end
        args = [{tmp, '-struct', 'S'}, flags];
        save(args{:});
    end

    % Atomic rename into place (overwrite if exists)
    move_ok = try_movefile_overwrite(tmp, fn);
    if ~move_ok
        % last resort: delete+move (not strictly atomic, but prevents failures)
        if exist(fn,'file') == 2
            try, delete(fn); catch, end
        end
        move_ok = movefile(tmp, fn);
        if ~move_ok
            error('atom_save: failed to move temp file into place: %s', fn);
        end
    end
end

% -------- helpers --------
function cleanup_temp(tmp)
    if exist(tmp,'file') == 2
        try, delete(tmp); catch, end
    end
end

function ok = try_movefile_overwrite(src, dst)
% Try overwrite-capable move (best effort across MATLAB versions / OSes)
    ok = false;
    try
        % Newer MATLABs: movefile(...,'f') overwrites
        [ok,~] = movefile(src, dst, 'f');
        if ok, return; end
    catch
        % fall through
    end
    % Try vanilla movefile; many OSes allow overwrite by rename
    try
        [ok,~] = movefile(src, dst);
        if ok, return; end
    catch
        % fall through
    end
    % On some systems, an existing dst blocks; try system rename forcibly (UNIX)
    if isunix
        [status,~] = system(sprintf('mv -f "%s" "%s"', src, dst));
        ok = (status == 0);
    end
end
