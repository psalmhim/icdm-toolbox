% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI" (submitted).
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

function val = getfield_default(S, field, defaultVal)
% GETFIELD_DEFAULT  Retrieve a struct field with a fallback default value.
%
%   val = getfield_default(S, field, defaultVal)
%
%   Returns S.(field) if S is a struct and the field exists; otherwise
%   returns defaultVal.  Used throughout the iCDM toolbox to provide safe
%   access to optional configuration parameters.
%
%   Inputs
%     S          : scalar struct (or non-struct, in which case default is used).
%     field      : char field name to look up.
%     defaultVal : value returned when the field is absent.
%
%   Output
%     val : the retrieved or default value.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_SUBJECT_VB, ICDM_ESTIMATE_BETA
if isstruct(S) && isfield(S, field), val = S.(field); else, val = defaultVal; end
end