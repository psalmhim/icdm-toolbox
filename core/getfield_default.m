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