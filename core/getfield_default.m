function val = getfield_default(S, field, defaultVal)
if isstruct(S) && isfield(S, field), val = S.(field); else, val = defaultVal; end
end