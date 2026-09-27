function x = icdm_design_vector(subj, opts)

spec = opts.design_spec;
x_list = {};

for i = 1:numel(spec)
    entry = spec{i};

    % ------------------ constant --------------------
    if ischar(entry) && strcmpi(entry,'const')
        x_list{end+1} = 1;
        continue;
    end

    % ------------------ normal variable --------------------
    if ischar(entry)
        name = entry;
        val = subj.property.(name);

        if strcmp(opts.design.type.(name),'continuous')
            val = (val - opts.design.mu.(name)) / opts.design.sd.(name);
        end

        x_list{end+1} = val;
        continue;
    end

    % ------------------ spline struct --------------------
    if isstruct(entry)
        basevar = entry.var;
        Kspl    = entry.K;

        fname = sprintf('%s_spline_%d', basevar, Kspl);
        fname = matlab.lang.makeValidName(fname);

        % z-scored base value
        v = subj.property.(basevar);
        v = (v - opts.design.mu.(basevar)) / opts.design.sd.(basevar);

        % build natural cubic spline basis
        bx = icdm_make_spline_basis(v, Kspl);

        x_list = [x_list, num2cell(bx)]; 
        continue;
    end
end

x = cell2mat(x_list);

end

function b = icdm_make_spline_basis(x, K)
% Natural cubic spline basis with K knots on [-2,2] (z-scored domain)
knots = linspace(-2,2,K);
b = zeros(1,K);
for i = 1:K
    b(i) = max(0, x - knots(i)).^3;
end
end
