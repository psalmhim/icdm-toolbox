function icdm_remap_dk84(inp, out)
% ICDM_REMAP_DK84  FreeSurfer aparc+aseg -> 84-target label image used for the iCDM counts.
%
%   icdm_remap_dk84(aparc_aseg_nii, out_nii)
%
% Output label value k is target channel k (see target_lookup_84.tsv):
%   1-34  left  Desikan-Killiany cortex  (FreeSurfer 1001-1003, 1005-1035; 1004 corpus callosum excluded)
%   35-68 right Desikan-Killiany cortex  (2001-2003, 2005-2035)
%   69-84 subcortical / cerebellar       (10 11 12 13 17 18 26 49 50 51 52 53 54 58 8 47)
% The 68-target analyses use channels 1-68 (target_lookup_68.tsv).
CORT = [1001 1002 1003 1005:1035 2001 2002 2003 2005:2035];
SUBC = [10 11 12 13 17 18 26 49 50 51 52 53 54 58 8 47];
order = [CORT SUBC];
assert(numel(order) == 84);
V = spm_vol(inp); A = spm_read_vols(V);
O = zeros(size(A));
for k = 1:84, O(round(A) == order(k)) = k; end
Vo = V; Vo.fname = out; Vo.dt = [4 0]; Vo.pinfo = [1;0;0];
spm_write_vol(Vo, O);
end
