function make_synthetic_visits(root)
% MAKE_SYNTHETIC_VISITS  Small synthetic longitudinal infant dataset laid out as unc4d_config.json expects.
%
% 5 participants x 3 visits (about 3, 6 and 12 months). Participants 1-4 form the reference sample;
% participant 5 is analysed under their prior without contributing to it.
% Writes under root:
%   data/sessions.tsv                                  participant_id, session_id, age_months, reference
%   data/sub-XX/ses-YY/dwi/{tracts.tck, atlas.nii}     tractogram (world mm) and target label image
%   data/sub-XX/ses-YY/xfm/y_to_monthNN.nii            subject -> age-specific template deformation
%   templates/monthNN.nii, templates/common.nii        stand-ins for the UNC-4D templates
%   templates/y_common_to_monthNN.nii                  common space -> age-specific template deformation
% All deformations are identities here; with real data they come from registration (e.g. SPM, ANTs
% converted to SPM y_*.nii format).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'minimal_synthetic'));
P = 5; months = [3 6 12]; names = {'month03','month06','month12'};
rng(3);
age_m = repmat(months, P, 1) + 0.3*randn(P, 3);              % participant x visit
gen = fullfile(root,'tmp_generate');
make_synthetic_cohort(gen, P*3, 2000, 3, 8 + 4*(reshape(age_m',[],1) - 3)/9);
td = fullfile(root,'templates'); if ~exist(td,'dir'), mkdir(td); end
copyfile(fullfile(gen,'template.nii'), fullfile(td,'common.nii'));
for b = 1:3
    copyfile(fullfile(gen,'template.nii'), fullfile(td,[names{b} '.nii']));
    make_identity_deformation(fullfile(td,'common.nii'), fullfile(td,['y_common_to_' names{b} '.nii']));
end
dd = fullfile(root,'data'); if ~exist(dd,'dir'), mkdir(dd); end
fid = fopen(fullfile(dd,'sessions.tsv'),'w');
fprintf(fid,'participant_id\tsession_id\tage_months\treference\n');
k = 0;
for p = 1:P
    for v = 1:3
        k = k + 1;
        sd = fullfile(dd, sprintf('sub-%02d',p), sprintf('ses-%02d',v));
        mkdir(fullfile(sd,'dwi')); mkdir(fullfile(sd,'xfm'));
        copyfile(fullfile(gen, sprintf('sub-%02d',k), 'tracts.tck'), fullfile(sd,'dwi','tracts.tck'));
        copyfile(fullfile(gen,'atlas.nii'), fullfile(sd,'dwi','atlas.nii'));
        make_identity_deformation(fullfile(td,[names{v} '.nii']), fullfile(sd,'xfm',['y_to_' names{v} '.nii']));
        fprintf(fid,'sub-%02d\tses-%02d\t%.1f\t%d\n', p, v, age_m(p,v), p <= 4);
    end
end
fclose(fid);
rmdir(gen, 's');
end
