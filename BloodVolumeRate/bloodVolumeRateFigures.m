function bloodVolumeRateFigures(Q_results, mask, name, M0_ff_video, xy_barycenter, systolesIndexes, sysIdx, diasIdx)

% 0. Initialise Variables

ToolBox = getGlobalToolBox;
path_png = ToolBox.path_png;
path_eps = ToolBox.path_eps;

if ~isfolder(fullfile(path_png, 'volumeRate'))
    mkdir(path_png, 'volumeRate')
    mkdir(path_eps, 'volumeRate')
end

params = ToolBox.getParams;

initial = name(1);
M0_ff_video = rescale(M0_ff_video);
M0_ff_img = rescale(mean(M0_ff_video, 3));
[numX, ~, numFrames] = size(M0_ff_video);
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

numSections = Q_results.numSections;
locs = Q_results.locs;
v_cell = Q_results.v_cell;
v_profiles_cell = Q_results.v_profiles_cell;
dv_profiles_cell = Q_results.dv_profiles_cell;
D_cell = Q_results.D_cell;
dD_cell = Q_results.dD_cell;
mask_mat = Q_results.mask_mat;
subImg_cell = Q_results.subImg_cell;
area_mat = Q_results.area_mat;
Q_mat = Q_results.Q_mat;
dQ_mat = Q_results.dQ_mat;

% 1. Sections Image
tic

if ~isempty(systolesIndexes)
    index_start = systolesIndexes(1);
    index_end = systolesIndexes(end);
else
    index_start = 1;
    index_end = numFrames;
end

try

    if params.json.BloodVolumeRateFigures.sectionImage
        sectionImage(M0_ff_img, mask_mat, initial)
    end

catch ME
    MEdisp(ME, ToolBox.path_dir)
end

subImageSize = checkSubImgSize(subImg_cell);

try

    if params.json.BloodVolumeRateFigures.widthImage
        widthImage(subImageSize, subImg_cell, numSections, name)
    end

catch ME
    MEdisp(ME, ToolBox.path_dir)
end

try

    if params.json.BloodVolumeRateFigures.crossSectionImages

        if ~isfolder(fullfile(ToolBox.path_png, 'volumeRate', 'sectionsImages'))
            mkdir(fullfile(path_png, 'volumeRate'), 'sectionsImages')
            mkdir(fullfile(path_eps, 'volumeRate'), 'sectionsImages')
            mkdir(fullfile(path_png, 'volumeRate', 'sectionsImages'), 'widths')
            mkdir(fullfile(path_eps, 'volumeRate', 'sectionsImages'), 'widths')
            mkdir(fullfile(path_png, 'volumeRate', 'sectionsImages'), 'num')
            mkdir(fullfile(path_eps, 'volumeRate', 'sectionsImages'), 'num')
            mkdir(fullfile(path_png, 'volumeRate', 'sectionsImages'), 'bvr')
            mkdir(fullfile(path_eps, 'volumeRate', 'sectionsImages'), 'bvr')
            mkdir(fullfile(path_png, 'volumeRate', 'sectionsImages'), 'vel')
            mkdir(fullfile(path_eps, 'volumeRate', 'sectionsImages'), 'vel')
        end

        crossSectionImages(M0_ff_img, xy_barycenter, area_mat, Q_mat, v_cell, mask_mat, locs, name)
    end

catch ME
    MEdisp(ME, ToolBox.path_dir)
end

try

    if params.json.BloodVolumeRateFigures.widthHistogram
        widthHistogram(D_cell, dD_cell, area_mat, name);
    end

catch ME
    MEdisp(ME, ToolBox.path_dir)
end

fprintf("    1. Sections Images Generation (%s) took %ds\n", name, round(toc))

% 2. Blood Volume Rate Figures
tic

numCircles = params.json.BloodVolumeRateAnalysis.NumberOfCircles;
px_size = params.px_size;
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
dr = (r2 - r1) / numCircles;
rad = linspace(r1, r2 - dr, numCircles) * numX / px_size;

[Q_t, dQ_t] = plotRadius(Q_mat, dQ_mat, t, rad, index_start, index_end, name);

try

    if params.json.BloodVolumeRateFigures.BloodFlowProfiles

        if ~isfolder(fullfile(ToolBox.path_png, 'volumeRate', 'velocityProfiles'))
            mkdir(fullfile(ToolBox.path_png, 'volumeRate', 'velocityProfiles'));
        end

        interpolatedBloodVelocityProfile(v_profiles_cell, dv_profiles_cell, sysIdx, diasIdx, numSections, name, rad)
    end

catch ME
    MEdisp(ME, ToolBox.path_dir)
end

% Call for arterial analysis
graphCombined(M0_ff_video, imdilate(mask, strel('disk', params.json.PulseAnalysis.LocalBackgroundWidth)), ...
    Q_t, dQ_t, xy_barycenter, sprintf('%s_vr', name), ...
    'etiquettes_locs', [], ...
    'etiquettes_values', [], ...
    'ylabl', 'Volume Rate (µL/min)', ...
    'xlabl', 'Time (s)', ...
    'fig_title', 'Blood Volume Rate', ...
    'unit', 'µL/min', ...
    'skip', ~params.exportVideos, ...
    'Color', name, ...
    'Visible', false);

fprintf("    2. Blood Volume Rate Figures (%s) took %ds\n", name, round(toc))

% 3. Arterial Indicators
tic

if strcmp(name, 'Artery')

    try

        if params.json.BloodVolumeRateFigures.ARIBVR
            ArterialResistivityIndex(t, Q_t, mask, 'BVR', 'volumeRate');
        end

    catch ME
        MEdisp(ME, ToolBox.path_dir)
    end

    try

        if params.json.BloodVolumeRateFigures.strokeAndTotalVolume && ~isempty(systolesIndexes)
            strokeAndTotalVolume(Q_t, dQ_t, systolesIndexes, t, 1000);
        end

    catch ME
        MEdisp(ME, ToolBox.path_dir)
    end

    fprintf("    3. Arterial Indicators Images Generation (%s) took %ds\n", name, round(toc))
end

end
