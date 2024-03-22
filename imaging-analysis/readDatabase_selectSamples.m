%% readDatabase_selectSamples.m
% 
% Script to select indices of time series from a .txt spreadsheet based on 
%  specified metadata values 
% Creates a struct, metaDat, each field is a column from the spreadsheet
% All of the i_ variables are logicals for whether particular metadata
%  values are met. Used to select particular subset of time series.
%
% Users: filename should point to your spreadsheet of info on each time
% series. Add i_ variables in same format for new conditions.
%
% 3/8/17 - HHY - Updated to include wavelength and ditherHold
%
% filename = '/Users/hyang/Documents/New Imaging Analysis/2p_imaging_log.txt';
formatSpec = '%s %f %s %f %s %s %f %f %f %f %f %f %f %f %f %s %f %f %f %f %f %s\r';

[seriesID, flyID, genotype, age, sex, stimcode, imFrames, rot, zdepth, ...
    layer, estFrameRate, wavelength, transmission, gain, offset, ...
    colorfilter, ND, pwm, ditherHold, isMoving, isResponding, comments] ...
    = textread(filename, formatSpec, 'headerlines', 1, 'delimiter', '\t');

metaDat.seriesID = seriesID;
metaDat.flyID = flyID;
metaDat.genotype = genotype;
metaDat.tissueLayer = layer;
metaDat.age = age;
metaDat.sex = sex;
metaDat.stimcode = stimcode;
metaDat.rotation = rot;
metaDat.zdepth = zdepth;
metaDat.wavelength = wavelength;
metaDat.estFrameRate = estFrameRate;
metaDat.colorfilter = colorfilter;
metaDat.NDfilter = ND; 
metaDat.pwm = pwm;
metaDat.ditherHold = ditherHold;
metaDat.isMoving = isMoving;

%% Also Works!
% fid = fopen(filename);
% D = textscan(fid, ... 
%     '%s %f %s %f %s %s %f %f %f %s %f %f %f %f %s %f',...
%     'headerlines', 1, 'delimiter','\t');
% fclose(fid);

%% 

i_ND = ND;
i_482_18 = strcmpi('482/18',colorfilter);
i_447_60 = strcmpi('447/60',colorfilter);

i_pwm = pwm;

% sparsen by fly ID
i_flyID = flyID;
i_male = strcmpi('M', sex);
i_female = strcmpi('F',sex);

i_300msSearchStimFlash = strcmpi('300ms_searchStimFlash', stimcode);
i_300msSearchStimFlashRGB = strcmpi('300ms_searchStimFlashRGB', stimcode);
i_2s_searchStimFlash = strcmpi('2s_searchStimFlash', stimcode);
i_2sFFF = strcmpi('2sFFF_template', stimcode);
i_300msFFF = strcmpi('300msFFF_binarycontrast', stimcode);
i_fullfield_LDflash20ms200ms_Gray750ms = strcmpi(...
    'fullfield_LDflash20ms200ms_Gray750ms', stimcode);
i_fullfield_LDflash20ms_Gray500ms = strcmpi(...
    'fullfield_LDflash20ms_Gray500ms', stimcode);
i_fullfield_LDflash20ms_c025Gray500ms = strcmpi(...
    'fullfield_LDflash20ms_c025Gray500ms', stimcode);
i_fullfield_LDflash10ms_Gray500ms = strcmpi(...
    'fullfield_LDflash10ms_Gray500ms', stimcode);
i_fullfield_LDflash30ms_Gray500ms = strcmpi(...
    'fullfield_LDflash30ms_Gray500ms', stimcode);
i_fullfield_LDflash40ms_Gray500ms = strcmpi(...
    'fullfield_LDflash40ms_Gray500ms', stimcode);
i_fullfield_Lc075c1flash20ms_c05Gray500ms = strcmpi(...
    'fullfield_L075-1flash20ms_05Gray500ms', stimcode);
i_fullfield_Lc05c075c1flash20ms_c025Gray500ms = strcmpi(...
    'fullfield_L05-075-1flash20ms_025Gray500ms', stimcode);
i_fullfield_Dc0c0125c025flash20ms_c05Gray500ms = strcmpi(...
    'fullfield_D0-0125-025flash20ms_05Gray500ms', stimcode);
i_fullfield_Dc0c0125flash20ms_c025Gray500ms = strcmpi(...
    'fullfield_D0-0125flash20ms_025Gray500ms', stimcode);
i_whitenoise_1D_rand9_20ms = strcmpi('whitenoise_1D_rand9_20ms', stimcode);
i_whitenoise_1D_rand9_50ms = strcmpi('whitenoise_1D_rand9_50ms', stimcode);
i_2s_2epoch_NatStim = strcmpi('2s_2epoch_NatStim', stimcode);
i_2s_1epoch2_NatStim = strcmpi('2s_1epoch-2_NatStim', stimcode);

i_movingbar_dark_c1_5deg_4dir_20dps = strcmpi(...
    'movingBar_darkC1_5deg_4dir_20ds', stimcode);
i_movingedge_4dir_20dps = strcmpi('movingEdge_4dir_20ds',stimcode);
i_staticmovinggrating_20deg_20dps_c1_12dir = strcmpi(...
    'staticMovingGrating_20deg_20ds_c1_12dir',stimcode);
i_movingbar_light_c1_5deg_4dir_20dps = strcmpi(...
    'movingBar_lightC1_5deg_4dir_20ds', stimcode);

% i_fff  = strcmpi('fff',(stimcode)); % not including use of ND
% i_fff2s3s = ~(cellfun(@isempty,regexpi(stimcode,'fff2s3s'))); % including use of arbitrary ND
% i_fff2s3sND0 = strcmpi('fff2s3s',(stimcode)); % not including use of ND
% i_fff4  = strcmpi('fff4s',(stimcode)); % not including use of ND
% i_fff3s_gray  = strcmpi('fff3s_gray',(stimcode)); % not including use of ND
% i_fffvar = strcmpi('fffvar',(stimcode));
% i_fffvarrand = strcmpi('fffvarrand',(stimcode));
% i_fff_partialC1 = strcmpi('fff_partialC1',(stimcode));
% i_fff_partialC2 = strcmpi('fff_partialC2',(stimcode));
% i_fff_partialC2per = strcmpi('fff_partialC2per',(stimcode));
% i_fff_partialC2per_brightgray = strcmpi('fff_partialC2per_brightgray',(stimcode));
% i_fff_partialC2per_darkgray = strcmpi('fff_partialC2per_darkgray',(stimcode));
% i_fff_partialC2per_darkgray_exp = strcmpi('fff_partialC2per_darkgray_exp',(stimcode));
% i_fff_partialC2per_brightgray_exp = strcmpi('fff_partialC2per_brightgray_exp',(stimcode));
% i_fff_partialC2per_darkgray_exp_rev = strcmpi('fff_partialC2per_darkgray_exp_rev',(stimcode));
% i_fff_partialC2per_brightgray_exp_rev = strcmpi('fff_partialC2per_brightgray_exp_rev',(stimcode));
% i_fff_partialC2per_darkgray_exp_fastCR = strcmpi('fff_partialC2per_darkgray_exp_fastCR',(stimcode));
% i_fff_partialC2per_brightgray_exp_fastCR = strcmpi('fff_partialC2per_brightgray_exp_fastCR',(stimcode));
% i_fff_partialC2per_darkgray_exp_revfastCR = strcmpi('fff_partialC2per_darkgray_exp_revfastCR',(stimcode));
% i_fff_partialC2per_brightgray_exp_revfastCR = strcmpi('fff_partialC2per_brightgray_exp_revfastCR',(stimcode));
% i_fff_partialC2per_p1_gray = strcmpi('fffpartial_p1_gray',(stimcode));
% i_fff_partialC2per_p2_gray = strcmpi('fffpartial_p2_gray',(stimcode));
% i_fff_partialC2per_p05_gray = strcmpi('fffpartial_p05_gray',(stimcode));
% i_fff_partialC2per_p025_gray = strcmpi('fffpartial_p025_gray',(stimcode));
% i_fff_partialC2per_p3_gray = strcmpi('fffpartial_p3_gray',(stimcode));
% i_fff_partialC2per_p10_gray = strcmpi('fffpartial_p10_gray',(stimcode)); % 1 second long pulses
% 
% i_circtor = ~(cellfun(@isempty,regexpi(stimcode,'circtor_cell')));
% i_circtor_gray = ~(cellfun(@isempty,regexpi(stimcode,'circtor_gray_cell')));
% i_circtor_graydark = ~(cellfun(@isempty,regexpi(stimcode,'circtor_graydark_LB_cell')));
% i_circtor_graybright = ~(cellfun(@isempty,regexpi(stimcode,'circtor_graybright_LB_cell')));
% i_circtor_morph_graydark = ~(cellfun(@isempty,regexpi(stimcode,'circtor_morph_graydark_cell')));
% i_circtor_morph_graybright = ~(cellfun(@isempty,regexpi(stimcode,'circtor_morph_graybright_cell')));
% i_circtor_morph_LB2_graydark = ~(cellfun(@isempty,regexpi(stimcode,'circtor_morph_graydark_LB2_cell')));
% i_circtor_morph_LB2_graybright = ~(cellfun(@isempty,regexpi(stimcode,'circtor_morph_graybright_LB2_cell')));
% i_circtor_morph_LB2_simul = ~(cellfun(@isempty,regexpi(stimcode,'circtor_morph_simul_LB2_cell')));
% i_circtor_morph_LB2_p2 = ~(cellfun(@isempty,regexpi(stimcode,'circtor_morph_simul_p2_LB2_cell')));
% 
% i_singbar_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_x_cell')));
% i_singbar_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_y_cell')));
% i_singbar_offset10_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_offset10_x_cell')));
% i_singbar_offset10_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_offset10_y_cell')));
% i_singbar_offset15_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_offset15_x_cell')));
% i_singbar_offset15_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_offset15_y_cell')));
% i_singbar_graydark_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graydark_x_cell')));
% i_singbar_graydark_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graydark_y_cell')));
% i_singbar_graydark_offset10_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graydark_offset10_x_cell')));
% i_singbar_graydark_offset10_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graydark_offset10_y_cell')));
% i_singbar_graydark_offset15_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graydark_offset15_x_cell')));
% i_singbar_graydark_offset15_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graydark_offset15_y_cell')));
% i_singbar_graydark_offset25_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graydark_offset25_x_cell')));
% i_singbar_graydark_offset25_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graydark_offset25_y_cell')));
% i_singbar_graydark_centered_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graydark_centered_x_cell')));
% i_singbar_graydark_centered_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graydark_centered_y_cell')));
% i_singbar_graybright_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graybright_x_cell')));
% i_singbar_graybright_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graybright_y_cell')));
% i_singbar_graybright_offset10_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graybright_offset10_x_cell')));
% i_singbar_graybright_offset10_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graybright_offset10_y_cell')));
% i_singbar_graybright_offset15_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graybright_offset15_x_cell')));
% i_singbar_graybright_offset15_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graybright_offset15_y_cell')));
% i_singbar_graybright_offset25_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graybright_offset25_x_cell')));
% i_singbar_graybright_offset25_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graybright_offset25_y_cell')));
% i_singbar_graybright_centered_x = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graybright_centered_x_cell')));
% i_singbar_graybright_centered_y = ~(cellfun(@isempty,regexpi(stimcode,'singbar_graybright_centered_y_cell')));
% i_barpair_x_dark_phi = ~(cellfun(@isempty,regexpi(stimcode,'barpair_graydark_x_phi_cell')));
% i_barpair_x_bright_phi = ~(cellfun(@isempty,regexpi(stimcode,'barpair_graybright_x_phi_cell')));
% i_barpair_x_bright_revphi = ~(cellfun(@isempty,regexpi(stimcode,'barpair_graybright_x_revphi_cell')));
% i_barpair_x_dark_revphi = ~(cellfun(@isempty,regexpi(stimcode,'barpair_graydark_x_revphi_cell')));
% i_barpair_y_dark_phi = ~(cellfun(@isempty,regexpi(stimcode,'barpair_graydark_y_phi_cell')));
% i_barpair_y_bright_phi = ~(cellfun(@isempty,regexpi(stimcode,'barpair_graybright_y_phi_cell')));
% i_barpair_y_bright_revphi = ~(cellfun(@isempty,regexpi(stimcode,'barpair_graybright_y_revphi_cell')));
% i_barpair_y_dark_revphi = ~(cellfun(@isempty,regexpi(stimcode,'barpair_graydark_y_revphi_cell')));
% 
% i_barpair_gray_rotalln = ~(cellfun(@isempty,regexpi(stimcode,'barpair_gray_rotalln')));
% 
% i_ND1 = ~(cellfun(@isempty,regexpi(stimcode,'ND1')));
% i_ND2 = ~(cellfun(@isempty,regexpi(stimcode,'ND2')));
% i_ND05 = ~(cellfun(@isempty,regexpi(stimcode,'ND0.5')));
% 
% i_mb = strcmpi('movebar',(stimcode));
% i_mb_rperm = strcmpi('movebar_rperm',(stimcode));
% i_mb_rpermcsl = strcmpi('movebar_rpermcsl',(stimcode));
% i_mb_cylstripe = strcmpi('movecyl_stripe',(stimcode));
% i_mb_widths = strcmpi('movebar_widths',(stimcode));
% i_mb_widthsi = strcmpi('movebar_widthsi',(stimcode));
% i_mb_nwidthsrperm = strcmpi('movebar_nwidthsrperm',(stimcode));
% i_mb_vnwidthsrperm = strcmpi('movebar_vertnwidthsrperm',(stimcode));
% 
% i_rf400 = strcmpi('rf400',(stimcode)); % will have to differentiate on date, eventually
% i_rf100n = strcmpi('rf100n',(stimcode)); % will have to differentiate on date, eventually
% 
% i_edges = strcmpi('edges2',stimcode);
% i_edges_rperm = ~(cellfun(@isempty,regexpi(stimcode,'edges2_rperm'))); % including use of arbitrary ND
% 
% i_gauss200 = strcmpi('gaussian200',stimcode);
% i_gauss500 = strcmpi('gaussian500',stimcode);
% 
% i_narrowV = ~(cellfun(@isempty,regexpi(stimcode,'vertgridnarrow'))); % including use of arbitrary ND
% i_narrowH = ~(cellfun(@isempty,regexpi(stimcode,'horizgridnarrow'))); % including use of arbitrary ND
% i_narrow5 = ~(cellfun(@isempty,regexpi(stimcode,'no5gridnarrow'))); % including use of arbitrary ND
% 
% i_horizthird = strcmpi('horizgridthird',stimcode);
% i_horizthirdn = strcmpi('horizgridthirdn',stimcode);
% i_horizs30w10 = strcmpi('horizgrids30w10',stimcode);
% i_horizs30w15 = strcmpi('horizgrids30w15',stimcode);
% i_horizs40w15 = strcmpi('horizgrids40w15',stimcode);
% i_horizs40w20 = strcmpi('horizgrids40w20',stimcode);
% i_horizs40w25 = strcmpi('horizgrids40w25',stimcode);
% i_horizs40w30 = strcmpi('horizgrids40w30',stimcode);
% i_horizs50w20 = strcmpi('horizgrids50w20',stimcode);
% i_horizs50w25 = strcmpi('horizgrids50w25',stimcode);
% i_horizs50w30 = strcmpi('horizgrids50w30',stimcode);
% i_horizs50w35 = strcmpi('horizgrids50w35',stimcode);
% i_horizs50w40 = strcmpi('horizgrids50w40',stimcode);
% i_horizs60w30 = strcmpi('horizgrids60w30',stimcode);
% i_horizs60w40 = strcmpi('horizgrids60w40',stimcode);
% i_horizs60w45 = strcmpi('horizgrids60w45',stimcode);
% 
% i_verts30w10 = strcmpi('vertgrids30w10',stimcode);
% i_vertthird = strcmpi('vertgridthird',stimcode);
% 
% i_movesinev1 = strcmpi('movesinev1',stimcode);
% i_movesinev2 = strcmpi('movesinev2',stimcode);
% i_movesinev3 = strcmpi('movesinev3',stimcode);
% i_movesinev4 = strcmpi('movesinev4',stimcode);
% i_movesinev4rperm = strcmpi('movesinev4_rperm',stimcode);
% i_cylsinewave_horiz = strcmpi('cyl_sinewave_intergray_horiz',stimcode);
% i_cylsinewave_vert = strcmpi('cyl_sinewave_intergray_vert',stimcode);
% i_cylsinewave2_horiz = strcmpi('cyl_sinewave_intergray_horiz2',stimcode);
% i_cylsinewave2_vert = strcmpi('cyl_sinewave_intergray_vert2',stimcode);
% i_cylsinewave_horiz_hc = strcmpi('cyl_sinewave_intergray_horiz_hc',stimcode);
% i_cylsinewave_vert_hc = strcmpi('cyl_sinewave_intergray_vert_hc',stimcode);
% i_cylsquarewave_horiz = strcmpi('cyl_squarewave_intergray_horiz',stimcode);
% i_cylsquarewave_vert = strcmpi('cyl_squarewave_intergray_vert',stimcode);
% i_cylsinewave_horiz2 = strcmpi('cyl_sinewave_intergray_horiz2',stimcode);
% i_cylsinewave_vert2 = strcmpi('cyl_sinewave_intergray_vert2',stimcode);
% i_cylsinewave_horiz3_cor = strcmpi('cyl_sinewave_intergray_horiz3_cor',stimcode);
% i_cylsinewave_vert2_cor = strcmpi('cyl_sinewave_intergray_vert2_cor',stimcode);
% i_fffsinverperiod = strcmpi('fffsinverperiod',stimcode);
% 
% i_no5barrandompos = strcmpi('no5barrandompos',stimcode);
% i_vertbarrandompos = strcmpi('vertbarrandompos',stimcode);
% i_no5barrandompos4s = strcmpi('no5barrandompos4s',stimcode);
% i_horizbarperm2s2ND = strcmpi('horizbarrandompos2perm_2ND',stimcode);
% i_horizbarperm4s2ND = strcmpi('horizbarrandompos4perm_2ND',stimcode);
% i_vertbarperm2sFD = strcmpi('vertbarrandompos2perm_FD',stimcode);
% 
% i_barrandompos2s_no5permcrct = strcmpi('barrandompos2s_no5permcrct',stimcode);
% i_barrandompos2s_horizrpermwide = strcmpi('barrandompos2s_horizrpermwide',stimcode);
% i_barrandompos2s_vertrpermwide = strcmpi('barrandompos2s_vertrpermwide',stimcode);
% i_barrandompos3s10d_graydark_horiz = strcmpi('barrand_3s_10s_graydark_horiz',stimcode);
% i_barrandompos3s10d_graydark_vert = strcmpi('barrand_3s_10s_graydark_vert',stimcode);
% 
% i_no5movebargridslow = strcmpi('no5movebargridslow',stimcode);
% i_vertmovebargridslow = strcmpi('vertmovebargridslow',stimcode);
% 
% i_movebarslow_2dir1 = strcmpi('movebarslow_2dir1',stimcode);
% i_movebarslow_2dir2 = strcmpi('movebarslow_2dir2',stimcode);
% i_movebarslow_2dir1rperm = strcmpi('movebarslow_2dir1rperm',stimcode);
% i_movebarslow_2dir2rperm = strcmpi('movebarslow_2dir2rperm',stimcode);
% i_movebarslow_2dir1rpermf = strcmpi('movebarslow_2dir1rpermf',stimcode);
% i_movebarslow_2dir2rpermf = strcmpi('movebarslow_2dir2rpermf',stimcode);
% i_movebarslow_widths1 = strcmpi('movebarslow_widths1',stimcode);
% i_movebarslow_widths2 = strcmpi('movebarslow_widths2',stimcode);
% i_movebarslow_widths1f = strcmpi('movebarslow_widths1f',stimcode);
% i_movebarslow_widths2f = strcmpi('movebarslow_widths2f',stimcode);
% i_movebarslow_widths1rpermfhoriz = strcmpi('movebarslow_widths1rpermfhoriz',stimcode);
% i_movebarslow_widths2rpermfhoriz = strcmpi('movebarslow_widths2rpermfhoriz',stimcode);
% i_movebarslow_widths1rpermfvert = strcmpi('movebarslow_widths1rpermfvert',stimcode);
% i_movebarslow_widths2rpermfvert = strcmpi('movebarslow_widths2rpermfvert',stimcode);
% i_movebarslow_widths1rpermf = strcmpi('movebarslow_widths1rpermf',stimcode);
% i_movebarslow_widths2rpermf = strcmpi('movebarslow_widths2rpermf',stimcode);
% 
% i_guass200 = strcmpi('gaussian200',stimcode);

% i_15bouts = stimbouts>=15;
% i_inverted = (inverted == 1);
% i_moving = (moving == 1);
% i_active = (quality ~= 0);

i_M2   = (layer == 2);
i_M1 = (layer == 1); 
i_M5 = (layer == 5);
% i_lamina = (layer == 0);
i_M10 = (layer == 10);
% i_Lo1 = (layer == 11);
% i_Lo3 = (layer == 13);
% i_M9 = (layer == 9);
% i_cellBody = (layer == -1);
% i_M8 = (layer == 8);

% i_21Dhh = strcmpi('21Dhh',genotype);
% i_21Dhh_all = ~(cellfun(@isempty,regexpi(genotype,'21Dhh')));
% i_21Dhh_DG = strcmpi('21Dhh(DG)',genotype);
% i_21Dhh_mg = strcmpi('21Dhh-mg',genotype);
% i_21Dhh_LA = strcmpi('21Dhh-test-LA',genotype);
% i_21Dhh_all_tests = ~(cellfun(@isempty,regexpi(genotype,'21Dhh-test'))); 
% i_21Dhh_MS = ~(cellfun(@isempty,regexpi(genotype,'MS'))); 
% i_21Dhh_test = strcmpi('21Dhh-test',genotype); 
% i_21Dhh_test2 = strcmpi('21Dhh-test2',genotype); 
% i_21Dhh_test3 = strcmpi('21Dhh-test3',genotype); 
% i_21Dhh_test4 = strcmpi('21Dhh-test4',genotype); 
% i_21Dhh_test5 = strcmpi('21Dhh-test5',genotype); 
% i_21Dhh_test6 = strcmpi('21Dhh-test6',genotype); 
% i_21Dhh_PCT = strcmpi('21Dhh+PCT',genotype); 
% i_21Dhh_PCT_CGP = strcmpi('21Dhh+PCT+CGP',genotype); 
% i_21Dhh_rh1_GARi_GBRi = strcmpi('21Dhh_Rh1+GARi+GBRi',genotype);
% i_21Dhh_GARi_GBRi = strcmpi('21Dhh+GARi+GBRi',genotype);
% i_21Dhh_rh1_DCR2_GARi_GBRi = strcmpi('Z.3_Rh1_DCR2_GARi_GBRi',genotype);
% i_21Dhh_DCR2_GARi_GBRi = strcmpi('Z.3_DCR2_GARi_GBRi',genotype);
% i_21Dhh_DCR2 = strcmpi('Z.3_DCR2',genotype);
% i_AChR_antag = ~(cellfun(@isempty,regexpi(genotype,'Z.3_MEC')));
% i_Z3 = strcmpi('Z.3',genotype);
% i_Z3_PCT_CGP = strcmpi('Z.3+PCT+CGP',genotype);
% i_Z3_rh1 = strcmpi('Z.3+Rh1',genotype);
% i_L1202a = strcmpi('L1202a',genotype);
% i_L1CHL = strcmpi('L1CHL',genotype);
% i_L15A_L1202a = strcmpi('"L1(5A),L1(202a)"',genotype);
% 
% i_21Dhh_sh = strcmpi('21Dhh+sh',genotype); % with shits expression and heating
% % i_21Dhh_sh2 = strcmpi('21Dhh+sh2',genotype); % with shits expression and
% % heating + inline - never done
% i_21Dhh_sh_ctrl = strcmpi('21Dhh-sh',genotype); % without shits expression, with heating
% i_21Dhh_sh_ctrl2 = strcmpi('21Dhh-sh2',genotype); % with shits expression, without heating
% % i_21Dhh_sh_ctrl3 = strcmpi('21Dhh-sh3',genotype); % without shits
% % expression, without heating + inline - this was never done
% i_21Dhh_L1CHL_sh = strcmpi('21Dhh+L1CHL+sh',genotype); % with shits expression and heating
% i_21Dhh_L1CHL_sh_ctrl = strcmpi('21Dhh+L1CHL-sh',genotype); % with shits expression, without heating
% i_21Dhh_L1201a_sh = strcmpi('21Dhh+L1(201a)+sh',genotype); % only heating prior to imaging
% i_21Dhh_L1201a_sh2 = strcmpi('21Dhh+L1(201a)+sh2',genotype); % heating + inline perfusion heating
% i_21Dhh_L1201a_sh_ctrl = strcmpi('21Dhh+L1(201a)-sh',genotype); % with shits expression, without heating
% 
% % RNAi lines
% i_GluCli = strcmpi('GluCli',genotype); % with shits expression, without heating
% i_GluCli_mg = strcmpi('GluCli-Mg',genotype); % with shits expression, without heating
% i_GBRi_2nd = strcmpi('GBRi_2nd',genotype); % with shits expression, without heating
% i_GBRi_2nd_mg = strcmpi('GBRi_2nd-Mg',genotype); % with shits expression, without heating
% i_PTX_4 = strcmpi('PTX.4-pb',genotype); % with shits expression, without heating
% i_PTX_4_Z3 = strcmpi('PTX.4_Z.3',genotype); % with shits expression, without heating
% i_GARi = strcmpi('GARi-pb',genotype); % with shits expression, without heating
% 
% % AChR RNAi lines
% i_Z3_101806 = strcmpi('Z.3+101806',genotype);
% i_Z3_48267=strcmpi('Z.3+48267',genotype);
% i_Z3_33824 = strcmpi('Z.3+33824',genotype);
% i_Z3_DCR2_33825 = strcmpi('Z.3+DCR2+33825',genotype);
% i_Z3_DCR2_33824 = strcmpi('Z.3+DCR2+33824',genotype);
% i_Z3_AChRi = strcmpi('Z.3+AChRi',genotype);
% i_Z3_48267_101806 = strcmpi('Z.3+48267+101806',genotype);
% 
% % L4 silencing
% i_Z3_987Q_Qshi = strcmpi('Z.3+987Q+Qshi',genotype);
% i_Z3_987Q = strcmpi('Z.3=987Q',genotype);

% i_490 = (wavelength == 490);
% i_575 = (wavelength == 575);
% i_562 = (wavelength == 562);

% new two-photon
% i_Z3n = strcmpi('Z.3n',genotype);
% i_Z3_101806n = strcmpi('Z.3+101806n',genotype);
% i_Z3_mec_tub = strcmpi('Z.3+mec+tub', genotype);
% i_Z3_987Q_Qshi_contam = strcmpi('Z.3+987Q+Qshi_contaminated',genotype);
% i_Z3new = strcmpi('Z.3new',genotype);
% i_Z3contam = strcmpi('Z.3contaminated',genotype);
% i_Z3_101806new = strcmpi('Z.3+101806new',genotype);
% i_Z3_48267n = strcmpi('Z.3+48267n',genotype);
% i_Z3_101806new2 = strcmpi('Z.3+101806new2',genotype);
% i_Z3_48267new = strcmpi('Z.3+48267new',genotype);
% 
% i_cylsinewave_horiz3_cor_ss5c100cf05 = strcmpi('cyl_sinewave_intergray_horiz3_cor_ss5c100cf05',stimcode);
% i_cylsinewave_vert2_cor_ss5c100cf05 = strcmpi('cyl_sinewave_intergray_vert2_cor_ss5c100cf05',stimcode);
% i_Z3130328 = strcmpi('Z.3130328',genotype);


% gcamp6

% L2
% i_21Dhh_gcamp6m = strcmpi('21Dhh-gcamp6m',genotype);
% i_21Dhh_gcamp6s = strcmpi('21Dhh-gcamp6s',genotype);

% i_21Dhh_gcamp6f = strcmpi('21Dhh-gcamp6f',genotype);
% i_21Dhh_gcamp6f_mec_tub = strcmpi('21Dhh-gcamp6f+mec+tub', genotype);
% i_21Dhh_gcamp6f_101806 = strcmpi('21Dhh-gcamp6f+101806', genotype);
% i_21Dhh_gcamp6f_48267 = strcmpi('21Dhh-gcamp6f+48267', genotype);
% i_21Dhh_gcamp6f_27251 = strcmpi('21Dhh-gcamp6f+27251', genotype);
% i_21Dhh_gcamp6f_48267_101806_27251 = ...
%     strcmpi('21Dhh-gcamp6f+248267+101806+27251', genotype);
% i_21Dhh_gcamp6f_picro = strcmpi('21Dhh-gcamp6f+picro', genotype);

% i_R53C02AD_R29G11DBD_gcamp6f = strcmpi('R53C02AD-R29G11DBD-gcamp6f', genotype);
% i_R82F12AD_R65H08DBD_gcamp6f = strcmpi('R82F12AD-R65H08DBD-gcamp6f', genotype);

% % ASAP2f and Red Ca Indicators
% i_21Dhh_A147SdA148_jRGECO1b = strcmpi('21Dhh-A147SdA148-jRGECO1b',genotype);
% i_21Dhh_A147SdA148_jRCaMP1b = strcmpi('21Dhh-A147SdA148-jRCaMP1b',genotype);
% i_21Dhh_A147SdA148_jRCaMP1a = strcmpi('21Dhh-A147SdA148-jRCaMP1a',genotype);
% 
% % FlpStop genotypes
% % paraFlpstop  
% i_paraFlpStop_male_R53C02AD_R29G11DBD_gcamp6f = strcmpi('paraFlpStop-male-R53C02AD-R29G11DBD-gcamp6f', genotype);
% i_paraFlpStop_male_12days_R53C02AD_R29G11DBD_gcamp6f = strcmpi('paraFlpStop-male-12days-R53C02AD-R29G11DBD-gcamp6f', genotype);
% i_paraFlpStopControl_NoUASFlp_male_R53C02AD_R29G11DBD_gcamp6f = strcmpi('paraFlpStopControl-NoUASFlp-male-R53C02AD-R29G11DBD-gcamp6f', genotype);
% i_paraFlpStopControl_female_R53C02AD_R29G11DBD_gcamp6f = strcmpi('paraFlpStopControl-female-R53C02AD-R29G11DBD-gcamp6f', genotype);
% 
% i_R53G02AD_R29G11DBD_ASAP2f_para_no_flp_control = strcmpi('R53G02AD-R29G11DBD-ASAP2f-para-no-flp-control', genotype);
% i_R53G02AD_R29G11DBD_ASAP2f_para_flp = strcmpi('R53G02AD-R29G11DBD-ASAP2f-para-flp', genotype);

% % cac Tm3 genotypes 
% i_GMR13E12_gcamp6f_tubGal80ts_cac217_Flp = strcmpi('GMR13E12-gcamp6f_tubGal80ts_cac217_Flp',genotype);
% i_GMR13E12_gcamp6f_tubGal80ts_cac217_maleControl = strcmpi('GMR13E12-gcamp6f_tubGal80ts_cac217_maleControl',genotype); % male conrol 
% i_GMR13E12_gcamp6f_tubGal80ts_cac217_Flp_femaleControl = strcmpi('GMR13E12-gcamp6f_tubGal80ts_cac217_Flp_femaleControl',genotype);
% 
% %L1
% % i_c202a_gcamp6f = strcmpi('c202a-gcamp6f',genotype);
% % i_c202a_gcamp6f_mec_tub = strcmpi('c202a-gcamp6f+mec+tub',genotype);
% 
% 
% % new stimuli
% % i_barrandompos3s10d_dgraydark_horiz_onsc = strcmpi('barrand_3s_10s_dgraydark_horiz_onsc',stimcode);
% % i_barrandompos3s10d_dgraylight_horiz_onsc = strcmpi('barrand_3s_10s_dgraylight_horiz_onsc',stimcode);
% % i_barrandompos3s10d_dgraylight_horiz_onsc_cor = strcmpi('barrand_3s_10s_dgraylight_horiz_onsc_cor',stimcode);
% % i_circle_morph_graydark = ~(cellfun(@isempty,regexpi(stimcode,'circle_morph_rperm_local_graydark_cell')));
% % i_circle_morph_graybright = ~(cellfun(@isempty,regexpi(stimcode,'circle_morph_rperm_local_graybright_cell')));
% i_fff2s = strcmpi('fff2s',(stimcode));
% i_fff2s_15deg_border = strcmpi('fff2s_15deg_border',(stimcode));
% i_fff05s_cor = strcmpi('fff0.5s_cor',(stimcode));
% i_fff01s = strcmpi('fff0.1s',(stimcode));
% i_fff0024s = strcmpi('fff0.024s',(stimcode));
% i_fff03s = strcmpi('fff0.3s',(stimcode));
% i_mb_rpermcslof = strcmpi('movebar_rpermcsl_of',(stimcode));
% i_mb_rpermcsl_60dps = strcmpi('moving_single_bar_grid_rperm_csl_60dps', ...
%     (stimcode));
% i_gauss25 = strcmpi('gaussian25',stimcode);
% i_fff025sL05sD = strcmpi('fff0.025sL0.5sD',stimcode);
% i_fff025sD05sL = strcmpi('fff0.025sD0.5sL',stimcode);
% i_fff025sDL05sG = strcmpi('fff0.025sDL0.5sG',stimcode);
% i_fff025sDL05sGM = strcmpi('fff0.025sDL0.5sGM',stimcode);
% i_fff025sDL05sGM0125c = strcmpi('fff0.025sDL0.5sGM0.125c',stimcode);
% i_fff025sDL05sGM025c = strcmpi('fff0.025sDL0.5sGM0.25c',stimcode);
% i_fff025sL1sD = strcmpi('fff0.025sL1sD',stimcode);
% i_fff025sD1sL = strcmpi('fff0.025sD1sL',stimcode);
% i_fff025sDL1sG = strcmpi('fff0.025sDL1sG',stimcode);
% i_fff025sDL15sG = strcmpi('fff0.025sDL1.5sG',stimcode); 
% i_fff025sDL15sGM = strcmpi('fff0.025sDL1.5sGM',stimcode);
% i_fff025sDL15sGM0125c = strcmpi('fff0.025sDL1.5sGM0.125c',stimcode);
% i_fff025sDL15sGM025c = strcmpi('fff0.025sDL1.5sGM0.25c',stimcode);
% i_fff008sL05sD = strcmpi('fff0.008sL0.5sD',stimcode);
% i_fff008sD05sL = strcmpi('fff0.008sD0.5sL',stimcode);
% i_fff008sL15sD = strcmpi('fff0.008sL1.5sD',stimcode);
% i_fff008sD15sL = strcmpi('fff0.008sD1.5sL',stimcode);
% i_fff008sDL05sGM = strcmpi('fff0.008sDL0.5sGM',stimcode);
% i_fffsinesweep1 = strcmpi('fff_sine_sweep1',stimcode);
% i_fff008sDL15sGM = strcmpi('fff0.008sDL1.5sGM',stimcode);
% 
% % Stimuli with different excitation wavelengths 
% i_fff2s_1050nm = strcmpi('fff2s_1050nm',(stimcode));
% i_fff300ms_900nm = strcmpi('fff0.3s_900nm',(stimcode));
% i_fff300ms_920nm = strcmpi('fff0.3s_920nm',(stimcode));
% i_fff300ms_1000nm = strcmpi('fff0.3s_1000nm',(stimcode));
% i_fff300ms_1050nm = strcmpi('fff0.3s_1050nm',(stimcode));
% 
% % impulse stimuli
% i_fffdiffLengthFlashes = strcmpi('fff_diffLengthFlashes',stimcode);
% i_fffdiffLengthFlashes_5sgrey = strcmpi('fff_diffLengthFlashes_5sgrey',stimcode);
% 
% % white noise stimuli
% i_fffrand9_16ms = strcmpi('fff_rand9_16ms', stimcode); 
% i_fffrand9_50ms = strcmpi('fff_rand9_50ms', stimcode);
% 
% % wavelengths
% i_562 = strcmpi('562',wavelength);
% i_562_40_514_30 = strcmpi('562/40+e514/30', wavelength);
% i_447_60 = strcmpi('447/60', wavelength);
% i_448_20 = strcmpi('448/20', wavelength);
% i_447_60_nd1 = strcmpi('447/60+ND1', wavelength);
% i_447_60_nd05 = strcmpi('447/60+ND0.5', wavelength);
% i_447_60_900 = strcmpi('447/60+ex900',wavelength);
% i_447_60_920 = strcmpi('447/60+ex920',wavelength);
% i_447_60_940 = strcmpi('447/60+ex940',wavelength);
% i_447_60_noND = strcmpi('447/60-noND',wavelength);
% i_562_40_514_30_ff = strcmpi('ff-562/40+e514/30', wavelength);
% i_447_60_nd05_new150622 = strcmpi('447/60-ND0.5new150622', wavelength);
% i_447_60_noND_new150622 = strcmpi('447/60-noNDnew150622', wavelength);
% i_482_18_ND05_newDet = strcmpi('482/18-ND0.5_newDet',wavelength);
% 
% 
% % genotypes
% i_21Dhh_I67T = strcmpi('21Dhh-I67T',genotype);
% i_21Dhh_I67H = strcmpi('21Dhh-I67H',genotype);
% i_21Dhh_R154Q = strcmpi('21Dhh-R154Q',genotype);
% i_g105_I67T = strcmpi('G105-I67T',genotype);
% i_21Dhh_macmCitrine = strcmpi('21Dhh-macmCitrine',genotype);
% i_21Dhh_I67T_tdTOM = strcmpi('21Dhh-I67T-tdTOM', genotype);
% i_R48A08AD_R66A01DBD_I67T = strcmpi('R48A08AD-R66A01DBD-I67T', genotype);
% i_R48A08AD_R66A01DBD_gcamp6f = strcmpi('R48A08AD-R66A01DBD-GCaMP6f', genotype);
% i_21Dhh_R415Q = strcmpi('21Dhh-R415Q',genotype);
% i_21Dhh_A147SdA148 = strcmpi('21Dhh-A147SdA148',genotype);
% i_21Dhh_A147SdA148R415Q = strcmpi('21Dhh-A147SdA148R415Q',genotype);
% i_21Dhh_ArcLight = strcmpi('21Dhh-ArcLight',genotype);
% i_R48A08AD_R66A01DBD_A147SdA148 = strcmpi('R48A08AD-R66A01DBD-A147SdA148',genotype);
% i_otd_A147SdA148_GMRGal80 = strcmpi('otd-A147SdA148_GMR-Gal80',genotype);
% i_GMR19F01_A147SdA148 = strcmpi('GMR19F01-A147SdA148',genotype);
% i_GMR74G01_A147SdA148 = strcmpi('GMR74G01-A147SdA148',genotype);
% i_GMR13E12_A147SdA148 = strcmpi('GMR13E12-A147SdA148',genotype);
% i_21Dhh_A147SdA148_VK00005 = strcmpi('21Dhh-A147SdA148_VK00005',genotype);
% i_21Dhh_R415Q_VK00005 = strcmpi('21Dhh-R415Q_VK00005',genotype);
% i_GMR19F01_gcamp6f = strcmpi('GMR19F01-gcamp6f',genotype);
% i_GMR74G01_gcamp6f = strcmpi('GMR74G01-gcamp6f',genotype);
% i_GMR74G01_gcamp6m = strcmpi('GMR74G01-gcamp6m',genotype);
% i_GMR13E12_gcamp6f = strcmpi('GMR13E12-gcamp6f',genotype);
% i_GMR13E12_gcamp6m = strcmpi('GMR13E12-gcamp6m',genotype);
% i_21Dhh_cyRFP = strcmpi('21Dhh-gcamp6f-cyRFP',genotype);
% i_21Dhh_mCD8cyRFP = strcmpi('21Dhh-gcamp6f-mCD8cyRFP',genotype);
% i_otd_gcamp6f_GMRGal80 = strcmpi('otd-gcamp6f_GMR-Gal80',genotype);
% i_Rh1_gcamp6f = strcmpi('Rh1-gcamp6f', genotype);
% i_Rh1_A147SdA148 = strcmpi('Rh1-A147SdA148', genotype);
% i_GMR74G01_I67T = strcmpi('GMR74G01-I67T', genotype);
% i_GMR13E12_I67T = strcmpi('GMR13E12-I67T', genotype);
% i_21Dhh_I67QQ397R = strcmpi('21Dhh-I67QQ397R',genotype);
% i_21Dhh_cpOPT2 = strcmpi('21Dhh-cpOPT2', genotype);
% 
% % Flp-Stop genotypes
% i_GMR19F01_gcamp6f_para216_Flp = strcmpi('GMR19F01-gcamp6f_para216_Flp',genotype);
% i_GMR19F01_gcamp6f_para216 = strcmpi('GMR19F01-gcamp6f_para216',genotype);
% i_GMR19F01_gcamp6f_cac217_Flp = strcmpi('GMR19F01-gcamp6f_cac217_Flp',genotype);
% i_GMR19F01_gcamp6f_cac217 = strcmpi('GMR19F01-gcamp6f_cac217',genotype);
% 
% % TTX
% i_R48A08AD_R66A01DBD_gcamp6f_TTX = strcmpi('R48A08AD-R66A01DBD-GCaMP6f_TTX',genotype);
% i_21Dhh_gcamp6f_TTX = strcmpi('21Dhh-gcamp6f_TTX',genotype);
% i_GMR19F01_gcamp6f_para216_Flp_TTX = strcmpi('GMR19F01-gcamp6f_para216_Flp_TTX',genotype);
% i_GMR19F01_gcamp6f_para216_TTX = strcmpi('GMR19F01-gcamp6f_para216_TTX',genotype);

% L2
i_21Dhh_asap2f = strcmpi('21Dhh-ASAP2f', genotype);

% ASAP2f and jRGECO1b
i_21Dhh_asap2f_jrgeco1b = strcmpi('21Dhh-ASAP2f-jRGECO1b', genotype);

% ort silencing genotypes
i_21Dhh_asap2f_ort_TNT = strcmpi('21Dhh-ASAP2f_ortC1-3lexA-lexAopTNT',...
    genotype);
i_21Dhh_asap2f_ort = strcmpi('21Dhh-ASAP2f_ortC1-3lexA',genotype);
i_21Dhh_asap2f_TNT = strcmpi('21Dhh-ASAP2f_lexAopTNT',genotype);

i_21Dhh_asap2f_ort_TNT_TTX = strcmpi(...
    '21Dhh-ASAP2f_ortC1-3lexA-lexAopTNT_TTX',genotype);

i_GMR16H03_asap2f_ort_shits = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_ort[C1-4]-Gal4_UAS-shi[ts]R(5)+37-1hr',...
    genotype);
i_GMR16H03_asap2f_shits = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_UAS-shi[ts]R(5)+37-1hr',genotype);
i_GMR16H03_asap2f_ort = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_ort[C1-4]-Gal4+37-1hr',genotype);
i_GMR16H03_asap2f_ort_shits_4AP = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_ort[C1-4]-Gal4_UAS-shi[ts]R(5)+37-1hr+5mM_4-AP',...
    genotype);

i_R48A08AD_R66A01DBD_asap2f_ort_TNT = strcmpi(...
    'R48A08AD-R66A01DBD-ASAP2f_ortC1-3lexA-lexAopTNT',genotype);
i_R48A08AD_R66A01DBD_asap2f_ort = strcmpi(...
    'R48A08AD-R66A01DBD-ASAP2f_ortC1-3lexA',genotype);
i_R48A08AD_R66A01DBD_asap2f_TNT = strcmpi(...
    'R48A08AD-R66A01DBD-ASAP2f_lexAopTNT',genotype);

% TTX controls
i_21Dhh_gcamp6f_noTTX = strcmpi('21Dhh_GCaMP6f_noTTX',genotype);
i_21Dhh_gcamp6f_TTX = strcmpi('21Dhh_GCaMP6f_TTX',genotype);

% L2 Ace2N-mNeon
i_21Dhh_ace2nmneon_jrgeco1b = strcmpi('21Dhh-Ace2NmNeon-jRGECO1b', ...
    genotype);


% T4/T5
i_R59E08AD_R42F06DBD_gcamp6f = strcmpi('R59E08AD-R42F06DBD-GCaMP6f',...
    genotype);
i_VT025965_gcamp6f = strcmpi('VT025965-GCaMP6f', genotype);
% para FlpStop
i_R59E08AD_R42F06DBD_gcamp6f_paraFlpStop_Flp = strcmpi(...
    'paraGDY216_GCaMP6f_R59E08AD-R42F06DBD_Flp', genotype);
i_R59E08AD_R42F06DBD_gcamp6f_paraFlpStop_Flp_het = strcmpi(...
    'paraGDY216/w_GCaMP6f_R59E08AD-R42F06DBD_Flp', genotype);
i_R59E08AD_R42F06DBD_gcamp6f_paraFlpStop = strcmpi(...
    'paraGDY216_GCaMP6f_R59E08AD-R42F06DBD', genotype);

% Rdl FlpStop
i_R59E08AD_R42F06DBD_gcamp6f_Rdl1_RdlFlpStop_Flp = strcmpi(...
    'R59E08AD-R42F06DBD_GCaMP6f_Flp_Rdl1_RdlGDY217', genotype);
i_R59E08AD_R42F06DBD_gcamp6f_Rdl1_RdlFlpStop = strcmpi(...
    'R59E08AD-R42F06DBD_GCaMP6f_Rdl1_RdlGDY217',genotype);
i_R59E08AD_R42F06DBD_gcamp6f_RdlFlpStop_Flp = strcmpi(...
    'R59E08AD-R42F06DBD_GCaMP6f_Flp_RdlGDY217', genotype);

% L2 ort rescue
i_21Dhh_asap2f_UASort_Df3RBSC809 = strcmpi(...
    '21Dhh-ASAP2f_UAS-ort_Df(3R)BSC809',genotype);
i_21Dhh_asap2f_ort1new_Df3RBSC809 = strcmpi(...
    '21Dhh-ASAP2f_ort1new_Df(3R)BSC809',genotype);
i_21Dhh_asap2f_UASort_ort1new = strcmpi(...
    '21Dhh-ASAP2f_UAS-ort_ort1new',genotype);
i_21Dhh_asap2f_UASort_ort1new_Df3RBSC809 = strcmpi(...
    '21Dhh-ASAP2f_UAS-ort_ort1new_Df(3R)BSC809',genotype);

% L2 ort rescue, L2 shi[ts]
i_21Dhh_asap2f_UASort_ort1new_Df3RBSC809_1h37 = strcmpi(...
    '21Dhh-ASAP2f_UAS-ort_ort1new_Df(3R)BSC809+1hr37', genotype);
i_21Dhh_asap2f_UASort_UASshits_ort1new_Df3RBSC809_1h37 = strcmpi(...
    '21Dhh-ASAP2f_UAS-ort_UAS-shits_ort1new_Df(3R)BSC809+1hr37', genotype);
i_21Dhh_asap2f_UASort_UASshits_ort1new_1h37 = strcmpi(...
    '21Dhh-ASAP2f_UAS-ort_UAS-shits_ort1new+1hr37', genotype);
i_21Dhh_asap2f_UASort_UASshits_Df3RBSC809_1h37 = strcmpi(...
    '21Dhh-ASAP2f_UAS-ort_UAS-shits_Df(3R)BSC809+1hr37', genotype);

% L2 lexA, lexAop ASAP2f; cell type-specific silencing
i_GMR16H03_asap2f = strcmpi('GMR16H03-lexA_lexAop-ASAP2f', genotype);
i_GMR16H03_asap2f_lexAopTNT = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_lexAop-TNT', genotype);
i_GMR16H03_asap2f_ort_TNT = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_ort[C1-4]-Gal4_UAS-TNT', genotype);
i_GMR16H03_asap2f_ort = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_ort[C1-4]-Gal4', genotype);
i_GMR16H03_asap2f_TNT = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_UAS-TNT', genotype);
i_GMR16H03_asap2f_R11D03AD_R19C10DBD_TNT = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_R11D03AD-R19C10DBD_UAS-TNT', genotype);
i_GMR16H03_asap2f_R11D03AD_R19C10DBD = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_R11D03AD-R19C10DBD', genotype);
i_GMR16H03_asap2f_R52H01AD_R17C11DBD_TNT = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_R52H01AD-R17C11DBD_UAS-TNT', genotype);
i_GMR16H03_asap2f_R52H01AD_R17C11DBD = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_R52H01AD-R17C11DBD', genotype);
i_GMR16H03_asap2f_R20C11AD_R48D11DBD_TNT = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_R20C11AD-R48D11DBD_UAS-TNT', genotype);
i_GMR16H03_asap2f_R20C11AD_R48D11DBD = strcmpi(...
    'GMR16H03-lexA_lexAop-ASAP2f_R20C11AD-R48D11DBD', genotype);


% L1
i_GMR37E04_asap2f = strcmpi('GMR37E04-ASAP2f', genotype);

% CDM
i_21Dhh_asap2f_CDM = strcmpi('21Dhh-ASAP2f-CDM20uM', genotype);
i_GMR37E04_asap2f_CDM = strcmpi('GMR37E04-ASAP2f-20uMCDM', genotype);

i_21Dhh_asap2f_ort_TNT_CDM = strcmpi(...
    '21Dhh-ASAP2f_ortC1-3lexA-lexAopTNT-CDM20uM',genotype);
i_21Dhh_asap2f_ort_CDM = strcmpi('21Dhh-ASAP2f_ortC1-3lexA-CDM20uM',...
    genotype);
i_21Dhh_asap2f_TNT_CDM = strcmpi('21Dhh-ASAP2f_lexAopTNT-CDM20uM',...
    genotype);

% narrow abdomen
i_nahar38_21Dhh_asap2f = strcmpi('na[har38]_21Dhh-ASAP2f', genotype);

% 21Dhh ASAP2f; L2-lexA silencing
i_21Dhh_asap2f_GMR16H03 = strcmpi('21Dhh-ASAP2f_GMR16H03lexA', genotype);
i_21Dhh_asap2f_GMR16H03_TNT = strcmpi(...
    '21Dhh-ASAP2f_GMR16H03lexA-lexAopTNT', genotype);