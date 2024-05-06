% OCT Data Reconstruction
% Double Phase

%% Data Preparation
clear
clc

fringe_num = 3;

sam1 = load('sam1.mat');
sam1_data = sam1.img;

[m, n, p] = size(sam1_data);

sam1_data = reshape(sam1_data, m, n*p);

% sam1_data = sam1_data(:, 3:end);

fringes = cell(1, fringe_num);
refs = cell(1, fringe_num);

for index = 1 : fringe_num
    
    fringe_name = sprintf('fringe%d.mat', index+1);
    ref_name = sprintf('ref%d.mat', index+1);

    fringes{index} = load(fringe_name).img;
    fringes{index} = reshape(fringes{index}, [], n*p);
%     figure; imshow(sam1_data, []); title(['Fringe', num2str(index)]);

    refs{index} = load(ref_name).img;
    refs{index} = reshape(refs{index}, m, n*p);
    
end

%% Truncation

trunc_point = 1;

fringes_used = cellfun(@(x) x(:, trunc_point:end), fringes, 'UniformOutput', false);
refs_used = cellfun(@(x) x(:, trunc_point:end), refs, 'UniformOutput', false);
sam1_used = sam1_data(:, trunc_point:end);

%% PSF Check before DC Subtraction

fringes_fft = cellfun(@fft, fringes_used, 'UniformOutput', false);

% figure
% for index = 1 : fringe_num
% 
%     subplot(fringe_num, 1, index);
%     plot(10.*log(abs(fringes_fft{index}(:, 20))));
%     title(['PSF of Fringe ', num2str(index+1), ' before DC Subtraction']);
%     ylabel('Log Intensity');
%     grid on
% end


%% DC Background Subtraction

sam1_mean = mean(sam1_used, 2);
fringes_mean = cellfun(@(x)  mean(x, 2), fringes_used, 'UniformOutput', false);
refs_mean = cellfun(@(x)  mean(x, 2), refs_used, 'UniformOutput', false);

sub_data = cellfun(@(x, y) x - sam1_mean - y, fringes_used, refs_mean, ...
    'UniformOutput', false);
sub_data = cellfun(@(x) x - mean(x, 1), sub_data, 'UniformOutput', false);

sub_fft = cellfun(@(x) fft(x, length(x)), sub_data, 'UniformOutput', false);

% figure
% for index = 1 : fringe_num
% 
%     subplot(fringe_num, 1, index);
%     plot(sub_data{index}(:, :));
%     title(['DC subtracted Fringe ', num2str(index+1)]);
%     xlabel('Pixel #');
%     ylabel('Intensity');
%     grid on
% end


% figure
% for index = 1 : fringe_num
% 
%     subplot(fringe_num, 1, index);
%     plot(10.*log(abs(sub_fft{index}(:, 20))));
%     title(['PSF of Fringe ', num2str(index+1), ' after DC Subtraction']);
%     ylabel('Log Intensity');
%     grid on
% end

%% Hilbert Transform & Phase Calculation

hilbert_factor_resam = size(fringes_used{1}, 1)/2;

hilbert_data = cellfun(@(cellData) arrayfun(@(col) MyHilbert(cellData(:, col), ...
    hilbert_factor_resam), 1:size(cellData, 2), 'UniformOutput', false), sub_data, ...
    'UniformOutput', false);
hilbert_data_mat = cellfun(@cell2mat, hilbert_data, 'UniformOutput', false);
phi = cellfun(@(x) unwrap(angle(x)), hilbert_data_mat, 'UniformOutput', false);
phi_mean = cellfun(@(x) mean(x, 2), phi, 'UniformOutput', false);

% figure
% for index = 1 : fringe_num
% 
%     subplot(fringe_num, 1, index);
%     plot(phi{index}(:, :));
%     title(['Phase for Fringe ', num2str(index+1)]);
%     xlabel('Pixel #');
%     ylabel('Phase');
% end

figure
for index = 1 : fringe_num

    plot(phi_mean{index}, 'LineWidth', 1.2, 'DisplayName', sprintf('Fringe %d', index+1));
    title('Average Phase before Resampling');
    xlabel('Pixel #');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');

%% Phase difference & Wavenumber

clc

deltaPhi32 = phi_mean{2} - phi_mean{1};
deltaPhi43 = phi_mean{3} - phi_mean{2};

figure
plot(deltaPhi32, 'LineWidth', 2);
hold on 
plot(deltaPhi43, 'LineWidth', 2);
title('Phase Differences before Resampling');
xlabel('# Pixels');
ylabel('Phase');
legend('\Delta\Phi_{32}', '\Delta\Phi_{43}')
grid on

%%

start_point = 137; %252
end_point = 2560;
trunc_x = start_point:end_point;

reliable_deltaPhi32 = deltaPhi32(start_point:end_point);
% reliable_deltaPhi43 = deltaPhi43(start_point:end_point);

trunc_x_axis = -length(trunc_x)/2 : length(trunc_x)/2-1;

reliable_deltaphi_coeffs = polyfit(trunc_x_axis, reliable_deltaPhi32, 7);
% reliable_deltaphi_coeffs = polyfit(trunc_x_axis, reliable_deltaPhi43, 7);
extrp_x = 1 : size(sam1_mean, 1);
extrp_x_axis = interp1(trunc_x, trunc_x_axis, extrp_x, 'linear', 'extrap');
extrp_deltaphi = polyval(reliable_deltaphi_coeffs, extrp_x_axis);

%% Phase extrapolation check

clc

figure
plot(deltaPhi32, 'LineWidth', 1.5);
hold on
plot(extrp_deltaphi, 'LineWidth', 1.5);
title('\Delta\Phi_{Actual}  vs \Delta\Phi_{Extrapolated}');
xlabel('Pixel #');
ylabel('Phase');
grid on
legend('Actual', 'Extrapolated');

phi_diffs = deltaPhi32 - extrp_deltaphi';
% phi_diffs = deltaPhi43 - extrp_deltaphi';

figure
plot(phi_diffs(:), 'LineWidth', 2);
title('(\Delta\Phi_{Actual} - \Delta\Phi_{Extrapolated})_{Before Resampling}');
xlabel('Pixel #');
ylabel('Phase');
grid on

wavenumber = linspace(min(extrp_deltaphi), max(extrp_deltaphi), length(extrp_deltaphi));

figure
plot(wavenumber, 'LineWidth', 1);
hold on
plot(deltaPhi32, 'LineWidth', 1);
title('Wavenumber vs \Delta\Phi');
legend('Wavenumber', '\Delta\Phi_{32}');
xlabel('Pixel #');
grid on
%% Resampling

clc

resam_Data = cellfun(@(cellData) arrayfun(@(col) interp1(extrp_deltaphi, ...
    cellData(:, col), wavenumber, 'spline'), 1:size(cellData, 2), 'UniformOutput', ...
    false), sub_data, 'UniformOutput', false);
resam_Data_tr = cellfun(@(x) cellfun(@(y) y', x, 'UniformOutput', false), ...
    resam_Data, 'UniformOutput', false);
resam_data = cellfun(@cell2mat, resam_Data_tr, 'UniformOutput', false);
resam_fft = cellfun(@(x) fft(x), resam_data, 'UniformOutput', false);

%% Plotting resampled fringes

clc

% figure
% for index = 1 : fringe_num
% 
%     subplot(fringe_num, 1, index);
%     plot(resam_data{index});
%     title(['Resampled Fringe ', num2str(index+1)]);
%     xlabel('Wavenumber');
%     ylabel('Intensity');
% end

figure
for index = 1 : fringe_num

    plot(10.*log(abs(resam_fft{index}(:,20))), 'LineWidth', 0.8, ...
        'DisplayName', sprintf('Fringe %d', index+1));
    title('PSF of Resampled data in one A-line');
    xlabel('Depth');
    ylabel('Log Intensity');
    xlim([0 2560]);
    grid on
    hold on
end
legend('show');


%% Fringes check after resampling

% for index = 1 : fringe_num
%     
%     figure
%     subplot(2, 1, 1);
%     plot(fringes{index}(:, 100));
%     title(['Before resampling Fringe ', num2str(index+1)]);
%     subplot(2, 1, 2); 
%     plot(resam_data{index}(:, 100));
%     title(['After resampling Fringe ', num2str(index+1)]);
% end

%% Phase calculation for resampled fringes

clc

hilbert_factor_dc = size(fringes_used{1}, 1)/4;

hilbert_resam = cellfun(@(cellData) arrayfun(@(col) MyHilbert(cellData(:, col), ...
    hilbert_factor_dc), 1:size(cellData, 2), 'UniformOutput', false), resam_data, ...
    'UniformOutput', false);
hilbert_resam = cellfun(@cell2mat, hilbert_resam, 'UniformOutput', false);

hilbert_resam_fft = cellfun(@(cellData) arrayfun(@(col) MyHilbert(cellData(:, col), ...
    hilbert_factor_dc), 1:size(cellData, 2), 'UniformOutput', false), resam_fft, ...
    'UniformOutput', false);
hilbert_resam_fft = cellfun(@cell2mat, hilbert_resam_fft, 'UniformOutput', false);
env_resam = cellfun(@(x) abs(x), hilbert_resam_fft, 'UniformOutput', false);

phi_resam = cellfun(@(x) unwrap(angle(x)), hilbert_resam, 'UniformOutput', false);
phi_resam_mean = cellfun(@(x) mean(x, 2), phi_resam, 'UniformOutput', false);

%% Phase plots for resampled fringes

figure
for index = 1 : fringe_num

    plot(phi_resam_mean{index}, 'LineWidth', 1, 'DisplayName', ...
        sprintf('Resampled Fringe %d', index+1))
    hold on
    title('Average Phases after Resampling');
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');

 
figure
for index = 1 : fringe_num

    plot(phi_resam_mean{index}, 'LineWidth', 1, 'DisplayName', ...
        sprintf('Resampled Fringe %d', index+1))
    hold on
    plot(phi_mean{index}, ':', 'LineWidth', 1, 'DisplayName', ...
        sprintf('Fringe %d', index+1));
    title('Average Phases before vs after Resampling');
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');

%% Checking the linearity of the resampled phase difference

clc

resam_deltaPhi32 = phi_resam_mean{2} - phi_resam_mean{1};
resam_deltaPhi43 = phi_resam_mean{3} - phi_resam_mean{2};

figure
plot(resam_deltaPhi32, 'LineWidth', 1);
hold on
plot(resam_deltaPhi43, 'LineWidth', 1);
legend('resam deltaPhi_3_2', 'resam deltaPhi_4_3');
title('Average Phase Differences after Resampling');
xlabel('Wavenumber');
ylabel('Phase');
grid on

sp = 1; %140
ep = 2560; %2540

resam_deltaPhi32_x = -length(deltaPhi32(sp:ep))/2 : length(deltaPhi32(sp:ep))/2 -1;
resam_deltaPhi43_x = -length(deltaPhi43(sp:ep))/2 : length(deltaPhi32(sp:ep))/2 -1;
resam_phi32_coeffs = polyfit(resam_deltaPhi32_x, resam_deltaPhi32(sp:ep), 2);
resam_phi43_coeffs = polyfit(resam_deltaPhi43_x, resam_deltaPhi43(sp:ep), 2);

%% Finding starting points for resampled phases

clc

figure
for index = 1 : fringe_num

    plot(phi_resam_mean{index}(1:700), 'LineWidth', 1, 'DisplayName', sprintf('Fringe%d', index+1));
    title('Averaged Phase After Resampling/First portion');
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');

figure
for index = 1 : fringe_num

    plot(phi_resam_mean{index}(2000:2560), 'LineWidth', 1, 'DisplayName', sprintf('Fringe%d', index+1));
    title('Averaged Phase after Resampling/Last portion');
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');

%% Resampled start and end points

clc

resam_start_points = 169; %296
resam_end_points = 2560;
resam_x = resam_start_points:resam_end_points;

%% Fitting on reliable resampled phase

clc

reliable_resam_x = -length(resam_x)/2:length(resam_x)/2-1;

reliable_resam_phase = cellfun(@(data) data(resam_start_points:resam_end_points), ...
    phi_resam_mean, 'UniformOutput', false);
% resam_coeffs_check = cellfun(@(y) polyfit(reliable_resam_x, y, 2), ...
%     reliable_resam_phase, 'UniformOutput', false);

fitting_order = 8;

resam_poly_coeffs = cellfun(@(y) polyfit(reliable_resam_x, y, fitting_order), ...
    reliable_resam_phase, 'UniformOutput', false);
resam_linpoly_coeffs = cellfun(@(y) polyfit(reliable_resam_x, y, 1), ...
    reliable_resam_phase, 'UniformOutput', false);
%%

clc

reliable_phase_fit_check = cellfun(@(p) polyval(p, reliable_resam_x), ...
    resam_poly_coeffs, 'UniformOutput', false);

reliable_phase_fit_diff = cellfun(@(x, y) x-y', reliable_resam_phase, ...
    reliable_phase_fit_check, 'UniformOutput', false);

reliable_phase_linfit_check = cellfun(@(p) polyval(p, reliable_resam_x), ...
    resam_linpoly_coeffs, 'UniformOutput', false);

reliable_phase_fitdif_lin_vs_nonlin = cellfun(@(x, y) x-y, reliable_phase_fit_check, ...
    reliable_phase_linfit_check, 'UniformOutput', false);

figure
for index = 1 : fringe_num

    plot(reliable_phase_fit_diff{index}, 'LineWidth', 1, ...
        'DisplayName', sprintf('Fringe %d', index+1))
    title('(Actual Phase - Fitted Phase)_{truncated}');
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');

figure
for index = 1 : fringe_num

    plot(reliable_resam_phase{index},':', 'LineWidth', 1, 'DisplayName', ...
        sprintf('Resampled Fringe %d', index+1));
    title('(Average Resampled Phase vs Fitted)_{truncated}');
    hold on
    plot(reliable_phase_fit_check{index}, 'LineWidth', 1, 'DisplayName', ...
        sprintf('Fitted Fringe %d', index+1));
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');

figure
for index = 1 : fringe_num

    plot(reliable_phase_fit_check{index},':', 'LineWidth', 1, 'DisplayName', ...
        sprintf('Higher order Fitted Fringe %d', index+1));
    title('(High Order Fit vs Linear Fit)_{truncated}');
    hold on
    plot(reliable_phase_linfit_check{index}, 'LineWidth', 1, 'DisplayName', ...
        sprintf('Linearly Fitted Fringe %d', index+1));
    hold on
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
end
legend('show');

figure
for index = 1 : fringe_num

    plot(reliable_phase_fitdif_lin_vs_nonlin{index}, 'LineWidth', ...
        1, 'DisplayName', sprintf('Fringe %d', index+1))
    title('(High Order Fit - Linear Fit)_{truncated}');
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');
%%
clc

extrp_x = 1 : size(sam1_mean, 1);
extrp_resam_x = interp1(resam_x, reliable_resam_x, extrp_x, 'linear', 'extrap');
extrp_resam_phi = cellfun(@(p) polyval(p, extrp_resam_x), resam_poly_coeffs, ...
    'UniformOutput', false);

figure
for index = 1 : fringe_num

    plot(phi_resam_mean{index},':', 'LineWidth', 1, 'DisplayName', ...
        sprintf('Resampled Fringe %d', index+1));
    hold on
    title('Average Resampled Phase vs Extrapolated');
    xlabel('Wavenumber');
    ylabel('Phase');
    plot(extrp_resam_phi{index}, 'LineWidth', 1, 'DisplayName', ...
        sprintf('Extrapolated Fringe %d', index+1));
    grid on
    hold on
end
legend('show');


extrp_phase_fit_diff = cellfun(@(x, y) x-y', phi_resam_mean, ...
    extrp_resam_phi, 'UniformOutput', false);

figure
for index = 1 : fringe_num

    plot(extrp_phase_fit_diff{index}(resam_start_points:end), 'LineWidth', 1, 'DisplayName', ...
        sprintf('Fringe %d', index+1))
    title('(Actual Phase - Fitted Phase)_{extrapolated}');
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');


%% 

clc

nonlin_coeffs = cellfun(@(x) [x(1:end-2), zeros(1, 2)], resam_poly_coeffs, ...
    'UniformOutput', false);

nonlin_phases = cellfun(@(p) polyval(p, extrp_resam_x), nonlin_coeffs, ...
    'UniformOutput', false);

%%

clc


figure
for index = 1 : fringe_num
    
    plot(nonlin_phases{index}, 'LineWidth', 1, 'DisplayName', ...
        sprintf('Fringe%d', index+1))
    title('Nonlinear Phases_{extrapolated}');
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');

figure
for index = 1 : fringe_num

    plot(nonlin_phases{index}(resam_start_points:end), 'LineWidth', 1, ...
        'DisplayName', sprintf('Fringe%d', index+1))
    title('Nonlinear Phases_{truncated}');
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');

%%

nonlinear_phase = (nonlin_phases{1}+nonlin_phases{2})/2;
nonlinear_phase = nonlinear_phase';

figure
plot(nonlinear_phase, 'LineWidth', 1);
title('Final nonlinear phase');
xlabel('Wavenumber');
ylabel('Phase');
grid on

figure
plot(nonlinear_phase(resam_start_points:end), 'LineWidth', 1);
title('Final nonlinear phase');
xlabel('Wavenumber');
ylabel('Phase');
grid on
%%

clc

compen_hilb = cellfun(@(x) x.*exp(-1i.*nonlinear_phase), hilbert_resam, ...
    'UniformOutput', false);
compen_data = cellfun(@real, compen_hilb, 'UniformOutput', false);
compen_fft = cellfun(@fft, compen_data, 'UniformOutput', false);
compen_phase = cellfun(@(x) unwrap(angle(x)), compen_hilb, 'UniformOutput', false);
compen_phase_mean = cellfun(@(x) mean(x, 2), compen_phase, 'UniformOutput', false);
%%

clc

figure
for index = 1 : fringe_num

    plot(compen_phase_mean{index}, 'LineWidth', 1, ...
        'DisplayName', sprintf('Fringe %d', index+1))
    title('Average Phase after Compensation');
    xlabel('Wavenumber');
    ylabel('Phase');
    grid on
    hold on
end
legend('show');


figure
for index = 1 : fringe_num

    plot(10.*log(abs(compen_fft{index}(:, 20))), 'LineWidth', 1, ...
        'DisplayName', sprintf('Fringe %d', index+1));
    title('PSF for Dispersion Compensated Fringes');
    xlabel('Depth');
    ylabel('Log Intensity');
    xlim([0 2560]);
    grid on
    hold on
end
legend('show');

%% PSF check separately

figure
plot(10.*log(abs(compen_fft{1}(:, 20))), 'LineWidth', 1); title('D.C. Fringe 2');
figure
plot(10.*log(abs(compen_fft{2}(:, 20))), 'LineWidth', 1); title('D.C. Fringe 3');
figure
plot(10.*log(abs(compen_fft{3}(:, 20))), 'LineWidth', 1); title('D.C. Fringe 4');


%% FWHM Check for Dispersion Compensated Fringes

clc

% FWHM for PSF

upsample_factor = 4;

max_compen_PSF = cellfun(@(x) max(abs(x)), compen_fft, 'UniformOutput', false);

compensated_PSF = cellfun(@(x) abs(fft(x, upsample_factor*length(sam1_mean))), ...
    compen_data, 'UniformOutput', false);

[compen_pks, compen_locs, compen_w, compen_p] = cellfun(@(cellData, max_val) ...
    arrayfun(@(col) findpeaks(cellData(1:2560, col), "WidthReference", ...
    "halfheight", "MinPeakHeight", max_val(col)-(1e-3), "MinPeakDistance", ...
    800), 1:size(cellData, 2), 'UniformOutput', false), compensated_PSF, ...
    max_compen_PSF, 'UniformOutput', false);

compen_w = cellfun(@(x) cell2mat(x), compen_w, 'UniformOutput', false);

A_line_num = 400;

for index = 1 : fringe_num

    fprintf(['FHWM of Compensated PSF for A_line number %d in Fringe %d:' ...
        ' %2.4f\n'], A_line_num, index+1, compen_w{index}(A_line_num));
end

%%

% FWHM for Envelope PSF

env_compen = cellfun(@abs, compen_hilb, 'UniformOutput', false);
env_compen_fft = cellfun(@(x) fftshift(fft(x)), env_compen, 'UniformOutput', false);

max_env_compen = cellfun(@(x) max(abs(x)), env_compen_fft, ...
    'UniformOutput', false);

env_compen_PSF = cellfun(@(x) abs(fftshift(fft(x, upsample_factor*length( ...
    sam1_mean)))), env_compen, 'UniformOutput', false);

[env_compen_pks, env_compen_locs, env_compen_w, env_compen_p] = cellfun(@(cellData, ...
    max_val) arrayfun(@(col) findpeaks(cellData(:, col), "WidthReference", ...
    "halfheight", "MinPeakHeight", max_val(col)-(1e-3), "MinPeakDistance", ...
    800), 1:size(cellData, 2), 'UniformOutput', false), env_compen_PSF, ...
    max_env_compen, 'UniformOutput', false);

env_compen_w = cellfun(@(x) cell2mat(x), env_compen_w, 'UniformOutput', false);

A_line_num = 400;

for index = 1 : fringe_num

    fprintf(['\nFHWM of Compensated Envelope PSF for A_line number %d in ' ...
        'Fringe %d: %2.4f'], A_line_num, index+1, env_compen_w{index}(A_line_num));
end

fprintf('\n')
%% Phase check for compensated data

clc

compen_poly_ceoffs = cellfun(@(data) polyfit(extrp_resam_x, data, 2), ...
    compen_phase_mean, 'UniformOutput', false);

compen_x = -length(compen_phase_mean{1})/2 : length(compen_phase_mean{1})/2-1;
compen_lin_coeffs = cellfun(@(x) polyfit(compen_x, x, 1), ...
    compen_phase_mean, 'UniformOutput', false);
compen_lin_fit = cellfun(@(p) polyval(p, compen_x), ...
    compen_lin_coeffs, 'UniformOutput', false);
compen_lin_fit_diff = cellfun(@(x, y) x-y', compen_phase_mean, ...
    compen_lin_fit, 'UniformOutput', false);

figure
for index = 1 : fringe_num

    plot(compen_phase_mean{index}, 'LineWidth', 1, 'LineStyle',':', ...
        'DisplayName', sprintf('Compensated Fringe %d', index+1));
    hold on
    plot(compen_lin_fit{index}, 'LineWidth', 1, 'DisplayName', ...
        sprintf('Linear Fit for Compensated Fringe %d', index+1));
    xlabel('Wavenumber');
    ylabel('Phase');
    title('Compensated Phase vs Linear Fit on Compensated Phase')
    grid on
    hold on
end
legend('Show');


figure
for index = 1 : fringe_num

    plot(compen_lin_fit_diff{index}, 'LineWidth', 1, ...
        'DisplayName', sprintf('Fringe %d', index+1));
    xlabel('Wavenumber');
    ylabel('Phase');
    title('Compensated Phase - Linear Fit')
    grid on
    hold on
end
legend('Show');
%%

clc


figure
for index = 1 : fringe_num

    plot(10.*log(abs(env_compen_fft{index}(:, 20))), 'LineWidth', 1, ...
        'DisplayName', sprintf('Fringe %d', index+1))
    title('Envelope of Compensated Fringes');
    xlabel('Depth');
    ylabel('Log Intensity');
    grid on
    hold on
end
legend('show');

%% Compensated PSF vs Envelope PSF

clc

freq = -size(env_compen_fft{1}, 1)/2 : size(env_compen_fft{1}, 1)/2-1;

figure
for index = 1 : fringe_num

    plot(freq, 10.*log(abs(env_compen_fft{index}(:, 20))), 'DisplayName', ...
        sprintf('Fringe %d', index+1))
    title('Envelope of Compensated Fringes');
    xlabel('Depth');
    ylabel('Log Intensity');
    grid on
    hold on
end
legend('show');


loc_shift = {compen_locs{1, 1}{1}, compen_locs{1, 2}{1}, compen_locs{1, 3}{1}};

for index = 1 : fringe_num
    
    figure
    plot(freq, 10.*log(abs(compen_fft{index}(:, 20))), 'LineWidth', 1, ...
        'DisplayName', sprintf('PSF of Fringe %d', index+1));
    hold on
    plot(freq+loc_shift{index}/upsample_factor+freq(1), ...
        10.*log(0.5.*abs(env_compen_fft{index}(:, 20))), 'LineWidth', 1, ...
        'DisplayName', sprintf('Envelope PSF of Fringe %d', index+1))
    title('PSF vs Envelope PSF of Compensated Fringes');
    xlabel('Depth');
    ylabel('Log Intensity');
    legend('show');
    grid on
end

%% Phase Comparison: Resampled vs Compensated 

clc

figure
for index = 1 : fringe_num

    plot(phi_resam_mean{index}, ':', 'LineWidth', 1, ...
        'DisplayName', sprintf('Resampled Fringe %d', index+1));
    hold on
    plot(compen_phase_mean{index}, 'LineWidth', 1, ...
        'DisplayName', sprintf('Compensated Fringe %d', index+1))
    xlabel('Wavenumber');
    ylabel('Phase');
    title('Compensated vs Resampled Phase')
    grid on
    hold on
end
legend('Show');

