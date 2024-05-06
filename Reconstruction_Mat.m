%% 

clear
clc


fringe_num = 4;
subband_num = 3;

Image = load('matlab.mat'); image_data = Image.img;

%%

figure; plot(image_data(:,10,10));

%%

clc

[pixels, Aline, Bline] = size(image_data);
selected_bline = 10;

oddIndices = 1:2:Bline;
evenIndices = 2:2:Bline;

pixel_start_points = [400, 2400, 3100];
pixel_end_points = [900-1, 3100-1, 3600-1];

images_subband = cell(1, subband_num);
for i = 1:subband_num
    images_subband{i} = image_data(pixel_start_points(i):pixel_end_points(i), :, :);
end

Images_Subband = cellfun(@(x) fft(x(:,:,:)), images_subband, 'UniformOutput', false);

%%

figure; imagesc_auto(log10(abs(Images_Subband{1}(:,:,selected_bline)).^2));clim([4.5 11]);
title(['Raw Image, 650 nm: frame ',num2str(selected_bline)]);

figure; imagesc_auto(log10(abs(Images_Subband{2}(:,:,selected_bline)).^2));clim([4.5 11]);
title(['Raw Image, 505 nm: frame ',num2str(selected_bline)]);

figure; imagesc_auto(log10(abs(Images_Subband{3}(:,:,selected_bline)).^2));clim([4.5 11]);
title(['Raw Image, 488 nm: frame ',num2str(selected_bline)]);

%% DC subtraction

clc

images_avg = cellfun(@(x) squeeze(mean(x(:,:,:), [2 3])), images_subband, 'UniformOutput', false);
interference = cellfun(@(x, y) x(:,:,:)-permute(repmat(y, [1,Bline,Aline]), [1 3 2]), ...
    images_subband, images_avg, 'UniformOutput', false);
interference_fft = cellfun(@(x) fft(x(:,:,:)), interference, 'UniformOutput', false);

%%

figure; imagesc_auto(log10(abs(interference_fft{1}(:,:,selected_bline)).^2));clim([4.5 7.5]);
title(['After DC Subtraction, 650 nm: frame ',num2str(selected_bline)]);

figure; imagesc_auto(log10(abs(interference_fft{2}(:,:,selected_bline)).^2));clim([4.3 8]);
title(['After DC Subtraction, 505 nm: frame ',num2str(selected_bline)]);

figure; imagesc_auto(log10(abs(interference_fft{3}(:,:,selected_bline)).^2));clim([4.5 7.5]);
title(['After DC Subtraction, 488 nm: frame ',num2str(selected_bline)]);

%% Resampling

clc

deltaphi = cell(1, subband_num);
wavenumber = cell(1, subband_num);

for i = 1:subband_num
    file_name_deltaphi = sprintf('deltaphi%d.mat', i);
    deltaphi{i} = load(file_name_deltaphi).extrp_deltaphi;
    file_name_wavenumber = sprintf('wavenumber%d.mat', i);
    wavenumber{i} = load(file_name_wavenumber).wavenumber;
end


%% 

figure
for i = 1:subband_num
    plot(deltaphi{i}, 'LineWidth', 1);
    hold on
end
legend('Red', 'Green', 'Cyan');

figure
for i = 1:subband_num
    plot(wavenumber{i}, 'LineWidth', 1);
    hold on
end
legend('Red', 'Green', 'Cyan');

%% Resampling 

resam_image = cellfun(@(phi,x,k) interp1(phi, x(:,:,:), k, "spline"), ...
    deltaphi, interference, wavenumber, 'UniformOutput', false);
Resam_Image = cellfun(@(x) fft(x(:,:,:)), resam_image, 'UniformOutput', false);

%% Resampling with FFTInterp

clc

resam_image_red = zeros(size(interference{1}, 1), size(interference{1}, ...
    2), size(interference{1}, 3));
resam_image_green = zeros(size(interference{2}, 1), size(interference{2}, ...
    2), size(interference{1}, 3));
resam_image_cyan = zeros(size(interference{3}, 1), size(interference{3}, ...
    2), size(interference{3}, 3));

pixshift_red = (wavenumber{1}-deltaphi{1})/(max(wavenumber{1}) ...
    -min(wavenumber{1}))*(length(wavenumber{1})-1);
pixshift_green = (wavenumber{2}-deltaphi{2})/(max(wavenumber{2}) ...
    -min(wavenumber{2}))*(length(wavenumber{2})-1);
pixshift_cyan = (wavenumber{3}-deltaphi{3})/(max(wavenumber{3}) ...
    -min(wavenumber{3}))*(length(wavenumber{3})-1);
%%

pixshift_red = pixshift_red';
pixshift_green = pixshift_green';
pixshift_cyan = pixshift_cyan';

interpMat_red = createInterpMatFromShifts(pixshift_red);
interpMat_green = createInterpMatFromShifts(pixshift_green);
interpMat_cyan = createInterpMatFromShifts(pixshift_cyan);

%%
clc

for i = 1 : size(interference{1}, 3)
    for j = 1 : size(interference{1}, 2)

        resam_image_red(:, j, i) = FFTInterp(deltaphi{1}, ...
            interference{1}(:, j, i), wavenumber{1}, interpMat_red);
    end
end

for i = 1 : size(interference{2}, 3)
    for j = 1 : size(interference{2}, 2)

        resam_image_green(:, j, i) = FFTInterp(deltaphi{2}, ...
            interference{2}(:, j, i), wavenumber{2}, interpMat_green);
    end
end

for i = 1 : size(interference{3}, 3)
    for j = 1 : size(interference{3}, 2)

        resam_image_cyan(:, j, i) = FFTInterp(deltaphi{3}, ...
            interference{3}(:, j, i), wavenumber{3}, interpMat_cyan);
    end
end

%%
clc

resam_red_fft = fft(resam_image_red(:,:,:));
resam_green_fft = fft(resam_image_green(:,:,:));
resam_cyan_fft = fft(resam_image_cyan(:,:,:));

%%
clc

figure; imagesc_auto(log10(abs(resam_red_fft(:,:,selected_bline)).^2));clim([4.5 7.5]);
title(['After Resampling, Red: frame ',num2str(selected_bline)]);latexall(0)

figure; imagesc_auto(log10(abs(resam_green_fft(:,:,selected_bline)).^2));clim([4.3 8]);
title(['After Resampling, Green: frame ',num2str(selected_bline)]);latexall(0)

figure; imagesc_auto(log10(abs(resam_cyan_fft(:,:,selected_bline)).^2));clim([4.5 7.5]);
title(['After Resampling, Cyan: frame ',num2str(selected_bline)]);latexall(0)

%%

figure; imagesc_auto(log10(abs(Resam_Image{1}(:,:,selected_bline)).^2));clim([4.5 7.5]);
title(['After Resampling, 650 nm: frame ',num2str(selected_bline)]);

figure; imagesc_auto(log10(abs(Resam_Image{2}(:,:,selected_bline)).^2));clim([4.3 8]);
title(['After Resampling, 505 nm: frame ',num2str(selected_bline)]);

figure; imagesc_auto(log10(abs(Resam_Image{3}(:,:,selected_bline)).^2));clim([4.5 7.5]);
title(['After Resampling, 488 nm: frame ',num2str(selected_bline)]);

%% Dispersion compensation: Finding coefficients

format long

start_point_red = 70; end_point_red = 105;
start_point_green = 90; end_point_green = 125;
start_point_cyan = 45; end_point_cyan = 65;

objectiveFunction_650 = @(x) sharpness(x, Resam_Image{1}(start_point:end_point,:,selected_bline));
objectiveFunction_505 = @(x) sharpness(x, Resam_Image{2}(start_point:end_point,:,selected_bline));
objectiveFunction_488 = @(x) sharpness(x, Resam_Image{3}(start_point:end_point,:,selected_bline));

initialGuess = [0.000000000 0.00000000];  %3rd and 2nd order
iterationData = cell(3, 1);

iterationData{1} = struct('sharpness', [], 'coeff_results', []);
options = optimset('Display', 'iter', 'PlotFcns', @optimplotfval, ...
    'OutputFcn', @(x, optimValues, state) saveIterationResults(x, optimValues, state, 1));
coeff_result_650 = fminsearch(objectiveFunction_650, initialGuess, options);

iterationData{2} = struct('sharpness', [], 'coeff_results', []);
options = optimset('Display', 'iter', 'PlotFcns', @optimplotfval, ...
    'OutputFcn', @(x, optimValues, state) saveIterationResults(x, optimValues, state, 2));
coeff_result_505 = fminsearch(objectiveFunction_505, initialGuess, options);

iterationData{3} = struct('sharpness', [], 'coeff_results', []);
options = optimset('Display', 'iter', 'PlotFcns', @optimplotfval, ...
    'OutputFcn', @(x, optimValues, state) saveIterationResults(x, optimValues, state, 3));
coeff_result_488 = fminsearch(objectiveFunction_488, initialGuess, options);

%% Dispersion compensation

format long

l = length(start_point:end_point);
x_fit = cellfun(@(x) linspace(-l/2, l/2-1, size(x,1)), images_subband, 'UniformOutput', false);
PhaseCoeffs = {[coeff_result_650 0 0], [coeff_result_505 0 0], [coeff_result_488 0 0]};
nonlin_phase = cellfun(@(p,x) polyval(p, x), PhaseCoeffs, x_fit, 'UniformOutput', false);

compen_image = cell(1, subband_num);

for i = 1 : subband_num
    for j = 1 : Bline

        compen_image{i}(:,:,j) = ifft(Resam_Image{i}(:,:,j)).*exp(-1i.*nonlin_phase{i}');
    end
end

Compen_Image = cellfun(@(x) fft(x(:,:,:)), compen_image, 'UniformOutput', false);

%%

clc

save('compen_image', "compen_image", '-v7.3');

%%

figure; imagesc_auto(log10(abs(Compen_Image{1}(:,:, selected_bline).^2)));
clim([4.5 7.5]); title(['After Dispersion Compensation , 650 nm: frame ', num2str(selected_bline)]);

figure; imagesc_auto(log10(abs(Compen_Image{2}(:,:, selected_bline).^2)));
clim([4.3 8]); title(['After Dispersion Compensation , 505 nm: frame ', num2str(selected_bline)]);

figure; imagesc_auto(log10(abs(Compen_Image{3}(:,:, selected_bline).^2)));
clim([4.5 7.5]); title(['After Dispersion Compensation , 488 nm: frame ', num2str(selected_bline)]);

%%

figure;
subplot(2,1,1);imagesc_auto(log10(abs(Resam_Image{1}(:,:,selected_bline)).^2));clim([4.5 7.5]);
title(['After Resampling, 650 nm: frame ',num2str(selected_bline)]);
subplot(2,1,2); imagesc_auto(log10(abs(Compen_Image{1}(:,:,selected_bline)).^2));clim([4 6.5]);
title(['After Dispersion Compensation, 650 nm: frame ',num2str(selected_bline)]);

figure;
subplot(2,1,1);imagesc_auto(log10(abs(Resam_Image{2}(:,:,selected_bline)).^2));clim([4.5 7.5]);
title(['After Resampling, 505 nm: frame ',num2str(selected_bline)]);
subplot(2,1,2); imagesc_auto(log10(abs(Compen_Image{2}(:,:,selected_bline)).^2));clim([4 6.5]);
title(['After Dispersion Compensation, 505 nm: frame ',num2str(selected_bline)]);

figure;
subplot(2,1,1);imagesc_auto(log10(abs(Resam_Image{3}(:,:,selected_bline)).^2));clim([4.5 7.5]);
title(['After Resampling, 488 nm: frame ',num2str(selected_bline)]);
subplot(2,1,2); imagesc_auto(log10(abs(Compen_Image{3}(:,:,selected_bline)).^2));clim([4 6.5]);
title(['After Dispersion Compensation, 488 nm: frame ',num2str(selected_bline)]);

%% Checking dispersion compensation for subbands

clc

Resam_Image505 = cell(1, size(Resam_Image{2},1)/50);
Compen_Image505 = cell(1, size(Compen_Image{2},1)/50);
coeff_subband505 = cell(1, size(Resam_Image{2},1)/50);
objFunc_subband505 = cell(1, size(Resam_Image{2},1)/50);

for i = 1 : size(Resam_Image505, 2)

    Resam_Image505{i} = Resam_Image{2}(50*i-49:50*i, :, :);
end

for i = 1 : size(Resam_Image505, 2)

    objFunc_subband505{i} = @(x) sharpness(x, Resam_Image505{i}(:,:,selected_bline));
    coeff_subband505{i} = fminsearch(objFunc_subband505{i}, initialGuess, options);
end

x_fit_subband = -size(Resam_Image505{1}, 1)/2 : size(Resam_Image505{1}, 1)/2-1;

for i = 1 : size(Resam_Image505, 2)

    PhaseCoeffs_subband = [coeff_subband505{i} 0 0];
    nonlinPhase_subband = polyval(PhaseCoeffs_subband, x_fit_subband);

    for j = 1 : Bline

        Compen_Image505{i}(:,:,j) = fft(ifft(Resam_Image505{i}(:,:,j)).*exp(-1i.*nonlinPhase_subband'));
    end
end


% for i = 1 : size(Resam_Image505, 2)
% 
%     figure;
%     subplot(2,1,1);imagesc_auto(log10(abs(Resam_Image505{i}(:,:,selected_bline)).^2));clim([4.5 7.5]);
%     title(['After Resampling, 505 nm: frame ',num2str(selected_bline), ', Subband ', num2str(i)]);
%     subplot(2,1,2); imagesc_auto(log10(abs(Compen_Image505{i}(:,:,selected_bline)).^2));clim([4 6.5]);
%     title(['After Dispersion Compensation, 505 nm: frame ',num2str(selected_bline), ', Subband ', num2str(i)]);
% end

%%

for i = 1 : size(Resam_Image505, 2)

    figure;
    imagesc_auto(log10(abs(Compen_Image505{i}(:,:,selected_bline)).^2));clim([4 6.5]);
    title(['After Dispersion Compensation, 505 nm: frame ',num2str(selected_bline), ', Subband ', num2str(i)]);
end


%%

results_650 = struct('coefficients', [], 'sharpness_output', []);
results_505 = struct('coefficients', [], 'sharpness_output', []);
results_488 = struct('coefficients', [], 'sharpness_output', []);

iter = 0;

while true
    iter = iter + 1;

    % For image 650
    [coeff_result_650, sharpness650, exitflag_650] = fminsearch(objectiveFunction_650, initialGuess, options);
    results_650(iter).coefficients = coeff_result_650;
    results_650(iter).sharpness_output = -sharpness650;

    % For image 505
    [coeff_result_505, sharpness505, exitflag_505] = fminsearch(objectiveFunction_505, initialGuess, options);
    results_505(iter).coefficients = coeff_result_505;
    results_505(iter).sharpness_output = -sharpness505;

    % For image 488
    [coeff_result_488, sharpness488, exitflag_488] = fminsearch(objectiveFunction_488, initialGuess, options);
    results_488(iter).coefficients = coeff_result_488;
    results_488(iter).sharpness_output = -sharpness488;

    % Check exit conditions for all images
    if exitflag_650 <= 0 && exitflag_505 <= 0 && exitflag_488 <= 0
        break;  % Exit the loop if fminsearch did not converge
    end
end
