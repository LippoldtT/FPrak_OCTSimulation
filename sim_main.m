%% OCT - Simulation

%created by Georg Schwartz and Tom Lippoldt

%This simple simulation will compute a simple OCT setup. With a given
%spectrum and and some parameters for the sample it will compute several
%important graphs etc. (SI-unities are choosen everywhere)
clear all
close all
%Options
    center = false ; %if set to 'true' enables the centering  of the spectrum
    window = false ; %if set to 'true' enables the window-function
    zero_pad = false; %if set to 'true' enables the zero-padding
%It is advised to use the  zero_padding with the centering option. The
%zero_padding will then add symmetrical zeros to the centered and cut spectrum.
%% Spectra Import and Simple Graphs
%At the start we have to import the measured spectra. Here we define their
%origin in the operating system 
    file_red = 'E:\Uni\Master\WiSe_2020_2021\Projektpraktikum\simulation\Spektren\Spektrum_roteLED.txt';
    file_white = 'E:\Uni\Master\WiSe_2020_2021\Projektpraktikum\simulation\Spektren\Spektrum_weißeLED.txt';
    file_green = 'E:\Uni\Master\WiSe_2020_2021\Projektpraktikum\simulation\Spektren\Spektrum_grünerlaser.txt';
%before one can use these spectra, the commas are replaced with dots.
    comma2point_overwrite(file_red);
    comma2point_overwrite(file_white);
    comma2point_overwrite(file_green);
%simple import function 'import_spectra' is used.
    [lambda_red, i_red] = import_spectra(file_red) ;
    [lambda_white, i_white] =import_spectra(file_white) ;  
    [lambda_green, i_green]= import_spectra(file_green) ; 
%conversion to SI-units
    lambda_red = lambda_red*10^-9;
    lambda_green = lambda_green*10^-9;
    lambda_white = lambda_white*10^-9;
%conversion from wavelength to wavenumber
    k_red = 2*pi./(lambda_red);
    k_white = 2*pi./(lambda_white);
    k_green = 2*pi./(lambda_green);
%exclusion of a plotting vector for i_red, i_white, i_green for final
%display at the end of the script
    i_red_disp_spectra = i_red;
    i_green_disp_spectra = i_green;
    i_white_disp_spectra = i_white;
%% generating fields/spectrums for comparison, simple gauss
%This section generates a simple gauss-spectrum for comaparison. The
%characteristics for the gauss cruve can be adjusted. Furthermore some
%constants are defined.
%constants 
    n = 1;                                                                                      %assumption for the refractive index
    c = 299792458;                                                                              %speed of ligth
%spectrum characteristics
    lambda_gauss=(450:0.01:550)*10^-9;                                                          %wavelength area and step-size
    lambda_0_gauss =  500*10^-9;                                                                %central wavelength of the gauss-curve
    d_lambda_1_gauss = 510*10^-9;                                                               %bandwith maximum as wavelentgh
    d_lambda_2_gauss = 490*10^-9;                                                               %bandwith minimum as wavelentgh
%conversion to other variables
    nu_gauss = c./(n.*lambda_gauss);                                                            %frequency
    w_gauss = 2*pi*nu_gauss;                                                                    %angular frequency
    k_gauss = transpose(2*pi./(lambda_gauss));                                                  %k-vector
    k_0_gauss=2*pi/lambda_0_gauss;                                                              %k_0-vector
    delta_k_gauss = 2*pi/(d_lambda_2_gauss)-2*pi/(d_lambda_1_gauss);                            %bandwith as k-vector
    i_gauss = 1/(delta_k_gauss.*sqrt(pi)) * exp(-((k_gauss-k_0_gauss)./delta_k_gauss).^2);      %electric field amplitude, gauss 
    i_gauss = 1/max(i_gauss) * i_gauss *10000;                                                  %normalisation to 10000 counts to match the other spectra
%exclusion of a plotting vector for i_gauss for final
%display at the end of the script
    i_gauss_disp_spectra = i_gauss;
%% interpolation of measured spectra on a equidistant k-vector 
%In this section the measured spectra are fitted over their k_vectors to
%create spectra which  can be related to an equidistant k-vector.
%creating fits
    fit_red = fit(k_red,i_red,'cubicinterp');                               %cubic fit of the red spectra over the varied k vector
    fit_green = fit(k_green,i_green,'cubicinterp');                         %cubic fit of the green spectra over the varied k vector
    fit_white = fit(k_white,i_white,'cubicinterp');                         %cubic fit of the white spectra over the varied k vector
    fit_gauss = fit(k_gauss,i_gauss,'cubicinterp');                         %cubic fit of the gaussian spectra over the varied k vector
%generating k-vector array and interpolating the function on the new
%uniform k-vector
    %red
        dis_k_red = (max(k_red)-min(k_red))/(length(k_red)-1);              %calculation of the new stepsize
        k_red_dis = min(k_red):dis_k_red:max(k_red);                        %creation of the new uniform k-vector array
        k_red_dis = k_red_dis.';                                            %transposing vector to match dimensions
        i_red = fit_red(k_red_dis);                                         %fit of the spectra to the new k-vector
    %green
        dis_k_green = (max(k_green)-min(k_green))/(length(k_green)-1);      %calculation of the new stepsize
        k_green_dis = min(k_green):dis_k_green:max(k_green);                %creation of the new uniform k-vector array
        k_green_dis = k_green_dis.';                                        %transposing vector to match dimensions
        i_green = fit_green(k_green_dis);                                   %fit of the spectra to the new k-vector
    %white
        dis_k_white = (max(k_white)-min(k_white))/(length(k_white)-1);      %calculation of the new stepsize
        k_white_dis = min(k_white):dis_k_white:max(k_white);                %creation of the new uniform k-vector array
        k_white_dis = k_white_dis.';                                        %transposing vector to match dimensions
        i_white = fit_white(k_white_dis);                                   %fit of the spectra to the new k-vector
    %gaussian
        dis_k_gauss = (max(k_gauss)-min(k_gauss))/(length(k_gauss)-1);      %calculation of the new stepsize
        k_gauss_dis = min(k_gauss):dis_k_gauss:max(k_gauss);                %creation of the new uniform k-vector array
        k_gauss_dis = k_gauss_dis.';                                        %transposing vector to match dimensions
        i_gauss = fit_gauss(k_gauss_dis);                                   %fit of the spectra to the new k-vector
%% Modifications:Centering of the function 
%To apply the following modifications correctly we first have to delete all
%unnecessary noise values. We are left with only the part of the spectra
%where the intensity counts are over a significant value. To determine this
%value one can use either a fixed number or a percentual factor.
if center
    mode_c = 'abs';                                                     %'abs' = absolute value, 'percent' = percentual factor
    absolute_c = 30;                                                        %absolute value where all intensity count values below are removed
    percentile_c = 0.001 ;                                                  %percentual factor to the maxium of the spectrum where all intensity counts below the resulting value are removed
    [k_red_dis,i_red]=centering(k_red_dis,i_red,mode_c,absolute_c,percentile_c);            
    [k_green_dis,i_green]=centering(k_green_dis,i_green,mode_c,absolute_c,percentile_c);
    [k_white_dis,i_white]=centering(k_white_dis,i_white,mode_c,absolute_c,percentile_c);
    [k_gauss_dis,i_gauss]=centering(k_gauss_dis,i_gauss,mode_c,absolute_c,percentile_c);
end
%% Modifications:Window function
%To minimize unwanted shapes in the spectrum you can multiplicate a
%normalized window function
if window
    sigma = 0.5;                                                            %sigma parameter for the gaussian window (can be neglected for all others)
    windowtype = 'tophat';                                                  %specification of the desired window function, chose between: 'tophat', 'welch', 'von_hahn', 'hamming', 'triangle', 'gaussian'
    
    i_red = window_f(i_red,windowtype);
    i_green = window_f(i_green,windowtype);
    i_white = window_f(i_white,windowtype);
    i_gauss = window_f(i_gauss,windowtype);
end
%% Modifications:Zero padding
%Simple zero padding function is used, where for either procentual or absolut values zeros are added. This function
%will stop adding zeros at 0.
if zero_pad
    mode_pad = 'abs';                                                       %choose between absolute 'abs' or percentile 'percent' zero padding
    percentile_pad = 0.1 ;                                                  %percentile factor of whole spectrum length, only possible from 0-1, added on both sides
    absolute_pad = 500;                                                     %absolute number of added zeros on both sides
    
    [k_red_dis,i_red]=zero_padding(k_red_dis,i_red,mode_pad,dis_k_red,absolute_pad,percentile_pad);
    [k_green_dis,i_green]=zero_padding(k_green_dis,i_green,mode_pad,dis_k_green,absolute_pad,percentile_pad);
    [k_white_dis,i_white]=zero_padding(k_white_dis,i_white,mode_pad,dis_k_white,absolute_pad,percentile_pad);
    [k_gauss_dis,i_gauss]=zero_padding(k_gauss_dis,i_gauss,mode_pad,dis_k_gauss,absolute_pad,percentile_pad);
end
%exclusion of a plotting vector for i_red, i_white, i_green for final
%display at the end of the script
    i_red_disp_modif = i_red;
    i_green_disp_modif = i_green;
    i_white_disp_modif = i_white;
    i_gauss_disp_modif = i_gauss;
%% depth-dependent electric field reflectrivity profile
% normaly continuous, resulting from the continuously varying refractive index of biological tissues and other samples
% HERE: assume a series of real deltafunction reflections
    %distances of the reflector and sample    
        z_r = 0.20;                                                         %distance of the reference from  the beamsplitter
        z_1= 0.20+100*10^-6;                                                 %distance of the first 'mirror' (in the sample) from  the beamsplitter
        delta_z =50*10^-6;                                                  %distance to the next 'mirror'
        z_2=z_1+delta_z;                                                    %distance of the second 'mirror' (in the sample) from  the beamsplitter
        z_d = 0.5;                                                     %distance of the spectrometer which would 'capture' the simulated signal
    %reflectivities
        r_r = 1;                                                            %refelctivity of the reference
        r_1 = 0.21;                                                         %reflectivity of z_1 (dirac shape of the sample lets us assume just one value at the specific distance)
        r_2 = 0.2;                                                          %reflectivity of z_2
    %detector
        rho = 1;                                                            %responsistivity of the detector
%% spectral interferrogram
%The incident light is emitted from the beamsplitter and hits the reference
%Reflector and the sample. We create one wave for every reflection, so
%three in total. All 3 Waves get added together. To obtain the Intesity we
%take the amount square with respect to our complex values.
%red spectrum
    %wave from the reference Reflector
        E_refl_red = sqrt(i_red)./sqrt(2) .* r_r .* exp(i .* 2 .* k_red_dis .* z_r) .* exp(i .* k_red_dis .* z_d);                                         %#ok<*IJCL>
    %wave from the sample
        E_obj_1_red = sqrt(i_red)./sqrt(2) .* r_1 .* exp(i .* 2 .* k_red_dis .* z_1) .* exp(i .* k_red_dis .* z_d);
        E_obj_2_red = sqrt(i_red)./sqrt(2) .* r_2 .* exp(i .* 2 .* k_red_dis .* z_2) .* exp(i .* k_red_dis .* z_d);
        E_sample_red = E_obj_1_red + E_obj_2_red;                           %addition of both sample waves
    %whole wave
        E_ges_red = E_refl_red + E_sample_red;                              %addition of reference and whole sample wave
    %calculation of the intesity
        I_Ausgabe_red = rho / 2 .* (E_ges_red).*conj(E_ges_red);            %amount square with respect to the complex values

%green spectrum
    %wave from the reference Reflector
        E_refl_green = sqrt(i_green)./sqrt(2) .* r_r .* exp(i .* 2 .* k_green_dis .* z_r) .* exp(i .* k_green_dis .* z_d);
    %wave from the sample
        E_obj_1_green = sqrt(i_green)./sqrt(2) .* r_1 .* exp(i .* 2 .* k_green_dis .* z_1) .* exp(i .* k_green_dis .* z_d);
        E_obj_2_green = sqrt(i_green)./sqrt(2) .* r_2 .* exp(i .* 2 .* k_green_dis .* z_2) .* exp(i .* k_green_dis .* z_d);
        E_sample_green = E_obj_1_green + E_obj_2_green;                     %addition of both sample waves  
    %whole wave
        E_ges_green = E_refl_green + E_sample_green;                        %addition of reference and whole sample wave
    %calculation of the intesity
        I_Ausgabe_green = rho / 2 .* (E_ges_green).*conj(E_ges_green);      %amount square with respect to the complex values

%white spectrum
    %wave from the reference Reflector
        E_refl_white = sqrt(i_white)./sqrt(2) .* r_r .* exp(i .* 2 .* k_white_dis .* z_r) .* exp(1i .* k_white_dis .* z_d);
    %wave from the sample
        E_obj_1_white = sqrt(i_white)./sqrt(2) .* r_1 .* exp(i .* 2 .* k_white_dis .* z_1) .* exp(i .* k_white_dis .* z_d);
        E_obj_2_white = sqrt(i_white)./sqrt(2) .* r_2 .* exp(i .* 2 .* k_white_dis .* z_2) .* exp(i .* k_white_dis .* z_d);
        E_sample_white = E_obj_1_white + E_obj_2_white;                     %addition of both sample waves
    %whole wave
        E_ges_white = E_refl_white + E_sample_white;                        %addition of reference and whole sample wave
    %calculation of the intesity
        I_Ausgabe_white = rho / 2 .* (E_ges_white).*conj(E_ges_white);      %amount square with respect to the complex values

%gaussian spectrum
    %wave from the reference Reflector
        E_refl_gauss = sqrt(i_gauss)./sqrt(2) .* r_r .* exp(i .* 2 .* k_gauss_dis .* z_r) .* exp(i .* k_gauss_dis .* z_d);
    %wave from the sample
        E_obj_1_gauss = sqrt(i_gauss)./sqrt(2) .* r_1 .* exp(i .* 2 .* k_gauss_dis .* z_1) .* exp(i .* k_gauss_dis .* z_d);
        E_obj_2_gauss = sqrt(i_gauss)./sqrt(2) .* r_2 .* exp(i .* 2 .* k_gauss_dis .* z_2) .* exp(i .* k_gauss_dis .* z_d);
        E_sample_gauss = E_obj_1_gauss + E_obj_2_gauss;                     %addition of both sample waves
    %whole wave
        E_ges_gauss = E_refl_gauss + E_sample_gauss;                        %addition of reference and whole sample wave
    %calculation of the intesity
        I_Ausgabe_gauss = rho / 2 .* (E_ges_gauss).*conj(E_ges_gauss);      %amount square with respect to the complex values
%% inverse fourier-transformation
%simple inverse fourier transformation with additional iFFTshift of the 
%function to match the zero frequency components
%iFFT with shift of the intesity vector
    IFFT_red=ifftshift(abs(ifft(I_Ausgabe_red)));
    IFFT_green=ifftshift(abs(ifft(I_Ausgabe_green)));
    IFFT_white=ifftshift(abs(ifft(I_Ausgabe_white)));
    IFFT_gauss=ifftshift(abs(ifft(I_Ausgabe_gauss)));

%definition of the z-axis to obtain the depth information via plotting
%added factor of 2pi to match the missing factor while using k instead of
%frequency
%added factor of 1/2 since OCT looks at roundtrip distance instead of
%single distance
    z_red = (0:length(IFFT_red)-1)/(length(IFFT_red).*dis_k_red)*2*pi/2 ;   
    z_green = (0:length(IFFT_green)-1)/(length(IFFT_green).*dis_k_green)*2*pi/2 ;
    z_white = (0:length(IFFT_white)-1)/(length(IFFT_white).*dis_k_white)*2*pi/2 ;
    z_gauss = (0:length(IFFT_gauss)-1)/(length(IFFT_gauss).*dis_k_gauss)*2*pi/2 ;

%zero frequency shift of the z-axis to match the shift of the intensity
%red spectrum
    z_red_shift = zeros(length(z_red),1);                                                               %allocating memory
    z_red_shift(1:round(length(z_red)/2)-1) = fliplr(- z_red(2:1:round(length(z_red)/2)));              %swap of the first half of the vector
    z_red_shift(round(length(z_red)/2)) = z_red(1);                                                     %definition of the zero position                     
    z_red_shift(round(length(z_red)/2)+1:1:length(z_red)) = z_red(2:1:round(length(z_red)/2));          %swap of the second half of the vector
    
%green spectrum
    z_green_shift = zeros(length(z_green),1);                                                           %allocating memory
    z_green_shift(1:round(length(z_green)/2)-1) = fliplr(- z_green(2:1:round(length(z_green)/2)));      %swap of the first half of the vector
    z_green_shift(round(length(z_green)/2)) = z_green(1);                                               %definition of the zero position
    z_green_shift(round(length(z_green)/2)+1:1:length(z_green)) = z_green(2:1:round(length(z_green)/2));%swap of the second half of the vector
    
%white spectrum
    z_white_shift = zeros(length(z_white),1);                                                           %allocating memory
    z_white_shift(1:round(length(z_white)/2)-1) = fliplr(- z_white(2:1:round(length(z_white)/2)));      %swap of the first half of the vector
    z_white_shift(round(length(z_white)/2)) = z_white(1);                                               %definition of the zero position
    z_white_shift(round(length(z_white)/2)+1:1:length(z_white)) = z_white(2:1:round(length(z_white)/2));%swap of the second half of the vector
    
%gaussian spectrum
    z_gauss_shift = zeros(length(z_gauss),1);                                                           %allocating memory
    z_gauss_shift(1:round(length(z_gauss)/2)-1) = fliplr(- z_gauss(2:1:round(length(z_gauss)/2)));      %swap of the first half of the vector
    z_gauss_shift(round(length(z_gauss)/2)) = z_gauss(1);                                               %definition of the zero position
    z_gauss_shift(round(length(z_gauss)/2)+1:1:length(z_gauss)) = z_gauss(2:1:round(length(z_gauss)/2));%swap of the second half of the vector
%% Display of results
%red spectrum
    fig_disp_1 = figure('NumberTitle', 'off', 'Name', 'Red LED');
%input spectra
    subplot(2,2,1);
    plot(lambda_red,i_red_disp_spectra);
        title('input spectrum', 'FontSize', 12);
        grid minor
        xlabel('wavelength $\lambda$ in m','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([338*10^-9 823*10^-9 0 13500]);
%modified spectrum (k-space, centering, windowfunctions, zeropadding)        
    subplot(2,2,2);
        plot(k_red_dis,i_red_disp_modif);
    hold off
        title('modified spectrum', 'FontSize', 12);
        grid minor
        xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([min(k_red_dis) max(k_red_dis) min(i_red_disp_spectra) max(i_red_disp_spectra)]);
%simulated signal at the detector
    subplot(2,2,3);
        plot(k_red_dis,I_Ausgabe_red);
        title('simulated signal', 'FontSize', 12);
        grid minor
        xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([min(k_red_dis) max(k_red_dis) min(I_Ausgabe_red) max(I_Ausgabe_red)]);
%OCT depth structure
    subplot(2,2,4);
        plot(z_red_shift,abs(IFFT_red));
        title('OCT depth structure', 'FontSize', 12);
        grid minor
        axis([min(z_red_shift) max(z_red_shift) min(abs(IFFT_red)) max(abs(IFFT_red))]);
        xlabel('z-axis  in $m$','Interpreter','latex');
        ylabel('intensity counts','Interpreter','latex');
    
%green spectrum     
    fig_disp_2 = figure('NumberTitle', 'off', 'Name', 'Green Laserdiode');
%input spectra
    subplot(2,2,1);
    plot(lambda_green,i_green_disp_spectra);
        title('input spectrum', 'FontSize', 12);
        grid minor
        xlabel('wavelength $\lambda$ in m','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([338*10^-9 823*10^-9 0 13500]);
%modified spectrum (k-space, centering, windowfunctions, zeropadding)        
    subplot(2,2,2);
        plot(k_green_dis,i_green_disp_modif);
        title('modified spectrum', 'FontSize', 12);
        grid minor
        xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([min(k_green_dis) max(k_green_dis) min(i_green_disp_spectra) max(i_green_disp_spectra)]);            
%simulated signal at the detector
    subplot(2,2,3);
        plot(k_green_dis,I_Ausgabe_green);
        title('simulated signal', 'FontSize', 12);
        grid minor
        xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([min(k_green_dis) max(k_green_dis) min(I_Ausgabe_green) max(I_Ausgabe_green)]);
%OCT depth structure
    subplot(2,2,4);
        plot(z_green_shift,abs(IFFT_green));
        title('OCT depth structure', 'FontSize', 12);
        grid minor
        axis([min(z_green_shift) max(z_green_shift) min(abs(IFFT_green)) max(abs(IFFT_green))]);
        xlabel('z-axis  in $m$','Interpreter','latex');
        ylabel('intensity counts','Interpreter','latex');


%white_spectrum     
    fig_disp_3 = figure('NumberTitle', 'off', 'Name', 'White broadband-LED');
%input spectra
    subplot(2,2,1);
    plot(lambda_white,i_white_disp_spectra);
        title('input spectrum', 'FontSize', 12);
        grid minor
        xlabel('wavelength $\lambda$ in m','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([338*10^-9 823*10^-9 0 13500]);
%modified spectrum (k-space, centering, windowfunctions, zeropadding)        
    subplot(2,2,2);
        plot(k_white_dis,i_white_disp_modif);
        title('modified spectrum', 'FontSize', 12);
        grid minor
        xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([min(k_white_dis) max(k_white_dis) min(i_white_disp_spectra) max(i_white_disp_spectra)]);            
%simulated signal at the detector
    subplot(2,2,3);
        plot(k_white_dis,I_Ausgabe_white);
        title('simulated signal', 'FontSize', 12);
        grid minor
        xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([min(k_white_dis) max(k_white_dis) min(I_Ausgabe_white) max(I_Ausgabe_white)]);
%OCT depth structure
    subplot(2,2,4);
        plot(z_white_shift,abs(IFFT_white));
        title('OCT depth structure', 'FontSize', 12);
        grid minor
        axis([min(z_white_shift) max(z_white_shift) min(abs(IFFT_white)) max(abs(IFFT_white))]);
        xlabel('z-axis  in $m$','Interpreter','latex');
        ylabel('intensity counts','Interpreter','latex');
        
%gaussian_spectrum     
    fig_disp_4 = figure('NumberTitle', 'off', 'Name', 'Simulated gaussian shape');
%input spectra
    subplot(2,2,1);
    plot(lambda_gauss,i_gauss_disp_spectra);
        title('input spectrum', 'FontSize', 12);
        grid minor
        xlabel('wavelength $\lambda$ in m','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([min(lambda_gauss) max(lambda_gauss) min(i_gauss_disp_spectra) max(i_gauss_disp_spectra)]);
%modified spectrum (k-space, centering, windowfunctions, zeropadding)        
    subplot(2,2,2);
        plot(k_gauss_dis,i_gauss_disp_modif);
        title('modified spectrum', 'FontSize', 12);
        grid minor
        xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([min(k_gauss_dis) max(k_gauss_dis) min(i_gauss_disp_spectra) max(i_gauss_disp_spectra)]);            
%simulated signal at the detector
    subplot(2,2,3);
        plot(k_gauss_dis,I_Ausgabe_gauss);
        title('simulated signal', 'FontSize', 12);
        grid minor
        xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
        ylabel(' intensity counts','Interpreter','latex');
        axis([min(k_gauss_dis) max(k_gauss_dis) min(I_Ausgabe_gauss) max(I_Ausgabe_gauss)]);
%OCT depth structure
    subplot(2,2,4);
        plot(z_gauss_shift,abs(IFFT_gauss));
        title('OCT depth structure', 'FontSize', 12);
        grid minor
        axis([min(z_gauss_shift)*0.5e-1 max(z_gauss_shift)*0.5e-1 min(abs(IFFT_gauss)) max(abs(IFFT_gauss))]);
        xlabel('z-axis  in $m$','Interpreter','latex');
        ylabel('intensity counts','Interpreter','latex');

