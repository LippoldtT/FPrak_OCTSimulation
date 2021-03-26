%% OCT - Simulation 
%This simple simulation will compute a simple OCT setup. With a given
%spectrum and and some parameters for the sample it will compute several
%important graphs etc. (SI-unities are choosen everywhere)
clear all
close all
%Options
center = false ; % if set to 'true' enables the centering  of the spectrum
window = false ; % if set to 'true' enables the window-function
zero_pad = false; % if set to 'true' enables the zero-padding
%It is advised to use the  zero_padding with the the centering. The
%zero_padding will then add symmetrical zeros to the centered and cut spectrum.
%% Spectra Import and Simple Graphs
%At the start we have to import the measured spectra. Here we define their
%origin in your system 
    file_red = 'E:\Uni\Master\WiSe_2020_2021\Projektpraktikum\simulation\Spektren\Spektrum_roteLED.txt';
    file_white = 'E:\Uni\Master\WiSe_2020_2021\Projektpraktikum\simulation\Spektren\Spektrum_weißeLED.txt';
    file_green = 'E:\Uni\Master\WiSe_2020_2021\Projektpraktikum\simulation\Spektren\Spektrum_grünerlaser.txt';
%before one can use these spectra, the commas are replaced with dots.
    comma2point_overwrite(file_red);
    comma2point_overwrite(file_white);
    comma2point_overwrite(file_green);
%for the next step a simple import function 'import_spectra' is used.
    [lambda_red, i_red] = import_spectra(file_red) ;
    [lambda_white, i_white] =import_spectra(file_white) ;  
    [lambda_green, i_green]= import_spectra(file_green) ; 
%conversion to SI-unities
    lambda_red = lambda_red*10^-9;
    lambda_green = lambda_green*10^-9;
    lambda_white = lambda_white*10^-9;
%conversion from [\lambda,I] to [k,I]
    k_red = 2*pi./(lambda_red);
    k_white = 2*pi./(lambda_white);
    k_green = 2*pi./(lambda_green);
% In the next step a simple graph for all spectra is created
    %red -intensity over  wavelength
     fig1=figure(1);
     plot(lambda_red,i_red);
     grid minor
     xlabel('wavelength $\lambda$ in m','Interpreter','latex');
     ylabel(' intensity counts','Interpreter','latex');
     axis([338*10^-9 823*10^-9 0 13500]);
    %green
     fig2=figure(2);
     plot(lambda_green,i_green);
     grid minor
     xlabel('wavelength $\lambda$ in m','Interpreter','latex');
     ylabel(' intensity counts','Interpreter','latex');
     axis([338*10^-9 823*10^-9 0 11500]);
    %white
     fig3=figure(3); 
     plot(lambda_white,i_white);
     grid minor
     xlabel('wavelength $\lambda$ in m','Interpreter','latex');
     ylabel(' intensity counts','Interpreter','latex');
     axis([338*10^-9 823*10^-9 0 13500]);
     % Normalisation and comparison picture
     fig4=figure(4);
     i_red = i_red./(max(i_red));
     i_green = i_green./(max(i_green));
     i_white = i_white./(max(i_white));
     plot(lambda_red, i_red,lambda_green, i_green./(max(i_green)),lambda_white, i_white./(max(i_white)))
     grid minor
     xlabel('wavelength $\lambda$ in nm','Interpreter','latex');
     ylabel('normalized intensity counts','Interpreter','latex');
     axis([338*10^-9 823*10^-9 0 1.05]);   
%% generating fields/spectrums for comparison, simple gauss
%This section generates a simple gauss-spectrum for comaparison. The
%characteristics for the gauss cruve can be adjusted. Furthermore some
%constants are defined.

%constants 
    n = 1;                                  %assumption for the refractive index
    c = 299792458;                          %speed of ligth
%spectrum characteristics
    lambda_gauss=(450:0.01:550)*10^-9;      %wavelength area and step-size
    lambda_0_gauss =  500*10^-9;            %central wavelength of the gauss-curve
    d_lambda_1_gauss = 501*10^-9;           %bandwith maximum as wavelentgh
    d_lambda_2_gauss = 499*10^-9;           %bandwith minimum as wavelentgh
%conversion to other variables
    nu_gauss = c./(n.*lambda_gauss);        %frequency
    w_gauss = 2*pi*nu_gauss;                %angular frequency
    k_gauss = 2*pi./(lambda_gauss);         %k-vector
    k_0_gauss=2*pi/lambda_0_gauss;          %k_0-vector
    delta_k_gauss = 2*pi/(d_lambda_2_gauss)-2*pi/(d_lambda_1_gauss);                            %bandwith as k-vector
    i_gauss = 1/(delta_k_gauss.*sqrt(pi)) * exp(-((k_gauss-k_0_gauss)./delta_k_gauss).^2);      %electric field amplitude, gauss 
    i_gauss = 1/max(i_gauss) * i_gauss ;                                                        %normalisation to 1
 %figure for the gauss
    fig5 = figure(5);
    plot(lambda_gauss,i_gauss)
    grid minor
    xlabel('wavelength $\lambda$ in m','Interpreter','latex');
    ylabel('normalized intensity ','Interpreter','latex');
%% interpolation of measured spectra on a equidistant k-vector 
% IN this section the measured spectra are fitted over their k_vectors to
% create spectra which  can be related to an equidistant k-vector.
%creating fits
    fit_red = fit(k_red,i_red,'cubicinterp');
    fit_green = fit(k_green,i_green,'cubicinterp');
    fit_white = fit(k_white,i_white,'cubicinterp');
%generating k-vector array and interpolating the function on the new
%k-vector
    %red
        dis_k_red = (max(k_red)-min(k_red))/1022;
        k_red_dis = min(k_red):dis_k_red:max(k_red);
        k_red_dis = k_red_dis.';
        i_red = fit_red(k_red_dis);
    %green
        dis_k_green = (max(k_green)-min(k_green))/1022;
        k_green_dis = min(k_green):dis_k_green:max(k_green);
        k_green_dis = k_green_dis.';
        i_green = fit_green(k_green_dis);
    %white
        dis_k_white = (max(k_white)-min(k_white))/1022;
        k_white_dis = min(k_white):dis_k_white:max(k_white);
        k_white_dis = k_white_dis.';
        i_white = fit_white(k_white_dis);             
%% depth-dependent electric field reflectrivity profil (TODO)
% normaly continuous, resulting from the continuously varying refractive index of biological tissues and other samples
% HERE: assume a series of real deltafunction reflections

    %distances of the reflector and sample    
        z_r = 0.20;                             %distance of the reference from  the beamsplitter
        z_1= 0.20+50*10^-6;                     %distance of the first 'mirror' (in the sample) from  the beamsplitter
        delta_z =50*10^-6;                      %distance to the next 'mirror'
        z_2=z_1+delta_z;                        %distance of the second 'mirror' (in the sample) from  the beamsplitter
        z_d = 50*10^-6;                         %distance of the spectrometer which would 'capture' the simulated signal
    %reflectivities
        r_r = 1;                                %refelctivity of the reference
        r_1 = 0.21;                             %reflectivity of z_1 (dirac shape of the sample lets us assume just one value at the specific distance)
        r_2 = 0.2;                              %reflectivity of z_2
%% centering of the function 
figure(42)
plot(k_white_dis,i_white)
if center
    mode_c = 'abs';
    percentile_c = 0.001 ;
    absolute_c = 0.002;
    [k_red_dis,i_red]=centering(k_red_dis,i_red,mode_c,absolute_c,percentile_c);
    [k_green_dis,i_green]=centering(k_green_dis,i_green,mode_c,absolute_c,percentile_c);
    [k_white_dis,i_white]=centering(k_white_dis,i_white,mode_c,absolute_c,percentile_c);
end
%creating plots for the 
figure(43)
plot(k_white_dis,i_white)
%% window function
% to minimize unwanted shapes in the spectrum you can multiplicate a
% normalized window function

if window
    sigma = 0.5;                          % sigma parameter for gaussian window
    windowtype = 'triangle';                  % specification of the desired window function,
    %choose between 'tophat','welch','von_hahn','hamming','triangle','gaussian'
    i_red = window_f(i_red,windowtype);
    i_green = window_f(i_green,windowtype);
    i_white = window_f(i_white,windowtype);
end
%% zero padding
%Simple zero padding function is used,where the spectrum is centred
%and for either procentual or absolut values zeros are added. This function
%will stop adding zeros at 0. Further info commented inside function.

if zero_pad
    % Options
    mode_pad = 'abs'; % choose between absolute 'abs' or percentile 'percent' padding
    percentile_pad = 0.1 ; %in percentile of whole spectrum, only possible from 0-1, added on both sides
    absolute_pad = 500; % absolute number of added zeros  on both sides
    % Change of Spectra
    [k_red_dis,i_red]=zero_padding(k_red_dis,i_red,mode_pad,absolute_pad,percentile_pad);
    [k_green_dis,i_green]=zero_padding(k_green_dis,i_green,mode_pad,absolute_pad,percentile_pad);
    [k_white_dis,i_white]=zero_padding(k_white_dis,i_white,mode_pad,absolute_pad,percentile_pad);
end
%% real simulation of OCT
responsivity = 1;

E_refl_red = sqrt(i_red)./sqrt(2) .* r_r .* exp(i .* 2 .* k_red_dis .* z_r) .* exp(i .* k_red_dis .* z_d);


E_obj_1_red = sqrt(i_red)./sqrt(2) .* r_1 .* exp(i .* 2 .* k_red_dis .* z_1) .* exp(i .* k_red_dis .* z_d);
E_obj_2_red = sqrt(i_red)./sqrt(2) .* r_2 .* exp(i .* 2 .* k_red_dis .* z_2) .* exp(i .* k_red_dis .* z_d);
E_sample_red = E_obj_1_red + E_obj_2_red;
E_ges_red = E_refl_red + E_sample_red;

I_Ausgabe_red = responsivity / 2 .* (E_ges_red).*conj(E_ges_red);

%green
E_refl_green = sqrt(i_green)./sqrt(2) .* r_r .* exp(i .* 2 .* k_green_dis .* z_r) .* exp(i .* k_green_dis .* z_d);


E_obj_1_green = sqrt(i_green)./sqrt(2) .* r_1 .* exp(i .* 2 .* k_green_dis .* z_1) .* exp(i .* k_green_dis .* z_d);
E_obj_2_green = sqrt(i_green)./sqrt(2) .* r_2 .* exp(i .* 2 .* k_green_dis .* z_2) .* exp(i .* k_green_dis .* z_d);
E_sample_green = E_obj_1_green + E_obj_2_green;
E_ges_green = E_refl_green + E_sample_green;

I_Ausgabe_green = responsivity / 2 .* (E_ges_green).*conj(E_ges_green);


%white
E_refl_white = sqrt(i_white)./sqrt(2) .* r_r .* exp(i .* 2 .* k_white_dis .* z_r) .* exp(i .* k_white_dis .* z_d);


E_obj_1_white = sqrt(i_white)./sqrt(2) .* r_1 .* exp(i .* 2 .* k_white_dis .* z_1) .* exp(i .* k_white_dis .* z_d);
E_obj_2_white = sqrt(i_white)./sqrt(2) .* r_2 .* exp(i .* 2 .* k_white_dis .* z_2) .* exp(i .* k_white_dis .* z_d);
E_sample_white = E_obj_1_white + E_obj_2_white;
E_ges_white = E_refl_white + E_sample_white;

I_Ausgabe_white = responsivity / 2 .* (E_ges_white).*conj(E_ges_white);



% % testausgabe
% I_1 = E_obj_1_red .* conj(E_obj_1_red);
% I_2 = E_obj_2_red .* conj(E_obj_2_red);
% I_ges = E_sample_red .* conj(E_sample_red);
% hold on
% plot(k_red_dis, I_1);
% plot(k_red_dis, I_2);
% plot(k_red_dis, I_ges);
% plot(k_red_dis, I_Ausgabe);
% plot(k_red_dis, i_d_red);
% legend
% hold off
%% spectral interferogram
%calculation of the intensity at the detector
rho = 1;            %responsivity of the detector (units Amperes/Watt),

i_d_gauss   = rho/4.*   (i_gauss.*(r_r+ r_1)) + rho/2.*(i_gauss.*((sqrt(r_r*r_1)*(cos(2.*k_gauss.*(z_r-z_1))))+(sqrt(r_r*r_2)*(cos(2.*k_gauss.*(z_r-z_2)))))) + rho/4 *(i_gauss.*sqrt(r_1*r_2).*cos(2.*k_gauss.*(z_1-z_2))) ;
i_d_red     = rho/4.*(i_red.*(r_r+r_1)) + rho/2.*(i_red.*((sqrt(r_r*r_1)*(cos(2.*k_red_dis.*(z_r-z_1))))+(sqrt(r_r*r_2)*(cos(2.*k_red_dis.*(z_r-z_2)))))) + rho/4 *(i_red.*sqrt(r_1*r_2).*cos(2.*k_red_dis.*(z_1-z_2))) ;
i_d_green   = rho/4.*(i_green.*(r_r+r_1)) + rho/2.*(i_green.*((sqrt(r_r*r_1)*(cos(2.*k_green_dis.*(z_r-z_1))))+(sqrt(r_r*r_2)*(cos(2.*k_green_dis.*(z_r-z_2)))))) + rho/4 *(i_green.*sqrt(r_1*r_2).*cos(2.*k_green_dis.*(z_1-z_2))) ;
i_d_white   = rho/4.*(i_white.*(r_r+r_1)) + rho/2.*(i_white.*((sqrt(r_r*r_1)*(cos(2.*k_white_dis.*(z_r-z_1))))+(sqrt(r_r*r_2)*(cos(2.*k_white_dis.*(z_r-z_2)))))) + rho/4 *(i_white.*sqrt(r_1*r_2).*cos(2.*k_white_dis.*(z_1-z_2))) ;


% plots for the itensity over k
%red
fig6=figure(6);
    plot(k_red_dis,i_d_red);
    grid minor
    axis([min(k_red_dis) max(k_red_dis) min(i_d_red) max(i_d_red)]);
    xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
    ylabel('intensity at the detector','Interpreter','latex');
%green    
fig7=figure(7);
    plot(k_green_dis,i_d_green);
    grid minor
    axis([min(k_green_dis) max(k_green_dis) min(i_d_green) max(i_d_green)]);
    xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
    ylabel('intensity at the detector','Interpreter','latex');
%white  
fig8=figure(8);
    plot(k_white_dis,i_d_white);
    grid minor
    axis([min(k_white_dis) max(k_white_dis) min(i_d_white) max(i_d_white)]);
    xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
    ylabel('intensity at the detector','Interpreter','latex');
%gauss
fig9=figure(9);
    plot(k_gauss,i_d_gauss);
    grid minor
    axis([min(k_gauss) max(k_gauss) min(i_d_gauss) max(i_d_gauss)]);
    xlabel('wave-vector $k$ in $\frac{1}{m}$','Interpreter','latex');
    ylabel('intensity at the detector','Interpreter','latex');    
%% inverse fourier-transformation
%simple inverse fourier transformation and the plot
%additional iFFTshift of the function to match the zero frequency components
IFFT_red=ifftshift(abs(ifft(I_Ausgabe_red)));
IFFT_green=ifftshift(abs(ifft(I_Ausgabe_green)));
IFFT_white=ifftshift(abs(ifft(I_Ausgabe_white)));
%IFFT_gauss=2*pi*abs(ifft(i_d_gauss));

% define the z-axis 
    z_red = (0:length(IFFT_red)-1)/(length(IFFT_red).*dis_k_red)*2*pi/2 ;
    z_green = (0:length(IFFT_green)-1)/(length(IFFT_green).*dis_k_green)*2*pi/2 ;
    z_white = (0:length(IFFT_white)-1)/(length(IFFT_white).*dis_k_white)*2*pi/2 ;
%z_gauss = (0:length(IFFT_gauss)-1)/(length(IFFT_gauss).*(k_gauss(2,1)-k_gaus(1,1)))*2*pi ;

% define the shifted z-axis

%red spectrum
    z_red_shift = zeros(length(z_red),1);   
    z_red_shift(1:round(length(z_red)/2)-1) = fliplr(- z_red(2:1:round(length(z_red)/2)));
    z_red_shift(round(length(z_red)/2)) = z_red(1);
    z_red_shift(round(length(z_red)/2)+1:1:length(z_red)) = z_red(2:1:round(length(z_red)/2));
    
%green spectrum
    z_green_shift = zeros(length(z_green),1);   
    z_green_shift(1:round(length(z_green)/2)-1) = fliplr(- z_green(2:1:round(length(z_green)/2)));
    z_green_shift(round(length(z_green)/2)) = z_green(1);
    z_green_shift(round(length(z_green)/2)+1:1:length(z_green)) = z_green(2:1:round(length(z_green)/2));
    
%white spectrum
    z_white_shift = zeros(length(z_white),1);   
    z_white_shift(1:round(length(z_white)/2)-1) = fliplr(- z_white(2:1:round(length(z_white)/2)));
    z_white_shift(round(length(z_white)/2)) = z_white(1);
    z_white_shift(round(length(z_white)/2)+1:1:length(z_white)) = z_white(2:1:round(length(z_white)/2));


%plot the inverse function
%red
fig10=figure(10);
    plot(z_red_shift,abs(IFFT_red));
    grid minor
    axis([min(z_red_shift) max(z_red_shift) min(abs(IFFT_red)) max(abs(IFFT_red))]);
    xlabel('z-axis  in $m$','Interpreter','latex');
    ylabel('','Interpreter','latex');
%green
fig11=figure(11);
    plot(z_green_shift,abs(IFFT_green));
    grid minor
    axis([min(z_green_shift) max(z_green_shift) min(abs(IFFT_green)) max(abs(IFFT_green))]);
    xlabel('z-axis  in $m$','Interpreter','latex');
    ylabel('','Interpreter','latex');
%white
fig12=figure(12);
    plot(z_white_shift,abs(IFFT_white));
    grid minor
    axis([min(z_white_shift) max(z_white_shift) min(abs(IFFT_white)) max(abs(IFFT_white))]);
    xlabel('z-axis  in $m$','Interpreter','latex');
    ylabel('','Interpreter','latex');
%gauss
% fig13=figure(13);
%     plot(z_gauss,abs(ifft(IFFT_gauss)));
%     grid minor
%     axis([min(z_gauss) max(z_gauss) min(IFFT_gauss) max(IFFT_gauss)]);
%     xlabel('z-axis  in $m$','Interpreter','latex');
%     ylabel('','Interpreter','latex');    
 
