function [Y]=window_f(y,mode,sigma)
%This function applies one of the predetermined window functions onto the
%input spectrum.
length_Window = length(y);                                                  %length of fitting array
Window = zeros(length_Window,1);                                            %allocating memory
Sample = 1:length(y);                                                       %sampling array
switch mode
    case 'tophat'                                                           %Tophat
        Window(Sample) = 1;
    case 'welch'                                                            %Welch (Kosinus)
        Window(Sample) = cos((pi.*Sample)./(length_Window) -(pi/2));
    case 'von_hahn'                                                         %Von Hann
        Window(Sample) = 0.5 * (1-cos(2*pi.*Sample./length_Window));
    case 'hamming'                                                          %Hamming
        Window(Sample) = 0.54 - 0.46*cos(2*pi.*Sample./length_Window);
    case 'triangle'                                                         %Triangle
        Window(Sample) = 2/length_Window * (length_Window/2 - abs(Sample - length_Window/2));
    case 'gaussian'                                                         %Gaussian
        Window(Sample) = exp(-0.5 *( (Sample - length_Window/2)./(sigma*length_Window/2)).^2);
    otherwise
        error('Die eingegebene Fensterfunktion wurde nicht erkannt!')
end
Y = y.* Window;                                                             %multiplikation of spectrum and window function