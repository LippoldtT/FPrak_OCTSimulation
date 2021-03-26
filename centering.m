function [X,Y]=centering(x,y,mode,abs,percentile)
%This simple centering function removes all near zero values from the given
%y-array and matches the corresponding x array. The threshold can be
%determined by a absolute value (mode 'abs') or by a percentile factor
%(mode 'percentile').
switch mode
    case 'abs'                                                              %selection of the absolute value mode
        counter = 0;                                                        %defining counting variables
        j=1;
        for k = 1:length(y)                                                 %checking loop to get the final length of the new array to predeterminde its size to save running time for long arrays
            if y(k)>abs
                counter = counter+1;
            end
        end
        Y = zeros(counter,1);                                               %allocating the new array for y
        X = zeros(counter,1);                                               %allocating the new array for x
        for k = 1:length(y)
            if y(k)>abs
                Y(j) = y(k);                                                %copiing all values of the input y array into the new array if the threshold is surpassed
                X(j) = x(k);                                                %copiing all values of the input x array aswell to keep them corresponding
                j=j+1;
            end
        end
        if mod(counter,2)==0                                                %making sure the array has uneven number of elements to prevent the future FFTshift functions from failing
            Y(counter+1) = Y(counter);                                      %if needed a additional character is created wich is equivalent with the last above threshold value
            X(counter+1) = X(counter);
        end  
        
    case 'percent'                                                          %selection of the percentile facotor mode
        counter = 0;                                                        %defining counting variables
        j=1;
        for k = 1:length(y)                                                 %checking loop to get the final length of the new array to predeterminde its size to save running time for long arrays
            if y(k)>(max(y)*percentile)
                counter = counter+1;
            end
        end
        Y = zeros(counter,1);                                               %allocating the new array for y
        X = zeros(counter,1);                                               %allocating the new array for x
        for k = 1:length(y)
            if y(k)>(max(y)*percentile)
                Y(j) = y(k);                                                %copiing all values of the input y array into the new array if the threshold is surpassed
                X(j) = x(k);                                                %copiing all values of the input x array aswell to keep them corresponding
                j=j+1;
            end
        end
        if mod(counter,2)==0                                                %making sure the array has uneven number of elements to prevent the future FFTshift functions from failing
            Y(counter+1) = Y(counter);                                      %if needed a additional character is created wich is equivalent with the last above threshold value
            X(counter+1) = X(counter);                                      
        end
end