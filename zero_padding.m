function [X,Y]=zero_padding(x,y,mode,delta_x,abs,percent)
%This function adds equal amounts of zero value entries to both sides of
%the input y function. To match the corresponding x vector the same amount
%of values are added, with the distance given in delta_x. The amount of
%zeros can be determined with a absolute value (mode 'abs') or as
%percentile factor of the total length of the input vector(mode 'percent').
switch mode
    case 'abs'                                                              %selection of the absolute value mode
        if x(1,1)-(delta_x*abs) < 0                                         %checking if the given number of zeros would result in a negative vector
            error('Die angegeben absolute Anzahl an Nullen, führt zu einer negativen x-Achse !')
        else
            X = zeros(length(y)+2*abs,1);                                   %allocating memory for the new x vector
            Y = zeros(length(x)+2*abs,1);                                   %allocating memory for the new y vector
            i=abs;
            while i < length(X)-abs                                         %adding of new zero value enties to the y vector in positive and negative direction
                Y(i,1)=y(i-abs+1,1);
                i=i+1;
            end
            i=2;
            X(1,1) = x(1,1)- (delta_x*abs);                                 %adding new entries to the x vector with stepsize delta_x in positive and negative direction
            while i < length(Y)+1
                X(i,1) = X(i-1,1) + delta_x;
                i=i+1;
            end
        end
        
    case 'percent'                                                          %selection of the percentile facotor mode
        if x(1,1)-(length(x)*delta_x*percent) < 0                           %checking if the given number of zeros would result in a negative vector
            error('Die angegeben prozentuale Anzahl an Nullen, führt zu einer negativen x-Achse !')
        else
            Y = zeros(length(x)+2*round(length(x)*percent),1);              %allocating memory for the new y vector
            X = zeros(length(x)+2*round(length(x)*percent),1);              %allocating memory for the new x vector
            i=round(length(x)*percent);
            while i < length(X)-round(length(x)*percent)                    %adding of new zero value enties to the y vector in positive and negative direction
                Y(i,1)=y(i-round(length(x)*percent)+1,1);
                i=i+1;
            end
            i=2;
            X(1,1) = x(1,1)- (length(x)*delta_x*percent);                   %adding new entries to the x vector with stepsize delta_x in positive and negative direction
            while i < length(Y)+1
                X(i,1) = X(i-1,1) + delta_x;
                i=i+1;
            end
        end
end
end
