function [x,y]=import_spectra(filename)
%% Import function for ????- spectra
fid = fopen(filename); 
text_line = fgetl(fid);
% test if file is empty
if text_line == -1
    fclose(fid);
    error('Datei war leer und wurde geschlossen.')
end



cell = textscan(fid,'%f64%f64','Delimiter',' ','MultipleDelimsAsOne',1,'headerlines',15);
x=cell2mat(cell(1,1));
y=cell2mat(cell(1,2));


fclose('all');