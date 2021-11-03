function X = readdat(f)
%Just the file name.
%
fileID = fopen(f);

%Skip header
IN = textscan(fileID,'%f %f %f\n',1);
%keyboard;
n = IN{1,2}; %numel(IN{1,1}); 
fgetl(fileID);

%Setup array
%X = zeros(3,n);

IN = textscan(fileID,'%f %f %f\n');
X(1,:) = IN{1,1}';
X(2,:) = IN{1,2}';
X(3,:) = IN{1,3}';
%X = flat(X);
fclose(fileID);

end
