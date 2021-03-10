function Generate_txt( data, Filename)

[rows, cols] = size(data);
fid = fopen( [Filename '.txt'], 'wt' );

%% Generate TXT
for j = 1:rows
    for i = 1:cols
        fprintf(fid,'%f\t',data(j,i));
    end
    fprintf(fid,'\n');
end

fclose(fid);
end