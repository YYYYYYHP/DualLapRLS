function save_results(filename,results)
fid = fopen(filename,'wt');
[m,n] = size(results);
for i = 1:1:m
    for j = 1:1:n
        if j==n
            fprintf(fid,'%f\n',results(i,j));
        else
            fprintf(fid,'%f,',results(i,j));
        end
    end
end
fclose(fid)