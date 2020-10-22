function increase_convergence(fname, amt)
fid = fopen(strcat(fname, '.inp'), 'r');
i = 1;
tline = fgetl(fid);
text{i} = tline;
while ischar(tline)
    i = i + 1;
    tline = fgetl(fid);
    text{i} = tline;
    if (strcmp(text{i}, "max_number_iterations"))
        target = i + 1;
    end
end
fclose(fid);
text{target} = num2str(amt);
fid = fopen(strcat(fname, '.inp'), 'w');
for i = 1:numel(text)
    if i+1 > numel(text)
        fprintf(fid, '%s', text{i});
        fclose(fid);
        break;
    else
        fprintf(fid, '%s\n', text{i});
    end
end
end