function  [result, FLAG] = export_mstm_intermediate(File)

%{
This function reads the MSTM intermediate file. The output of this file is 

result = [execution time (min), number of itterations, error or solution]

IF the simulation did not finish or there was an error the output will be

result = [NaN, NaN, NaN]

%}
FLAG = 0;
fid = fopen(strcat(File, "_output_intermediate.dat") ,'r');
if fid == -1
    disp('Intermediate file does not exist!')
    %keyboard;
    result = [NaN, NaN, NaN];
    FLAG = 1;
    return;
end
line = 1;
while ~feof(fid)
    text{line} = fgetl(fid);
    line = line+1;
end
fclose(fid);
text = text';
L = length(text);

if strfind(text{L}, 'execution time:')
    %dummy = erase(text{L}, 'execution time:     ');
    dummy = erase(text{L}, 'execution time:');
    dummy = erase(dummy, 'min');
    time = str2num(dummy);
    %if isempty(time)
    if length(time) ~= 1
        dummy = erase(dummy, 'sec');
        time = str2num(dummy) / 60;
    end
    if length(time) ~= 1
        dummy = erase(dummy, 'hours');
        time = str2num(dummy) * 60;
    end

    dummy = erase(text{L-1}, 'max iterations, soln error:');
    dummy = str2num(dummy);
    result = [time, dummy];
else
    result = [NaN, NaN, NaN];
    FLAG = 1;
    disp('Simulation didnt finish! Not converged!')
end


if length(result) ~= 3
    keyboard;
end
end




