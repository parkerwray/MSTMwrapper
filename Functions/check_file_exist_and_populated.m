function FLAG = check_file_exist_and_populated(file)

FLAG = 0;
if ~isfile(file)
     FLAG = 1;
     disp('File does not exist!')
else

    line = 1;
    fid = fopen(file,'r');
    text = fgetl(fid);
    fclose(fid);
    if text == -1
        FLAG = 1;
        disp('File is empty!')
    end
end


end












