% [File, FileName, PathName] = get_file();

fid = fopen(File,'r');

txt = str2num(fgetl(fid));
RHScond = txt(1);
nord = txt(2);
nordRHS = txt(3);

txt = str2num(fgetl(fid));
nT = txt(1);
xv = txt(2);



data_lkqnmp = [];
l = [];
k = [];
q = [];
n = [];
m = [];
Qe_l = [];
Qa_l = [];
Qs_l = [];

% Grab the data from the txt file and create a matrix that has columns 
% [output order, output degree, output mode, input order, input degree, input mode, value] 
tic
while ~feof(fid)
txt = fgetl(fid);
txt = str2num(txt);
if length(txt) == 3
    l = [l,txt(1)];
    k = [k,txt(2)];
    q = [q,txt(3)];
end
if length(txt) == 6
    n = [n,txt(1)];
    m = [m,txt(2)];
    data_lkqnmp = [data_lkqnmp; [l(end),k(end),q(end),n(end),m(end),1,txt(3)+1i.*txt(4)]];
    data_lkqnmp = [data_lkqnmp; [l(end),k(end),q(end),n(end),m(end),2,txt(5)+1i.*txt(6)]];
    
%     if RHScond == 0 && n(end) < l(end)
%         ikm = (-1).^(m(end)+k(end));
%         data_lkqnmp = [data_lkqnmp; [l(end),-k(end),1,n(end),-m(end),q(end),(txt(3)+1i.*txt(4)).*ikm]];
%         data_lkqnmp = [data_lkqnmp; [l(end),-k(end),2,n(end),-m(end),q(end),(txt(5)+1i.*txt(6)).*ikm]];
%         k = [k,-k(end)];
%         m = [m,-m(end)];
%         q = [q,1];
%         q = [q,2];
%     end
     
     
end
if length(txt) == 4
    Qe_l = [Qe_l,txt(2)];
    Qa_l = [Qa_l,txt(3)];    
    Qs_l = [Qs_l,txt(4)];    
end



end

fclose(fid);
toc


%% 
% Rearrange the data to store in a large 2x2 matrix for viewing.
tic
clearvars T EH_map
row = 1;
col = 1;
counter = 0;
max_order = 1;
for q = 1:2
    for l = 1:max_order  
        for k = -l:l
            for p = 1:2   
                for n = 1:max_order
                    for m = -n:n
            
                       dummy = ismember(data_lkqnmp(:,1:6),[l,k,q,n,m,p],'rows');
                       location = find(dummy == 1);
                       if length(location)>1
                           % The code is signaling that there is a repeated
                           % data set. I.e., multiple locations in the
                           % array have the same l,k,q,n,m,p 
                            keyboard;
                       elseif isempty(location) 
                           % To save file space, the MSTM code relies on 
                           % known symetery relations. You need to 
                           % repopulate those values to represent the 
                           % entire T-matrix.
                           
                            % Find the symmetric data location
                            dummy = ismember(data_lkqnmp(:,1:6),[n,m,p,l,k,q],'rows');
                            sym_location = find(dummy == 1);
                            
                            
                            %dummy = ismember(data_lkqnmp(:,1:6),[n,m,p,l,k,q],'rows');
                            %sym_location = find(dummy == 1);
                            
                            
                            
                            
                            if isempty(sym_location)
                                % ERROR! The data and/or its symmetric
                                % analog can not be found.
                                keyboard;
                            else 
                                % Populate the symmetric version of the
                                % data ** NEED TO CHECK IF CORRECT!!!!!!
                                T(row,col) =  (data_lkqnmp(sym_location,7)).*(-1)^(m+k);
                                EH_map(row, col) = map_EH(q,p);                                 
                            end

                       else
                            T(row,col) = data_lkqnmp(location,7);
                            EH_map(row, col) = map_EH(q,p);                           
                       end
                        row = row+1;
                        
                    end
                end
            end
            col = col+1;
            row = 1;
        end
    end
end
toc
%%
figure,
imagesc(log10(abs(T)))
%imagesc(EH_map)
caxis([-2,0])
colormap(flipud(hot))
colorbar



%%
function dummy = map_EH(q,p)
                                
% Create a matrix which maps the EE, EH,
% HE, and HH regions of the T-matrix
switch q
    case 1
        switch p
            case 1
                dummy = 1;
            case 2
                dummy = 2;
        end
    case 2
        switch p
            case 1
                dummy = 3;
            case 2
                dummy = 4;
        end
end

end




