function [im_data, x, y,dx,dy] = sp_load_file_for_script_cryo_filename(folder, filename)

%file name template is, e.g. TaS2_190K_s3_bwd.gsf
%filebasename is TaS2_190K
%suffix is z, p, s3
%dir is fwd, bwd

%COMPLETELY FUCKED THIS UP SO THIS IS LEGIT JUST FOR THE CRYO GSF FILES
%THAT WE PROCESSED
%AND TAKES IN THE ENTIRE FILE NAME FOR NOW


%read in the header which is in the first few lines to get the resolution
%and then the header is a constant number so skip to next ones

%filepath=folder;%delete this 
%disp(filepath);

%we're just gonna take in THE ENTIRE FUCKING FILE NAME
% 
% if exist('extension','var')
%     filepath = strcat(folder, "\", filebasename, "-", suffix, "-", dir, ".",extension);
% else
%     
%     filepath = strcat(folder, "\", filebasename, "-", suffix, "-", dir, ".asc");
% end
folder='C:\Data\Github\ET_hBN\DataFiles\221121_junhe'
filename="hBN_ET_50K_1450_S4_fwd_1"
filepath=strcat(folder,"\",filename,".gsf");
disp(filepath)

% Open the file
file_ID = fopen(filepath,'r','l');

%skip everything until the last '=' and then place the cursor back
flag=1;
line=0;
lineend=fgets(file_ID);
%throwaway=fgetl(file_ID);
while flag > 0
    bytesprev=ftell(file_ID);
    throwaway=fgetl(file_ID);
    flag=strfind(throwaway,' = ');
    line=line+1;
    
    if (strfind(throwaway, 'XRes ='))>0
        idx=strfind(throwaway,'=');
        xres=str2num(throwaway(idx+1:end));
    elseif (strfind(throwaway, 'YRes ='))>0
        idx=strfind(throwaway,'=');
        yres=str2num(throwaway(idx+1:end));
    elseif (strfind(throwaway, 'XReal ='))>0
        idx=strfind(throwaway,'=');
        xext=str2num(throwaway(idx+1:end));
    elseif (strfind(throwaway, 'YReal ='))>0
        idx=strfind(throwaway,'=');
        yext=str2num(throwaway(idx+1:end));
    end
    
end
%return to the start of the previous line

% set x/y axis
dx=xext/xres;
dy=yext/yres;

x = (dx:dx:xext)';
y = (dy:dy:yext)';
im_data=zeros(xres,yres);

% Now use ftell to find the exact position of the end of the header
header_end = bytesprev;

% Skip any missing bits
if(sum(ismember(throwaway, char([13])))==0)
    header_end_real = ceil(header_end/4)*4; %because of CRLF issues... 
    if mod(header_end,4)==0
        header_end_real=header_end_real+4;
    end
else
    header_end_real = ceil(header_end/8)*8 ;%this works for CRLF i guess
end

% Now skip the correct number of bytes
bytesToSkip = int32(header_end_real);
fseek(file_ID, bytesToSkip, 'bof');

% % Now, additionally, skip past (sliceNumber - 1) slices
% % so that we end up at the beginning byte of the slice that we want.
% bytesToSkip = bytesToSkip + int32(x_size * y_size *
% stHeader.BytesPerVoxel * (sliceNumber - 1));
% fseek(fileHandle, bytesToSkip, 'bof');

im_data = fread(file_ID,[xres,yres],'float32','l');

im_data=flipud(im_data');

dim_out=[xext, yext];

fclose('all');