% image/matrix 'smoothing' using either 'nearest neighbour', 'next nearest
% neighbour', or 'gaussian'
% require the first 3 inputs

function [smoothed] = image_average(toSmooth,mode,output_flag,gaussNum)

    % toSmooth should be a m (rows) x n (columns) array
    % toSmooth(m,:) would output row m
    % toSmooth(:,n) would output column n
    
    % mode is either: NN, NNN, or gauss (case insensitive)
    
    % output_flag is just gonna be y,yes,1 or it gets ignored

% number of rows and columns of the array to smooth
row_num=size(toSmooth,1);
col_num=size(toSmooth,2);
smoothed=zeros(row_num,col_num);
weight=0.35;


%(median filter?)
if strcmpi(mode,'NN') % NEAREST NEIGHBOUR


    for i=2:row_num-1
        for j=2:col_num-1
            %sprintf("i,j=%d,%d",i,j)
            
            smoothed(i,j)=(toSmooth(i,j)+(toSmooth(i-1,j)+toSmooth(i+1,j)+toSmooth(i,j-1)+toSmooth(i,j+1))*weight)/(1+4*weight);
        end
    end

    %for first row and last row
    smoothed(1,1)=(toSmooth(1,1)+(toSmooth(2,1)+toSmooth(1,2))*weight)/(1+2*weight); %1,1 corner
    smoothed(1,col_num)=(toSmooth(1,col_num)+(toSmooth(2,col_num)+toSmooth(1,col_num-1))*weight)/(1+2*weight); %1,last corner

    for j=2:col_num-1
        smoothed(1,j)=(toSmooth(1,j)+(toSmooth(1,j+1)+toSmooth(1,j-1)+toSmooth(2,j))*weight)/(1+3*weight); %first row
        smoothed(row_num,j)=(toSmooth(row_num,j)+(toSmooth(row_num,j+1)+toSmooth(row_num,j-1)+toSmooth(row_num-1,j))*weight)/(1+3*weight); %last row
    end

    for i=2:row_num-1
        smoothed(i,1)=(toSmooth(i,1)+(toSmooth(i-1,1)+toSmooth(i+1,1)+toSmooth(i,2))*weight)/(1+3*weight); %first column
        smoothed(i,col_num)=(toSmooth(i,col_num)+(toSmooth(i-1,col_num)+toSmooth(i+1,col_num)+toSmooth(i,col_num-1))*weight)/(1+3*weight); %last column
    end

   % last,1 corner and last,last corner
    smoothed(row_num,1)=(toSmooth(row_num,1)+(toSmooth(row_num-1,1)+toSmooth(row_num,2))*weight)/(1+2*weight); %1,1 corner
    smoothed(row_num,col_num)=(toSmooth(row_num,col_num)+(toSmooth(row_num-1,col_num)+toSmooth(row_num,col_num-1))*weight)/(1+2*weight); %1,last corner


    %disp(smoothed);


    % in terms of weighting, maybe take 1/4th? from each one?
    % or 1/2? 
    % so we would get: NEW= (OLD + 1/2*left + 1/2*right + 1/2*up +
    % 1/2*down)/3
    % or NEW=(OLD + (left + right + up + down)/2)/2


elseif strcmpi(mode,'NNN') %next nearest neighbour (or rather just include diagonals?)
    %let's... not do this for now lol

elseif strcmpi(mode,'gauss')
    %use a gaussian averaging

    if exist('gaussNum','var')
        %uh use gaussNum i guess?
    else
        gaussNum=1; %using just the 4 nearest i guess?
    end


else
    disp("Invalid method")

end


%% Potentially output in a seperate figure
if isnumeric(output_flag)
    str2num(output_flag);
end

if strcmpi(output_flag,'Yes') | strcmpi(output_flag,'y') | strcmpi(output_flag,'1')
    figure();
    s=pcolor(smoothed);
    set(s,'EdgeColor','none');
    colormap('gray');

end
% otherwise just don't output if it's not one of these lol
    


%% this was the one used for the hyperspectral
% %HS_raw=HS_normed_ETref_avg;
% sizeHS=size(HS_raw,2);
% avgHS(:,1)=(HS_raw(:,1)+HS_raw(:,2))/2;
% for j=2:sizeHS-1
%     avgHS(:,j)=(HS_raw(:,j-1)+HS_raw(:,j)+HS_raw(:,j+1))/3;
% end
% avgHS(:,sizeHS)=(HS_raw(:,sizeHS-1)+HS_raw(:,sizeHS))/2;
% 
% figure();
% s=pcolor(xs,wns,avgHS);set(s,'EdgeColor','none');
% ylim([1330 1550]);
% caxis([0.7 1.8]);
% xlabel('x (um)');
% ylabel('wn (cm-1)');
% title("Oversampled Norm S2 hBN/ET to ET",'Interpreter', 'none');
% set(gcf,'Position',[450 300 600 250]);

