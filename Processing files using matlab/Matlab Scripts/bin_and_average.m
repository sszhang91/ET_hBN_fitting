function [avged_array]=bin_and_average(data_in);

sizeData1=size(data_in,1); %the number of rows
sizeData2=size(data_in,2); %the number of columns
avgD=ones(sizeData1,sizeData2); %initiate as an array of ones with sizeData1 rows, and sizeData2 columns

avgD(:,1)=(data_in(:,1)+data_in(:,2))/2; %first column avg of first 2 col
avgD(1,:)=(data_in(1,:)+data_in(1,:))/2; %first row avg of first 2 rows
for j=2:sizeData1-1 %this is the number of columns
    for k=2:sizeData2-1 %and going down through every row in that column
%         avgD(j,k)= mean([data_in(j,k)...
%             data_in(j,k-1) data_in(j,k+1) ...
%             data_in(j-1,k-1) data_in(j-1,k+1) ...
%             data_in(j+1,k-1) data_in(j+1,k+1) ...
%             data_in(j-1,k) data_in(j+1,k)],'all');

        avgD(j,k)= mean([data_in(j,k)...
            data_in(j,k-1) data_in(j,k+1) ...
            data_in(j-1,k) data_in(j+1,k)],'all');
    end
end
avgD(sizeData1,:)=(data_in(sizeData1-1,:)+data_in(sizeData1,:))/2;
avgD(:,sizeData2)=(data_in(:,sizeData2-1)+data_in(:,sizeData2))/2;

avged_array=avgD;