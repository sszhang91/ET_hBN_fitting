%% Preamble things


addpath 'D:\BASOV\Matlab scripts'; % the path to your folder with relevant matlab scripts
datafolder='.\DataFiles\221121_junhe\'; % a relative path to where the .gsf files are stored
% in a folder named DataFiles inside the same directory the .m file is in

% these are just some colormaps, and are in the Matlab scripts folder
load('sky_cmap.mat');
load('colororder_igor.mat');
%to apply the colororder_igor, do colororder(color_order_mult) before
%the plot. reset to the default colororder with colororder(default);

%% Data loading variables
temps=[50]; % 
%wns=[1450,1471,1472]; %if we just want to do some test 
wns=readmatrix("file_wns.csv");
sets=[1,1,1]; %for now (for this to be used, we could similarly 
% have a csv file that listed the number of sets for each wavenumber, and
% load them in individually. this is currently being overwritten by the
% below variables

% only using one temperature and one set for now
temp=50; %forcing it to only be 50K files
setnum=1; %and only the one labelled as _1
dir='fwd'; %and only the forward direction
chn='S4'; %and only channel S4


%% Load Data Files

% Workflow for this cell:
% create a new file load script that takes the complete file name in
% in other words, the file path is *outside* this function
% this includes any 'set' or 'fwd' options, however, the function will
% need to account for the file format itself (e.g. x-res: vs X-Res= and
% number of lines

% Create file name to open the file

% Output it into storage cell arrays: S4_cell, x_cell, y_cell, xres_cell,
% yres_cell etc. alternatively, S4_cell, 
% xrange_cell, yrange_cell, xres_cell, yres_cell
% 

% Plot if needed


% the file names right now have a format of "hBN_ET_50K_1450_S4_fwd_1.gsf"
% we can write this as "hBN_ET_[temp]K_[wn]_[S4/z/S3]_[dir]_[set].gsf"
% therefore, the file location to open will be a
% strcat(datafolder,filename)
% we will then pass to the load_file function the entire filepath
% and then get out all the different parts we need
% for now, we'll just only use set 1, and manually get the list of
% temperatures

%allocate variables
S4_cell=cell([1,length(wns)]); %cell array to store the scattering amp.
xs_cell=cell([1,length(wns)]); %cell array to store the x-axis array
ys_cell=cell([1,length(wns)]); %cell array to store the y-axis array
res_cell=cell([1,length(wns)]);%cell array to store the x, y pixle number 
dim_cell=cell([1,length(wns)]);%cell array to store the x,y real dimensions

for j=1:length(wns)
    filename=sprintf("hBN_ET_%gK_%g_%s_%s_%g.gsf",temp,wns(j),chn,dir,setnum);
    filepath=strcat(datafolder,filename) % do not suppress this print (check the file path directory)
    [S4_cell{j}, xs_cell{j}, ys_cell{j},res_cell{j},dim_cell{j}] ...
        = sp_load_file_filename(filepath); %load using this script that is in the Matlab Scripts folder
end

%% Shift each line to line them up
% Correct for physical thermal drift

% Assuming we now have at least three cell arrays of S4_cell, x_cell, and
% y_cell
% potentially also load in the z files, and do edge detection
% and then find the maximum around that location, presumably within
% 100nms(?)
% store this '0' value in an index_cell 

% for now we are just going to create a new cell array of S4_shift
% and in this one, we will store the shifted file

% if needed, we'll also plot the scans in a single tiled layout
% comment all parts out for figures if not using
fig_S4_raw=figure('Name','S4');
p=tiledlayout('flow','TileSpacing','tight','Padding','tight');


S4_shift_cell=cell([1,length(wns)]);
shift_ind_cell=cell(1,length(wns));
shift_zeros=NaN(1,length(wns)); %store the location of '0'


smoothnum=3; % in order to remove outliers and make this automatic matching
% more reliable, we do some smoothing using nearest neighbours
% note: the current script  of "image_average" is an incomplete one
% the only part that works is nearest neighbour averaging


for k=1:length(wns)
    % put the S4 into a temp array
    % put the x, y res into a temp values
    unshifted_mat=S4_cell{k};
    xyres=res_cell{k};
    shift_ind=NaN(1,xyres(2)); % this should be the number of rows
    shifted_mat=NaN(size(unshifted_mat));
    
    %get the dimensional info
    xs_raw=xs_cell{k}*1e6; %convert to um
    ys_raw=ys_cell{k}*1e6;
    xydim=dim_cell{k}*1e6; 
    
    %do the plotting here, if we want
    nexttile(p);
    s4=pcolor(xs_raw,ys_raw,unshifted_mat); %s4 is a temporary variable i'll regret eventually
    set(s4,'EdgeColor','none'); colormap(sky_cmap); %remove lines and change colormap
    pbaspect([xydim(1) xydim(2) 1]); %xydim(1) is x length, xydim(2) is y length
    set(gca,'Yticklabel',[],'Xticklabel',[]); %do not plot the numbers on axis
    title(sprintf('%gcm-1',wns(j)));
    %make a scale bar of 0.5um
    %since we know the x dimensions, and this has the correct dimensions
    %for aspect ratio: plot a line using these coordinates
    line([0.15 0.65],[xydim(2)/5 xydim(2)/5],'Color','White','LineWidth',1,'Marker','|');
    text(0.3,xydim(2)/5*2,'0.5um','Color','White','FontWeight','bold');
    
   

    %do a bit of smoothing to hopefully remove outliers
    %unshifted_smoothed=image_average(unshifted_mat,'NN','n'); %ignore other flags
    unshifted_smoothed=image_average(unshifted_mat,'NN','n');
    if smoothnum == 1
        remsooth=unshifted_smoothed; %no smoothing
    elseif smoothnum > 1 
        resmooth=image_average(unshifted_mat,'NN','n'); %do an iterative smoothing of the image
        for renum=1:smoothnum
            if renum==smoothnum
                resmooth=image_average(resmooth,'NN','n'); %so the fact that this is the same as below
            else
                resmooth=image_average(resmooth,'NN','n'); %is for there to be options for other variables (just leave for now)
            end
        end
        
    end

    % get the index of the maximum line for each row
    for m=1:xyres(2)
        %unshifted_line=unshifted_smoothed(m,:);
        unshifted_line=resmooth(m,:); % a temporary variable for the 'current' line
        [~,shift_ind(m)]=max(unshifted_line); %finding the maximum value (throw away) and the index

    end

    %shift using circshift (see matlab doc.)
    shift_zeros(k)=round(mean(shift_ind));
    for m = 1:size(unshifted_smoothed,1)      
        shifted_mat(m,:) = circshift(unshifted_smoothed(m,:),-shift_ind(m) + shift_zeros(k),2);
        %the main problem here is that the parts at the ends that get
        %wrapped around are essentially garbage
        %but since we're not using them anyway, we can ignore them
    end

    % store the indices into the cell
    shift_ind_cell{k}=shift_ind;
    S4_shift_cell{k}=shifted_mat;


    %plot it, essentially same as above but now with the shifted plots
    
    s=pcolor(shifted_mat);
    set(s,'EdgeColor','none');
    colormap('gray');

    nexttile(q);
    s4_shift=pcolor(shifted_mat); %s4 is a temporary variable i'll regret eventually
    set(s4_shift,'EdgeColor','none'); colormap(sky_cmap); %remove lines and change colormap
    pbaspect([xydim(1) xydim(2) 1]); %xydim(1) is x length, xydim(2) is y length
    set(gca,'Yticklabel',[],'Xticklabel',[]); %do not plot the numbers on axis
    title(sprintf('%gcm-1',wns(j)));
    %make a scale bar of 0.5um
    %since we know the x dimensions, and this has the correct dimensions
    %for aspect ratio: plot using dimensions
    line([0.15 0.65],[xydim(2)/5 xydim(2)/5],'Color','White','LineWidth',1,'Marker','|');
    text(0.3,xydim(2)/5*2,'0.5um','Color','White','FontWeight','bold');
    
   
    

end

%% Take line average

% now that each of these are lined up, we average each line
% we may need to clean up the data prior to this
% if we have 'NA's, just do 'mean' which should ignore nan
% this can be done with: mean(array,'omitnan') | this is only if we expand
% the array instead of doing a cric shift

xs_shift_cell=cell(size(wns)); %create a blank shell

figure(); %we actually do want to plot this
hold on;
for a=1:length(wns)
    % average all line-cuts: hopefully by this point everything looks good
    % at this point, we can possibly do an interpolation and double the
    % number of points? try with both
    shifted_mat=S4_shift_cell{a};
    averaged_line=mean(shifted_mat);
    shifted_line{a}=averaged_line; %store it
    

    % create the shifted xs
    xydim=dim_cell{a};
    xyres=res_cell{a};
    dx=xydim(1)/xyres(1);
    xs_shift=zeros(size(xs_cell{a}));

    for xs_ind=1:xyres(1) %number of x elements?
        xs_shift(xs_ind)=(xs_ind-shift_zeros(a))*dx;
    end
    xs_shift_cell{a}=xs_shift; %store it

    plot(xs_shift,averaged_line);
end
leg_text=(num2str(wns')); %makes an array of each wavenumber but as a string
legend(leg_text);
hold off;

%% Plot the colour plots seperately; do not need if already done above
% 
% fig_S4_raw=figure('Name','S4');
% p=tiledlayout('flow','TileSpacing','tight','Padding','tight');
% 
% % and to plot the shifted ones
% fig_S4_shift=figure('Name','Shifted S4');
% q=tiledlayout('flow','TileSpacing','tight','Padding','tight');
% for a=1:length(wns)
%     %retrieve values
%     xs_raw=xs_cell{a}*1e6; %convert to um
%     ys_raw=ys_cell{a}*1e6;
%     xydim=dim_cell{a}*1e6; 
%     
%     unshifted_mat=S4_cell{a};
%     shifted_mat=S4_shift_cell{a};
%     
%     
%     nexttile(p);hold on;
%     s4=pcolor(xs_raw,ys_raw,unshifted_mat); %s4 is a temporary variable i'll regret eventually
%     set(s4,'EdgeColor','none'); colormap(sky_cmap); %remove lines and change colormap
%     pbaspect([xydim(1) xydim(2) 1]); %xydim(1) is x length, xydim(2) is y length
%     set(gca,'Yticklabel',[],'Xticklabel',[]); %do not plot the numbers on axis
%     title(sprintf('%gcm-1',wns(a)));
%     %make a scale bar of 0.5um
%     line([0.15 0.65],[xydim(2)/5 xydim(2)/5],'Color','White','LineWidth',1,'Marker','|');
%     %text(0.3,xydim(2)/5*2,'0.5um','Color','White','FontWeight','bold');
%     hold off;
%     
%     xs_shift=xs_shift_cell{a}*1e6; %retrieve it
%     nexttile(q);hold on;
%     s4_shift=pcolor(xs_shift,ys_raw,shifted_mat); %s4 is a temporary variable i'll regret eventually
%     set(s4_shift,'EdgeColor','none'); colormap(sky_cmap); %remove lines and change colormap
%     pbaspect([xydim(1) xydim(2) 1]); %xydim(1) is x length, xydim(2) is y length
%     set(gca,'Yticklabel',[],'Xticklabel',[]); %do not plot the numbers on axis
%     title(sprintf('%gcm-1',wns(a)));
%     %make a scale bar of 0.5um
%     line([xs_shift(1)+0.15 xs_shift(1)+0.65],[xydim(2)/5 xydim(2)/5],'Color','White','LineWidth',1,'Marker','|');
%     %text(0.3,xydim(2)/5*2,'0.5um','Color','White','FontWeight','bold');
%     xlim([xs_shift(1) xs_shift(length(xs_shift))]);
%     hold off;
% end


%% Define a '0' distance

% up to this point we'll be using an x and y pixel size, instead of
% absolute values of x and y (i.e. all distances will be relative)
% now that we have a single line, we define the '0' value index as x=0
% and then create the line_x_units array
% potentially, we can now redimension it such that all curves will have the
% same x resolution, and extract just the parts of the averaged curve such
% that it includes 0.5 um past the edge, and 3um (values can be changed)
% down. TRY TO DO THIS USING DISTANCES INSTEAD OF NUMBER OF PIXELS
% this can also be done after, I guess?

%% Normalise each line
% Have it centered in y at a '0'
% rather than normalising the curve amplitude, we just want to be able to
% waterfall this nicely (...for now)

figure();
hold on;
for b=1:length(wns)
    to_norm=shifted_line{b};
    norm_val=mean(to_norm(20:shift_zeros(b)));
    line_normed=to_norm/norm_val;


    xs_shift=xs_shift_cell{b};
    line_normed_cell{b}=line_normed;

    % plot them each w an offset
    plot(xs_shift,line_normed+0.3*b); %the 0.3 scaling factor is just trial and error to see how to best offset it

    
end
leg_text=(num2str(wns));
legend(leg_text);
hold off;
%% Waterfall each plot
% and then we can compare where the peaks change, hopefully

% just plot some of them...
% makes it easier to handle
wns_to_plot=[1461,1465,1467,1468,1469,1470,14705,1471,14715,1472,1473,1474,1475,1478,1479,1480,1481];
figure();
hold on;
index=0;
for j=1:length(wns)

    %don't show the wonky ones
    if ismember(wns(j),wns_to_plot)==0
        continue
    end
    xs_shift=xs_shift_cell{j}*1e6;
    line_normed=line_normed_cell{j};

    % plot them each w an offset
    plot(xs_shift,line_normed+0.3*index);
    index=index+1;
end
xlim([-1.2, 0.2]);
leg_text=(num2str(wns_to_plot'));
legend(leg_text);
hold off;

