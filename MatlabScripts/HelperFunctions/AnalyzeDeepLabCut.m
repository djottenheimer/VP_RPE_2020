function [timestamps,xcoordinates,ycoordinates]=AnalyzeDeepLabCut(file)


data=h5read(file,'/df_with_missing/table');
leftlight_x=data.values_block_0(1,:)';
leftlight_y=data.values_block_0(2,:)';
leftlight_l=data.values_block_0(3,:)';
rightlight_x=data.values_block_0(4,:)';
rightlight_y=data.values_block_0(5,:)';
rightlight_l=data.values_block_0(6,:)';
headcap_x=data.values_block_0(7,:)';
headcap_y=data.values_block_0(8,:)';
headcap_l=data.values_block_0(9,:)';

%only look at frames where deeplabcut is above a certain likelihood cutoff
likelihood_cutoff=0.95;
xcoord=headcap_x(headcap_l>likelihood_cutoff);
ycoord=headcap_y(headcap_l>likelihood_cutoff);
frames=1:length(headcap_x);
incframes=frames(headcap_l>likelihood_cutoff);
incframetimes=(incframes-1)/30;

%remove outliers
wrongx = isoutlier(xcoord,'movmedian',30);%,'samplepoints',incframetimes);
wrongy = isoutlier(ycoord,'movmedian',30);%,'samplepoints',incframetimes);

timestamps=incframetimes(wrongx==0 & wrongy==0);
xcoordinates=xcoord((wrongx==0 & wrongy==0));
ycoordinates=ycoord((wrongx==0 & wrongy==0));

end

