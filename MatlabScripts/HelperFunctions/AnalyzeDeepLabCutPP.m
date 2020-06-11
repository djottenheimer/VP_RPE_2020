function [timestamps,xcoordinates,ycoordinates]=AnalyzeDeepLabCutPP(h5file,videofile,ffmpeg_path)

data=h5read(h5file,'/df_with_missing/table');
headcap_x=data.values_block_0(1,:)';
headcap_y=data.values_block_0(2,:)';
headcap_l=data.values_block_0(3,:)';

%only look at frames where deeplabcut is above a certain likelihood cutoff
likelihood_cutoff=0.95;
xcoord=headcap_x(headcap_l>likelihood_cutoff);
ycoord=headcap_y(headcap_l>likelihood_cutoff);

%get timestamps of frames from video using downloaded function and ffmpeg
ts = videoframets(ffmpeg_path,videofile);

incframetimes=ts(headcap_l>likelihood_cutoff);


%remove outliers
wrongx = isoutlier(xcoord,'movmedian',30);%,'samplepoints',incframetimes);
wrongy = isoutlier(ycoord,'movmedian',30);%,'samplepoints',incframetimes);

timestamps=incframetimes(wrongx==0 & wrongy==0);
xcoordinates=xcoord((wrongx==0 & wrongy==0));
ycoordinates=ycoord((wrongx==0 & wrongy==0));

end

