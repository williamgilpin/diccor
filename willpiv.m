% DICCOR
%
% An analysis routine for DIC analysis of microscopy data that generates a
% displacement field and a corresponding strain map.
%
% Dependencies
% ------------
% All code tested on MATLAB R2013a or greater
% Image Processing Toolbox
%
%
% This code was developed by the Brangwynne Group at Princeton University.
% If modifying or using this code, please cite the following article:
%
% W.Gilpin, S. Uppaluri, C. Brangwynne ?Worms under pressure: 
% bulk mechanical properties of C. elegans are independent of the cuticle? 
% Biophysical Journal, 2015.
%
% The "disp_on_blocks.m" MATLAB function has an unknown origin. If you are
% aware of the source of this file please let the authors know so that we
% can provide appropriate attribution.

%% calculate displacements in a set of images
close all;
pics=dir('sample_images/*.tif');          % just look in working directory
N = numel(pics);% N=500;
subpix=1;           % greater than 0 is best
%  blocksize=68/4;
blocksize=17;
%blocksize=23; % has to be odd, ideally grains move less than block radius from frame to frame
stdcut=.01;         % b/w 0 and 1, larger is more stringent requirement for coarse structure
im2=im2double(mat2gray(imread(['sample_images/' pics(1).name])));
im2=im2(:,:,1);
imagesize=size(im2);
xsize=floor(imagesize(2)/blocksize);
ysize=floor(imagesize(1)/blocksize);
x=0.5*(blocksize)+blocksize*(0:(xsize));
y=0.5*(blocksize)+blocksize*(0:(ysize));
[x y]=meshgrid(x,y);




%% the slow step where I build displacement array
vel=zeros(N,ysize+1,xsize+1,2);
for ii=1:(N-1)
    im1=im2double(mat2gray(imread(['sample_images/' pics(ii).name])));
    im2=im2double(mat2gray(imread(['sample_images/' pics(ii+1).name])));
    im1=im1(:,:,1);
    im2=im2(:,:,1);
    [bx,by,corsh]=disp_on_blocks(im1,im2,blocksize,subpix);
    stdim1 = blkproc(im1,[blocksize blocksize],'std2(x)');
    stdim2 = blkproc(im2,[blocksize blocksize],'std2(x)');
    stdgeo=(stdim1.*stdim2).^(0.5);
    k=find(stdgeo < stdcut);
    bx(k)=0;
    by(k)=0;
    %test=size(bx);
    vel(ii,:,:,1)=bx;
    vel(ii,:,:,2)=by;
    disp(ii);
end;



%% Check for outliers and use this to guess the upper and lower limits
% of acceptable velocities
alldat=vel(:,:,:,:);
figure(); hist(alldat(:),400)



%% manually remove outliers from vel and subtract mean
% topcut=-.35; %for 10x image
% botcut=-.65;
topcut=1; %60x images
botcut=-1;
blkarea=[60 60];
fun=@(x) nanmean(x(:));
for ii=1:N
    for jj=1:2      % do both components
        bb=squeeze(vel(ii,:,:,jj));
        nzerobb=(bb~=0);
        bb(~nzerobb)=NaN;
        loc_otls=(bb>topcut)|(bb<botcut);
        bb(logical(loc_otls.*nzerobb))=NaN;
        meds=nlfilter(bb,blkarea,fun);       
        bb=bb-meds;
        vel(ii,:,:,jj)=bb;
    end
end



%% spatially smooth field
sz_smooth = 3;
for ii=1:N
    for jj=1:2      % do both components
        bb=squeeze(vel(ii,:,:,jj));          %
        g_ker = fspecial('gaussian',[sz_smooth sz_smooth],2);
        bb=conv2(bb,g_ker,'same');
        vel(ii,:,:,jj)=bb;
    end
end



%% temporally average and save images of average velocity fields

window=2; % the number of frames to average into velocity field
iter=1; % number of frames to skip when generating image sequence (reduces frame rate)
dosave=1; % toggle to save images
dsamp=0;
close all;
figvec=figure();
figstr=figure();

if dosave
    mkdir(['sample_output/vec_field']);
    mkdir(['sample_output/strain_map']);
end

for ii=(floor(window/2)+1):iter:(N-floor(window/2))
    disp(ii);
    swapfig(figvec);
    pic=im2double(mat2gray(imread(['sample_images/' pics(ii).name])));
    imshow(imresize(pic,1));
    hold on;
    range=(ii-floor(window/2)):(ii+floor(window/2));
    bx=squeeze(nanmedian(vel(range,:,:,1)));
    by=squeeze(nanmedian(vel(range,:,:,2)));
    % don't display all vectors
    if dsamp==1
        sf=2;
        x2=downsample(downsample(x,sf)',sf);
        y2=downsample(downsample(y,sf)',sf);
        bx=downsample(downsample(bx,sf)',sf);
        by=downsample(downsample(by,sf)',sf);
    else
        x2=x;
        y2=y;
    end

    %bx(end-1:end,:)=NaN; % adjust for each image
    %by(end-1:end,:)=NaN; % adjust for each image
    
    quiver(x2,y2,bx,by,3,'LineWidth',3,'AutoScale','Off')
    axis off
    if dosave
        export_fig(fullfile('sample_output/vec_field',['overlay' num2str(ii)]),gcf)
    end
    
    % find the Cauchy strain %div = divergence(x,y,bx,by);
    e_xx=diff(bx,1,1); e_xx(end+1,:)=0;
    e_xy=diff(bx,1,2); e_xy(:,end+1)=0;
    e_yx=diff(by,1,1); e_yx(end+1,:)=0;
    e_yy=diff(by,1,2); e_yy(:,end+1)=0;
    e_xy=.5*(e_xy+e_yx);
    e_yx=e_xy;
    ee_size=size(e_xx);
    allc=cat(2,e_xx(:), e_xy(:),e_yx(:), e_yy(:));
    div=zeros(length(allc),1);
    for jj=1:length(allc)
        eig1=.5*((allc(jj,1) +  allc(jj,4))-sqrt((allc(jj,1)-allc(jj,4))^2+4*allc(jj,2)*allc(jj,3)));
        eig2=.5*((allc(jj,1) +  allc(jj,4))+sqrt((allc(jj,1)-allc(jj,4))^2+4*allc(jj,2)*allc(jj,3)));
        div(jj)=max(abs([eig1 eig2]));
        if isnan(eig1)|isnan(eig2)
            div(jj)=NaN;
        end
    end
    div=reshape(div,ee_size(1),ee_size(2));
    div=real(div);  % remove local rotation
    
    hold off;
    
    swapfig(figstr);
    where_nan=isnan(div);       % set so that nans are middle of map
    div(where_nan)=0;
    %div(:,end-1:end)=0; % adjust for each image
    %div(end-6:end,:)=0; % adjust for each image
    imagesc(div)
    colormap(parula(1024))
    if ii==(floor(window/2)+1)
        ord_div=sort(div(:));
        bound=ord_div(end-10);  % change this range to reduce saturation
    end
    caxis([-bound bound])
    axis off
    
    pbaspect([xsize ysize 1]);
    if dosave
        export_fig(fullfile('sample_output/strain_map',['map' num2str(ii)]),gcf)
    end
    if ~dosave
        pause(.4)
    end
end