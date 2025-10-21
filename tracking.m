clear;
clc;
%Read the movie and get all the frames.
vidObj = VideoReader("1xVWF_2.mov");
vidframes = read(vidObj,[1 Inf]);

%play the video
%videoPlayer = vision.VideoPlayer;
%while hasFrame(vidObj)
%   frame = readFrame(vidObj);
%   videoPlayer(frame);
%end

% change the RGB frame to gray scale
numFrames = size(vidframes, 4); 
grayFrames = zeros(size(vidframes, 1), size(vidframes, 2), numFrames); % Pre-allocate for grayscale frames
for i = 1:numFrames
    grayFrames(:,:,i) = rgb2gray(vidframes(:,:,:,i));
end

% Method 2: Using im2gray (more efficient)
%grayFrames = im2gray(vidframes);

% Display the first frame as an example
imshow(grayFrames(:,:,150),[0,150])

%Recover the motion field with SPTRACK.
param = [];
x1 = [218 442 442 218 218] ;
x2 = [894 1118 1118 894 894];
y1 = [50 50 841 841 50];
y2 = [50 50 841 841 50];
%Create the mask specifying the size of the image.
bw1 = poly2mask(x1,y1,size(grayFrames,1),size(grayFrames,2));
bw2 = poly2mask(x2,y2,size(grayFrames,1),size(grayFrames,2));
%Display the mask, drawing a line around the polygon.
figure
imshow(bw1+bw2)   %(bw1.*grayFrames(:,:,150)+grayFrames(:,:,150).*bw2,[0,150])
hold on
plot(x1,y1,'b','LineWidth',2)
hold on
plot(x2,y2,'b','LineWidth',2)
hold off

count = 0;
param.winsize = [64 64; 32 32; 16 16];
pad = max(param.winsize(:))/2;
param.iminc = 1; % image increment
for roiMask = {bw1, bw2}
    param.ROI = roiMask{1};    
    %safeROI = imerode(param.ROI, strel('rectangle', [pad pad]));
    %param. ROI = safeROI;
    %param.ROI = logical(bw1+bw2);
    

    count = count + 1;
    Video_name = sprintf('%s%02d%s', "speckeltracking16smoothHalf_ROI", count, '.mp4');
    vidfile = VideoWriter(Video_name,'MPEG-4');
    vidfile.FrameRate = vidObj.FrameRate;
    open(vidfile);

    Uf = []; Vf = [];  % for temporal smoothing (optional)
    alpha = 0.5;
    
    for i=1:1:size(grayFrames,3)-4
        %% Gentle denoising + contrast normalization.
        %G = grayFrames(:,:,i:i+4);
        % G = imgaussfilt3(G, 1);                  % light 3D (x,y,t) smoothing
        % for k=1:size(G,3)
        %     G(:,:,k) = adapthisteq(uint8(G(:,:,k)),'clipLimit',0.02,'Distribution','rayleigh');   % stabilize local contrast
        % end
        % [Di, Dj, id, jd] = sptrack(G, param);

        %% Build a small 5-frame stack with light smoothing + contrast norm
         G = grayFrames(:,:,i:i+4);
         G = imgaussfilt3(G, 0.6);
         for k=1:size(G,3), G(:,:,k) = adapthisteq(uint8(G(:,:,k))); end
    
        %% Track
        [Di, Dj, id, jd] = sptrack(G, param);   % Di=row, Dj=col displacements
        % [Di,Dj,id,jd] = sptrack(grayFrames(:,:,i:i+4),param); 
        % Try arrows as-is first:
        U = Dj; V = Di;
        % If directions look flipped, use:
        % U = -Dj; V = -Di

        %% Optional: temporal smoothing
        if isempty(Uf), Uf = Dj; Vf = Di; else
            Uf = alpha*Dj + (1-alpha)*Uf;
            Vf = alpha*Di + (1-alpha)*Vf;
        end
            
        %% Outlier rejection & thinning (de-noise)
        U = Uf; V = Vf;
        mag = hypot(U, V);
        minMag = 0.1; % tune for the pixels/frame
        keep = mag >= minMag;
        
        % Median-based outlier rejection (robust)
        u_med = medfilt2(U, [3 3], 'symmetric' );
        v_med = medfilt2(V, [3 3], 'symmetric');
        dev = hypot(U - u_med, V - v_med);
        th = prctile(dev(keep), 75); % robust threshold
        keep = keep & (dev <= max(th, 0.05));
        
        %% Optional: decimate vectors so the plot isn't cluttered
        stride = 2;
        K = false(size(keep));
        K(1:stride:end, 1:stride:end) = true;
        keep = keep & K;
        
        % Apply mask before plotting
        U = U(keep); V=V(keep); XX=jd(keep);YY=id(keep);
        
        %% plotting:
        %Display the motion field.
        imshow(grayFrames(:,:,i), []); hold on
        set(gca,'YDir','reverse');        % same convention as sptrack docs (axis ij)
        h = quiver(XX, YY, U, V, 2, 'r', 'LineWidth', 1, 'MaxHeadSize', 3);
            
        %image(grayFrames(:,:,i));
        %colormap gray
        %hold on
        %h = quiver(jd,id,Dj,Di,3,'r');
        %set(h, 'LineWidth',1)
        title(sprintf('Speckle flow – frame %d', i));
        %axis equal off ij;
        F = getframe(gca); % fig_handle is the handle to the figure, Capturing the axes avoiding borders
        % img_data = F.cdata; % Extract the image data from the captured frame
        % im = im2double(img_data); %gcf;
        %imshow(im)
        writeVideo(vidfile, F);
        hold off
    end
    close(vidfile)
end