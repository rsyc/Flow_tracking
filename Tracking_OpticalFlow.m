clc;
clear;
% --------------------------------------------------------
% Optical Flow Estimation in a Microtube Fluid Flow Video
% --------------------------------------------------------

%% Parameters
videoFile = '1xVWF_2.mov';   % filename
% corners of each of the flow channels
x1 = [218 442 442 218 218] ;
x2 = [894 1118 1118 894 894];
y1 = [50 50 841 841 50];
y2 = [50 50 841 841 50];


%% Choose optical flow method
% Options: opticalFlowHS (Horn-Schunck), opticalFlowLK (Lucas-Kanade)
opticFlow = opticalFlowRAFT;%('Smoothness',5, 'MaxIteration',100,'VelocityDifference',0);   %opticalFlowRAFT;  %opticalFlowHS;   % You can also try opticalFlowLK or opticalFlowRAFT(deep learning algorithm)

%% Process frames one by one
figure; 
for roiMask = 1:2
    % Read the video
    videoReader = VideoReader(videoFile);
    frameNum = 0;
    %vidframes = read(videoReader,[1 Inf]);
    
    % Initialize a video file to save the output resul
    Video_name = sprintf('%s%02d%s', "OpticalFlow_HS_Smoothness5_MaxIt100_VelDiff0_ROI", roiMask, '.mp4');
    vidfile = VideoWriter(Video_name,'MPEG-4');
    vidfile.FrameRate = videoReader.FrameRate;
    open(vidfile);
    
    while hasFrame(videoReader)
        frameNum = frameNum +1;
        
        % Read current frame
        frameRGB = readFrame(videoReader); 
        
        % Convert to grayscale (optical flow works on intensity, not color)
        frameGray = rgb2gray(frameRGB);

        %Create the mask specifying the size of the image.
        if roiMask == 1
            bw = poly2mask(x1,y1,size(frameGray,1),size(frameGray,2));
        else 
            bw = poly2mask(x2,y2,size(frameGray,1),size(frameGray,2));        
        end

        frameGray = uint8(bw).*frameGray;
        % Step 4: Estimate flow
        flow = estimateFlow(opticFlow, frameGray);
        
        % Step 5: Visualize results
        imshow(frameGray); hold on;
        
        % Plot flow vectors every N pixels (to avoid clutter)
        stepSize = 15;  
        plot(flow, 'DecimationFactor', [stepSize stepSize], 'ScaleFactor', 2, 'Parent', gca);
        
        title(sprintf('Optical Flow of Fluid – frame %d', frameNum));
        hold off;
        drawnow;       
        F = getframe(gca); % fig_handle is the handle to the figure, Capturing the axes avoiding borders
        writeVideo(vidfile, F);
    end
    close(vidfile)
    %reset(flowModel);
end