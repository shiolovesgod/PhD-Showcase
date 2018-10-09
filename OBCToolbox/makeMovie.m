function makeMovie(saveFN, imgFNs, fps, quality)

%the save FN should be 'filename.avi' or 'DIR/filename.avi'

%STEP 1: create the video object
myVideo = VideoWriter(saveFN);

%STEP 2: Make optional adjustments
%you can adjust the: FrameRate, VideoBitsPerPixel, VideoFormat, Video
%Compression Method, or Quality Properties

%--------------------------------------------------------------------------
if nargin > 2
    myVideo.FrameRate = fps; %Default 30
end

if nargin > 3
    myVideo.Quality = quality; %Default 75
end
%--------------------------------------------------------------------------

%STEP 3: Open the file
open(myVideo);

%Note: After you call open, you cannot modify frame rate or quality settings

%STEP 4: Write frames, still images, or an existing movie
for iFileName = imgFNs
    iFN = iFileName{1};
    
    if ischar(iFN)
        [~,~,ext] = fileparts(iFN);
        
        switch ext
            case '.tifb'
                if exist('GetTif','file')
                    %Reading Tifb images (must be between 0 and 1 for
                    %writeVideo fcn
                    img = mat2gray(GetTif(iFN)); %use imread for regular image files
                    
                else
                    error('No function found to read tifb images.  Please add the ''GetTif'' function to your currend directory.');
                end
            otherwise
                img = imread(iFN);
        end
    else
        img = mat2gray(iFN);
    end
    
    writeVideo(myVideo, img);
   
   %for intensity adjusted plots, you may want to use getframe
end

%STEP 5: Close the file

close(myVideo);

fprintf(1,'\nMovie saved in: %s\n\n', saveFN);