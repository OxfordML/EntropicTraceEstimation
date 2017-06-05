% saveToVideo

% feat = data.features{183};
load ../../expData/MHAD_data_whole;
feat = data{25+55};

% Prepare the new file.
vidObj = VideoWriter('handWave2.avi');
set(vidObj,'FrameRate',10);
open(vidObj);

% Create an animation.
% Z = peaks; surf(Z);
% show_skel_MSRA(feat(:,1));
clf;
show_skel_MHAD(feat(:,1));
axis tight
set(gca,'nextplot','replacechildren');
nFrames = size(feat,2);

for k = 2:nFrames
%     surf(sin(2*pi*k/20)*Z,Z)
%     show_skel_MSRA(feat(:,k));
    show_skel_MHAD(feat(:,k));
    
    % Write each frame to the file.
    currFrame = getframe(gca);
    writeVideo(vidObj,currFrame);
end

% Close the file.
close(vidObj);