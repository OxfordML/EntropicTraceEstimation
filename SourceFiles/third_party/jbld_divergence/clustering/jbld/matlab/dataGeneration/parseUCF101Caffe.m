% parse caffe feature

function data = parseUCF101Caffe(dataPath)

N = 9216;
maxNumOfVideos = 20;

data = cell(1, maxNumOfVideos);
count = 1;
aPathList = dir(dataPath);
for ai = 1:length(aPathList)
    if isInvalid(aPathList(ai)), continue; end
    actionPath = aPathList(ai).name;
    vPathList = dir(fullfile(dataPath, actionPath));
    for vi = 1:length(vPathList)
        if isInvalid(vPathList(vi)), continue; end
        videoPath = vPathList(vi).name;
        file = dir(fullfile(dataPath, actionPath, videoPath, '*.mat'));
        feat = zeros(N, length(file));
        for fi = 1:length(file)
            load(fullfile(dataPath, actionPath, videoPath, file(fi).name));
            feat(:, fi) = deepcaffe;
        end
        data{count} = feat;
        count = count + 1;
    end
end

data(count:end) = [];

end

function invalid = isInvalid(p)

invalid = ( strcmp(p.name, '.') || strcmp(p.name, '..') || ...
    strcmp(p.name, '.DS_Store') || ~(p.isdir) );

end