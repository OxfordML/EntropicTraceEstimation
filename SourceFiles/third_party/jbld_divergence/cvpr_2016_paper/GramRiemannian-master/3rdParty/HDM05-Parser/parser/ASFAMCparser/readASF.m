function skel = readASF(filename)
% skel = readASF(filename)

cl1 = clock;
%%%%%%%%%%%%%% open ASF file
fid = fopen(filename,'rt');
if fid ~= -1
    skel = emptySkeleton;
    
    idxBackSlash = findstr(filename, '\');
    if not(isempty(idxBackSlash))
        skel.filename = filename(idxBackSlash(end)+1:end);
    else
        skel.filename = filename;
    end
    
    if strncmp(skel.filename, 'TOM', 3) == 1
        % for the database of tomohiko assign another asf-filetype as
        % for the HDM05-ASF-filem because the tomokiko-files contain
        % different joints. The filetype is used in the constructNameMap
        % file where for the tomohiko asf-files a different name map has to
        % be constructed.
        skel.fileType = 'TOM.ASF';
    else
        skel.fileType = 'ASF';
    end
    
    [result1, skel] = readVersion(skel,fid);
    [result2, skel] = readName(skel,fid);
    [result3, skel] = readUnits(skel,fid);
    [result4, skel] = readDocumentation(skel,fid);
    [result5, skel] = readRootASF(skel,fid);
    [result6, skel] = readBonedata(skel,fid);
    [result7, hierarchy] = readHierarchyASF(fid);
    [result8, skel] = readSkin(skel,fid);
    if ~(result1 & result2 & result4 & result5 & result6 & result7)
        fclose(fid);
        error(['ASF: Error parsing file ' filename '!']);
    end
    
    fclose(fid);        
else
    error(['Could not open ASF file "' filename '"!']);
end

%%%%%%%%% build skeleton data structure
[result,skel] = constructSkeleton(skel, 1, hierarchy, 1); 
if ~result
    error(['ASF: Error constructing skeleton topology from hierarchy data in ' filename '!']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% traverse kinematic chain to generate a minimal edge-disjoint path covering used in the rendering process
skel = constructPaths(skel, 1, 1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% traverse kinematic chain filling in skel's "name" and "animated/unanimated" arrays
skel = constructAuxiliaryArrays(skel, 1); 
skel.animated = sort(skel.animated);
skel.unanimated = sort(skel.unanimated);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

skel.nameMap = constructNameMap(skel);

t=etime(clock,cl1);
disp(['Read ASF file ' filename ' in ' num2str(t) ' seconds.']);


