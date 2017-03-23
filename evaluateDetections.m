function detResults=evaluateDetections(seqmap,resDir,dataDir)
%% evaluate Detections using Dollars toolbox
% concatenate ALL sequences and evaluate as one!
%
% SETUP:
%
% define directories for tracking results...
% resDir = fullfile('res','data',filesep);
% ... and the actual sequences
% dataDir = fullfile('..','data','MOT17Det','train',filesep);
%
%

addpath(genpath('.'));

% Evaluation parameters
cls = [2,7,8,12]; %% ambiguous classes
minvis = 0.5; % min. visibility of a person
ref=0:.1:1; % reference recall
showEachRef=1;

% read sequence map
seqmapFile=fullfile('seqmaps',seqmap);
allSeq = parseSequences(seqmapFile);

fprintf('Sequences: \n');
disp(allSeq')


% Find out the length of each sequence
% and concatenate ground truth
gtInfoSingle=[];
gtAll={};
detInfoSingle=[];
detAll={};
seqCnt=0;
allFrCnt=0;
evalMethod=1;
for s=allSeq
    seqCnt=seqCnt+1;
    seqName = char(s);
%     [seqName, seqFolder, imgFolder, imgExt, F, dirImages] ...
%         = getSeqInfo(seqName, dataDir);
    [seqName, seqFolder, imgFolder, frameRate, F, imWidth, imHeight, imgExt] ...
        = getSeqInfoFromFile(seqName, dataDir);
    
    assert(isdir(seqFolder),'Sequence folder %s missing',seqFolder);
    gtFile = fullfile(dataDir,seqName,'gt','gt.txt');
    gtRaw = dlmread(gtFile);
    

    % if something (a result) is missing, we cannot evaluate this tracker
    resFile = fullfile(resDir,'data',[seqName '.txt']);
    if ~exist(resFile,'file')
        fprintf('WARNING: results for %s not complete, no result file\n', seqName);
        evalMethod=0;
        break;
    end
        
    % if MOT16, preprocess (clean)
    if ~isempty(strfind(seqFolder,'MOT16'))
        resFile = preprocessResult(resFile, seqName,dataDir);
    end
        
    detRaw=dlmread(resFile);
    
    % 
    gtOne= {};
    detOne = {};
    for t=1:F
        allFrCnt=allFrCnt+1;
        
        exgt=find(gtRaw(:,1)==t & gtRaw(:,8)==1 & gtRaw(:,9)>=minvis);
        gtAll{allFrCnt}=[gtRaw(exgt,3:6) zeros(length(exgt),1)];
        gtOne{t}=[gtRaw(exgt,3:6) zeros(length(exgt),1)];

        exdet=find(detRaw(:,1)==t);
        bbox=detRaw(exdet,3:7);
        detAll{allFrCnt}=bbox;
        detOne{t}=bbox;

    end    

    allFgt(seqCnt) = F;    
    gtInfoSingle(seqCnt).gt=gtOne;
    detInfoSingle(seqCnt).det = detOne;
    
end


detResults=[];

mcnt=1;


fprintf('Evaluating ... \n');


% evalMethod=1;

% flags for entire benchmark
% if one seq missing, evaluation impossible
assert(evalMethod==1,'Evaluation not possible. Missing results?');
eval2D=1;

seqCnt=0;

% iterate over each sequence
for s=allSeq

    seqCnt=seqCnt+1;
    seqName = char(s);

    fprintf('\t... %s\n',seqName);


    gt0=gtInfoSingle(seqCnt).gt;
    dt0=detInfoSingle(seqCnt).det;
    [gt,dt]=bbGt('evalRes',gt0,dt0);
    [rc,pr,scores,refprcn] = bbGt('compRoc',gt,dt,0,ref);
    
    AP = mean(refprcn);
    detResults(mcnt).mets(seqCnt).rc=rc;
    detResults(mcnt).mets(seqCnt).pr=pr;
    detResults(mcnt).mets(seqCnt).ref=refprcn;
    detResults(mcnt).mets(seqCnt).AP=AP;
    detResults(mcnt).mets(seqCnt).name=seqName;
    
    fprintf('Recall:    ')
    for r=1:showEachRef:length(ref)
        fprintf('%6.3f',ref(r)); 
    end
    fprintf('\n')
    fprintf('Precision: ')
    for r=1:showEachRef:length(ref)
        fprintf('%6.3f',refprcn(r)); 
    end
    fprintf('\n');
    fprintf('Average Precision: %.2f\n',AP);      
    fprintf('\n');

    
end

gt0=gtAll;
dt0=detAll;
[gt,dt]=bbGt('evalRes',gt0,dt0);
[rc,pr,scores,refprcn] = bbGt('compRoc',gt,dt,0,ref);    
AP=mean(refprcn);
detResults(mcnt).rc=rc;
detResults(mcnt).pr=pr;
detResults(mcnt).ref=refprcn;
detResults(mcnt).AP=AP;

fprintf('\n');
fprintf(' ********************* Your Benchmark Results (2D) ***********************\n');
fprintf('Recall:    ')
for r=1:showEachRef:length(ref)
    fprintf('%6.3f',ref(r)); 
end
fprintf('\n')
fprintf('Precision: ')
for r=1:showEachRef:length(ref)
    fprintf('%6.3f',refprcn(r)); 
end
fprintf('\n');
fprintf('Average Precision: %.2f\n',AP);

