function varargout = report(BIDS, Type, Subj, Ses, Run, ReadNII)
% Create a text summary of a BIDS dataset
% FORMAT R = bids.report(BIDS, Type, Subj, Ses, Run, ReadNII)
% BIDS     - directory formatted according to BIDS [Default: pwd]
%
% Type     - Specfies the type of report {'short','article','detailed'}
%            [Default: 'short']
% Subj     - Specifies which subject(s) to use as template [Default: 1]
% Ses      - Specifies which session(s) to use as template [Default: 1]
%            Can be a vector. Set to 0 to use all sessions.
% Run      - Specifies which BOLD run(s) to use as template [Default: 1]
% ReadNII  - Specifies whether or not to read NIfTI headers [Default: true]
%
% R        - string containing report. If not specified, the text report is
%            displayed on screen.
%
% Unless specified, the function will only read the data from the first
% subject, session, and run (for each task of BOLD). This can be an issue
% if different subjects/sessions contain very different data.
%
% See also:
% bids.layout

% Copyright (C) 2018, Remi Gau
% Copyright (C) 2018--, BIDS-MATLAB developers


%-Check inputs
%--------------------------------------------------------------------------
if ~nargin
    BIDS = pwd;
end
if nargin < 2 || isempty(Type)
    Type = 'short';
end
if nargin < 3 || isempty(Subj)
    Subj = 1;
end
if nargin < 4 || isempty(Ses)
    Ses = 1;
end
if nargin < 5 || isempty(Run)
    Run = 1;
end
if nargin < 6
    ReadNII = true;
end

%-Parse the BIDS dataset, if necessary
%--------------------------------------------------------------------------
BIDS = bids.layout(BIDS);

%-Get sessions
%--------------------------------------------------------------------------
subjs_ls = bids.query(BIDS,'subjects');
sess_ls = bids.query(BIDS,'sessions', 'sub', subjs_ls(Subj));
if isempty(sess_ls)
    sess_ls = {''};
end
if Ses == 0
    Ses = 1:numel(sess_ls);
end

%-Text report
%==========================================================================
switch lower(Type)
    case 'short'
        R = text_report_short(BIDS);
    case 'detailed'
        warning('Option not yet supported.');
        R = '';
    case 'article'
        R = text_report_article(BIDS,subjs_ls,Subj,sess_ls,Ses,Run,ReadNII);
    otherwise
        error('Unknown option.');
end

%-Return, display or save as file
%--------------------------------------------------------------------------
if ~nargout
    fprintf('%s',R);
elseif false
    fid = fopen(output_file,'wt');
    if fid == -1
        error('Could not open file for writing.');
    end
    fprintf(fid,'%s',R);
    fclose(fid);
else
    varargout = { R };
end


%==========================================================================
function R = text_report_short(BIDS)
d = dict;
t = bids.query(BIDS,'types');
t = cellfun(@(x) d.datatypes(x),t,'UniformOutput',false);
t =  sprintf('%s, ',t{:});
m = bids.query(BIDS,'modalities');
m = cellfun(@(x) d.modality(x),m,'UniformOutput',false);
m = sprintf('%s, ',m{:});
R = sprintf(['BIDS dataset: %s\n',...
    '%d subjects with %d sessions and %d runs.\n',...
    'Data types are: %s\n',...
    'Modalities are: %s\n'],...
    BIDS.dir,...
    numel(bids.query(BIDS,'subjects')),...
    numel(bids.query(BIDS,'sessions')),...
    numel(bids.query(BIDS,'sessions')),...
    t(1:end-2),...
    m(1:end-2));


%==========================================================================
function R = text_report_article(BIDS,subjs_ls,Subj,sess_ls,Ses,Run,ReadNII)

R = '';

%-Scanner details
%--------------------------------------------------------------------------
% str = 'MR data were acquired using a {tesla}-Tesla {manu} {model} MRI scanner.';

%-Loop through all the required sessions
%--------------------------------------------------------------------------
for iSess = Ses
    
    types_ls = bids.query(BIDS,'types', 'sub', subjs_ls(Subj), 'ses', sess_ls(iSess));
    tasks_ls = bids.query(BIDS,'tasks', 'sub', subjs_ls(Subj), 'ses', sess_ls(iSess));
    % mods_ls = bids.query(BIDS,'modalities');
    
    for iType = 1:numel(types_ls)
        
        switch types_ls{iType}
            
            %-Anatomical
            %--------------------------------------------------------------
            case {'T1w', 'inplaneT2', 'T1map', 'FLASH'}
                
                % anat text template
                anat_tpl = fullfile(fileparts(mfilename('fullpath')),...
                    'config','anat.tpl');
                
                % get parameters
                acq_param = get_acq_param(BIDS, subjs_ls{Subj}, sess_ls{iSess}, ...
                    types_ls{iType}, '', '', ReadNII);
                
                % render template
                R = bids.internal.parse_template(anat_tpl,acq_param);
                
            %-Functional
            %--------------------------------------------------------------
            case 'bold'
                
                % func text template
                func_text = cat(2, ...
                    '%s run(s) of %s %s %s fMRI data were collected (%s slices acquired in a %s fashion; repetition time, TR= %s ms; \n', ...
                    'echo time, TE= %s ms; flip angle, FA= %s deg; field of view, FOV= %s mm; matrix size= %s; \n', ...
                    'voxel size= %s mm; multiband factor= %s; in-plane acceleration factor= %s). Each run was %s minutes in length, during which \n', ...
                    '%s functional volumes were acquired. \n\n');
                
                % loop through the tasks
                for iTask = 1:numel(tasks_ls)
                    
                    runs_ls = bids.query(BIDS,'runs','sub', subjs_ls{Subj}, 'ses', sess_ls{iSess}, ...
                        'type', 'bold', 'task', tasks_ls{iTask});
                    
                    if isempty(runs_ls)
                        % get the parameters for that task
                        acq_param = get_acq_param(BIDS, subjs_ls{Subj}, sess_ls{iSess}, ...
                            'bold', tasks_ls{iTask}, '', ReadNII);
                        
                        % compute the number of BOLD run for that task
                        acq_param.run_str = '1';
                        
                    else % if there is more than 1 run
                        % get the parameters for that task
                        acq_param = get_acq_param(BIDS, subjs_ls{Subj}, sess_ls{iSess}, ...
                            'bold', tasks_ls{iTask}, runs_ls{Run}, ReadNII);
                        % compute the number of BOLD run for that task
                        acq_param.run_str = num2str(numel(runs_ls));
                    end
                    
                    % set run duration
                    if ~strcmp(acq_param.tr,'[XXXX]') && ~strcmp(acq_param.n_vols,'[XXXX]')
                        acq_param.length = num2str(str2double(acq_param.tr)/1000 * str2double(acq_param.n_vols) / 60);
                    end
                    
                    % print output
                    fprintf(func_text,...
                        acq_param.run_str, acq_param.task, acq_param.variants, acq_param.seqs, ...
                        acq_param.n_slices, acq_param.so_str, acq_param.tr, ...
                        acq_param.te, acq_param.fa, ...
                        acq_param.fov, acq_param.ms, ...
                        acq_param.vs, acq_param.mb_str, acq_param.pr_str, ...
                        acq_param.length, ...
                        acq_param.n_vols);
                    fprintf('\n\n')
                end
                
            %-Fieldmap
            %--------------------------------------------------------------
            case 'phasediff'
                %warning('fmap not supported yet');
                
                % fmap text template
                fmap_text = cat(2, ...
                    'A %s %s field map (phase encoding: %s; %s slices; repetition time, TR= %s ms; \n',...
                    'echo time 1 / 2, TE 1/2= %s ms; flip angle, FA= %s deg; field of view, FOV= %s mm; matrix size= %s; \n',...
                    'voxel size= %s mm) was acquired %s. \n\n');
                
            %-DWI
            %--------------------------------------------------------------
            case 'dwi'
                %warning('dwi not supported yet');
                
                % dwi text template
                dwi_text = cat(2, ...
                    'One run of %s %s diffusion-weighted (dMRI) data were collected (%s  slices %s ; repetition time, TR= %s ms \n', ...
                    'echo time, TE= %s ms; flip angle, FA= %s deg; field of view, FOV= %s mm; matrix size= %s ; voxel size= %s mm \n', ...
                    'b-values of %s acquired; %s diffusion directions; multiband factor= %s ). \n\n');
                
            %-Physio
            %--------------------------------------------------------------
            case 'physio'
                %warning('physio not supported yet');
                
            %-M/EEG
            %--------------------------------------------------------------
            case {'eeg', 'meg', 'channels', 'headshape'}
                %warning('MEEG not supported yet');
                
            %-Events
            %--------------------------------------------------------------
            case 'events'
                %warning('events not supported yet');
                
        end
        
    end
    
end


%==========================================================================
function acq_param = get_acq_param(BIDS, subj, sess, type, task, run, ReadNII)
% Will get info from acquisition parameters from the BIDS structure or from
% the NIfTI files

% to return dummy values in case nothing was specified
acq_param.type = type;
acq_param.variants = '[XXXX]';
acq_param.seqs = '[XXXX]';

acq_param.tr = '[XXXX]';
acq_param.te = '[XXXX]';
acq_param.fa = '[XXXX]';

acq_param.task  = task;

acq_param.run_str  = '[XXXX]'; % number of runs (dealt with outside this function but initialized here)
acq_param.so_str  = '[XXXX]'; % slice order string
acq_param.mb_str  = '[XXXX]'; % multiband
acq_param.pr_str  = '[XXXX]'; % parallel imaging
acq_param.length  = '[XXXX]';

acq_param.for_str = '[XXXX]'; % for fmap: for which run this fmap is for.
acq_param.phs_enc_dir = '[XXXX]'; % phase encoding direction.

acq_param.bval_str = '[XXXX]';
acq_param.n_vecs = '[XXXX]';

acq_param.fov = '[XXXX]';
acq_param.n_slices = '[XXXX]';
acq_param.ms = '[XXXX]'; % matrix size
acq_param.vs = '[XXXX]'; % voxel size
acq_param.n_vols  = '[XXXX]';


%-Look into the metadata sub-structure for BOLD data
%--------------------------------------------------------------------------
if ismember(type, {'T1w' 'inplaneT2' 'T1map' 'FLASH' 'dwi'})
    
    filename = bids.query(BIDS, 'data', 'sub', subj, 'ses', sess, 'type', type);
    metadata = bids.query(BIDS, 'metadata', 'sub', subj, 'ses', sess, 'type', type);
    
elseif strcmp(type, 'bold')
    
    filename = bids.query(BIDS, 'data', 'sub', subj, 'ses', sess, 'type', type, ...
        'task', task, 'run', run);
    metadata = bids.query(BIDS, 'metadata', 'sub', subj, 'ses', sess, 'type', type, ...
        'task', task, 'run', run);
    
elseif strcmp(type, 'phasediff')
    
    filename = bids.query(BIDS, 'data', 'sub', subj, 'ses', sess, 'type', type, 'run', run);
    metadata = bids.query(BIDS, 'metadata', 'sub', subj, 'ses', sess, 'type', type, 'run', run);
    
end

%fprintf(' - %s\n', filename{1})

if isfield(metadata, 'EchoTime')
    acq_param.te = num2str(metadata.EchoTime*1000);
elseif isfield(metadata, 'EchoTime1') && isfield(metadata, 'EchoTime2')
    acq_param.te = [num2str(metadata.EchoTime1*1000) ' / '  num2str(metadata.EchoTime2*1000)];
end

if isfield(metadata, 'RepetitionTime')
    acq_param.tr = num2str(metadata.RepetitionTime*1000);
end

if isfield(metadata, 'FlipAngle')
    acq_param.fa = num2str(metadata.FlipAngle);
end

if isfield(metadata, 'SliceTiming')
    acq_param.so_str = define_slice_timing(metadata.SliceTiming);
end

if isfield(metadata, 'PhaseEncodingDirection')
    acq_param.phs_enc_dir = metadata.PhaseEncodingDirection;
end

if isfield(metadata, 'IntendedFor')
    acq_param.for_str = metadata.IntendedFor;
end

%-Try to read the relevant NIfTI file to get more info from it
%--------------------------------------------------------------------------
if ReadNII
    fprintf('  Opening file %s.\n',filename{1})
    try
        % read the header of the NIfTI file
        hdr = read_nifti_header(filename{1});
        acq_param.n_vols = hdr.dim(4); % nb volumes
        
        dim = hdr.dim;
        acq_param.n_slices = sprintf('%i', dim(3)); % nb slices
        acq_param.ms = sprintf('%i X %i', dim(1), dim(2)); %matrix size
        
        vs = abs(diag(hdr.mat));
        acq_param.vs = sprintf('%.2f X %.2f X %.2f', vs(1), vs(2), vs(3)); % voxel size
        
        acq_param.fov = sprintf('%.2f X %.2f', vs(1)*dim(1), vs(2)*dim(2)); % field of view
        
    catch
        warning('Could not read the header from file %s.\n', filename{1});
    end
end


%==========================================================================
function ST_def = define_slice_timing(SliceTiming)
% Try to figure out the order the slices were acquired from their timing
if iscell(SliceTiming)
    SliceTiming = cell2mat(SliceTiming);
end
[~,I] = sort(SliceTiming);
if all(I==(1:numel(I))')
    ST_def = 'ascending';
elseif all(I==(numel(I):-1:1)')
    ST_def = 'descending';
elseif I(1)<I(2)
    ST_def = 'interleaved ascending';
elseif I(1)>I(2)
    ST_def = 'interleaved descending';
else
    ST_def = '????';
end


%==========================================================================
function hdr = read_nifti_header(filename)
% Read NIfTI-1 header
hdr = [];
fp  = fopen(filename,'r','native'); % TODO % detect gz files
if fp==-1, return; end

sts = fseek(fp,40,'bof');
if sts==-1, fclose(fp); return; end
dim = fread(fp,8,'*int16');
fclose(fp);
if isempty(dim), return; end
if dim(1)<1 || dim(1)>7
    dim = swapbytes(dim);
end
dim = double(dim(2:(dim(1)+1)))';
dim = [dim 1 1 1 1 1];

mat = eye(4); % TODO % read qform, sform

hdr = struct('dim',dim,'mat',mat);


%==========================================================================
function d = dict
% Term dictionary
persistent D
if ~isempty(D), d = D; return; end

d.modality = containers.Map;
d.modality('anat') = 'Anatomy imaging data';
d.modality('func') = 'Task (including resting state) imaging data';
d.modality('dwi')  = 'Diffusion imaging data';
d.modality('fmap') = 'Fieldmap data';
d.modality('meg')  = 'MEG recording data';
d.modality('eeg')  = 'EEG recording data';
d.modality('ieeg') = 'iEEG recording data';
d.modality('beh')  = 'Behavioral recording data';

d.datatypes = containers.Map;
d.datatypes('T1w')         = 'T1 weighted';
d.datatypes('T2w')         = 'T2 weighted';
d.datatypes('T1rho')       = 'T1 Rho map';
d.datatypes('T1map')       = 'Quantitative T1 map';
d.datatypes('T2 map')      = 'Quantitative T2 map';
d.datatypes('T2star')      = 'T2*';
d.datatypes('FLAIR')       = 'FLAIR';
d.datatypes('FLASH')       = 'FLASH';
d.datatypes('MEFLASH')     = 'MEFLASH'; % not in specs
d.datatypes('PD')          = 'Proton density';
d.datatypes('PDmap')       = 'Proton density map';
d.datatypes('PDT2')        = 'Combined PD/T2';
d.datatypes('inplaneT1')   = 'Inplane T1';
d.datatypes('inplaneT2')   = 'Inplane T2';
d.datatypes('angio')       = 'Angiography';
d.datatypes('dwi')         = 'Diffusion weighted';
d.datatypes('bold')        = 'BOLD';
d.datatypes('cbv')         = 'CBV';
d.datatypes('phase')       = 'Phase';
d.datatypes('phase1')      = 'Phase';
d.datatypes('phase2')      = 'Phase';
d.datatypes('phasediff')   = 'Phase difference image';
d.datatypes('magnitude')   = 'Magnitude image';
d.datatypes('magnitude1')  = 'Magnitude image';
d.datatypes('magnitude2')  = 'Magnitude image';
d.datatypes('fieldmap')    = 'Fieldmap image';
d.datatypes('epi')         = 'PEpolar';
d.datatypes('coordsystem') = 'Coordinate system';
d.datatypes('events')      = 'Events';
d.datatypes('headshape')   = 'Headshape';
d.datatypes('channels')    = 'Channels';
d.datatypes('electrodes')  = 'Electrodes';
d.datatypes('eeg')         = 'EEG';
d.datatypes('ieeg')        = 'iEEG';
d.datatypes('meg')         = 'MEG';
d.datatypes('physio')      = 'Physiological recording';
d.datatypes('stim')        = 'Stimulus';
d.datatypes('photo')       = 'Photo';

d.dir = containers.Map;
d.dir('i')  = 'left to right';
d.dir('i-') = 'right to left';
d.dir('j')  = 'posterior to anterior';
d.dir('j-') = 'anterior to posterior';
d.dir('k')  = 'inferior to superior';
d.dir('k-') = 'superior to inferior';

d.seq = containers.Map;
d.seq('EP') = 'echo planar';
d.seq('GR') = 'gradient recalled';
d.seq('IR') = 'inversion recovery';
d.seq('RM') = 'research mode';
d.seq('SE') = 'spin echo';

d.seqvar = containers.Map;
d.seqvar('MP')   = 'MAG prepared';
d.seqvar('MTC')  = 'magnetization transfer contrast';
d.seqvar('NONE') = 'no sequence variant';
d.seqvar('OSP')  = 'oversampling phase';
d.seqvar('SK')   = 'segmented k-space';
d.seqvar('SP')   = 'spoiled';
d.seqvar('SS')   = 'steady state';
d.seqvar('TRSS') = 'time reversed steady state';

D = d;
