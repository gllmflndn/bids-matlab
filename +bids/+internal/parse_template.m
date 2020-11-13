function tpl = parse_template(tpl,data,tplist)
% Data-driven template for generating textual output
% FORMAT tpl = parse_template(tpl,fm,tplist)
% tpl    - template filename or scalar cellstr containing template text
% data   - data structure
% tplist - template dictionary
%
% Implement a small subset of the Go template package:
% https://golang.org/pkg/text/template/

% Copyright (C) 2020 Guillaume Flandin
% Copyright (C) 2020-- BIDS-MATLAB developers

%-Load template
%==========================================================================
if ischar(tpl)
    fid = fopen(tpl,'rt','native','UTF-8');
    if fid == -1
        error('Cannot open template file "%s".',tpl);
    end
    tpl = fscanf(fid,'%c');
    fclose(fid);
elseif iscellstr(tpl)
    tpl  = char(tpl);
else
    error('Invalid template.');
end

%-Split template into actions {{...}} and text
%==========================================================================
[startidx,endidx] = regexp(tpl,'\{\{[^}]*\}\}');
cur = 1;
boom = {}; % split template
while cur < numel(tpl)
    if isempty(startidx)
        % end of file
        boom{end+1} = tpl(cur:end);
        break;
    elseif startidx(1) > cur
        % text
        boom{end+1} = tpl(cur:startidx(1)-1);
        cur = startidx(1);
    else
        % action
        boom{end+1} = tpl(startidx(1):endidx(1));
        cur = endidx(1)+1;
        startidx(1) = []; endidx(1) = [];
    end
end
bt = cellfun(@(x) strncmp(x,'{{',2), boom); % is action?

%-Execute Actions
%==========================================================================

%-Initialise cursor, ie dot, '.'
%--------------------------------------------------------------------------
dot = data;

%-Loop over actions (doing it twice to make sure '-' is handled properly:
%TODO this has to be simplified)
%--------------------------------------------------------------------------
for i=find(bt)
    %-Current action
    %----------------------------------------------------------------------
    act = boom{i}(3:end-2);
    if isempty(act), continue; end
    
    %-Detect delimiters '{{- ' and ' -}}'
    %----------------------------------------------------------------------
    if strncmp(act,'- ',2)
        if i>1 && ~bt(i-1)
            boom{i-1} = deblank(boom{i-1});
        end
        act = act(2:end);
    end
    if isempty(act), continue; end
    if strncmp(fliplr(act),'- ',2)
        if i<numel(bt) && ~bt(i+1)
            boom{i+1} = fliplr(deblank(fliplr(boom{i+1})));
        end
        act = act(1:end-2);
    end
    act = strtrim(act);
    if isempty(act), continue; end
end

for i=find(bt)
    %-Current action
    %----------------------------------------------------------------------
    act = boom{i}(3:end-2);
    if isempty(act), continue; end
    
    %-Detect delimiters '{{- ' and ' -}}'
    %----------------------------------------------------------------------
    if strncmp(act,'- ',2)
        if i>1 && ~bt(i-1)
            boom{i-1} = deblank(boom{i-1});
        end
        act = act(2:end);
    end
    if isempty(act), continue; end
    if strncmp(fliplr(act),'- ',2)
        if i<numel(bt) && ~bt(i+1)
            boom{i+1} = fliplr(deblank(fliplr(boom{i+1})));
        end
        act = act(1:end-2);
    end
    act = strtrim(act);
    if isempty(act), continue; end
    
    %-Action
    %----------------------------------------------------------------------
    [t,r] = strtok(act);
    r = strtrim(r);
    boom{i} = '';
    switch lower(t)
        case '/*'
            % {{/* a comment */}}
        case 'if'
            % {{if pipeline}} T1 {{end}}
            % {{if pipeline}} T1 {{else}} T0 {{end}}
            % {{if pipeline}} T1 {{else if pipeline}} T0 {{end}}
            [el,en] = findelseend(boom,bt,i);
            tf = evaluate(r,dot);
            if isstruct(tf) || iscell(tf) || ischar(tf)
                tf = ~isempty(tf);
            end
            if ~tf
                if el
                    [boom{i+1:el}] = deal(''); % until {{else}}
                else
                    [boom{i+1:en}] = deal(''); % until {{end}}
                end
            else
                if el
                    [boom{el:en}] = deal(''); % {{else}} to {{end}}
                else
                    boom{en} = '';
                end
            end
        case 'range'
            % {{range pipeline}} T1 {{end}}
            % {{range pipeline}} T1 {{else}} T0 {{end}}
            dot_backup = dot;
            p = evaluate(r,dot);
            [el,en] = findelseend(boom,bt,i);
            if isempty(p)
                if el
                    [boom{i+1:el}] = deal('');
                else
                    [boom{i+1:en}] = deal('');
                end
            else
                if el
                    [boom{el:en}] = deal('');
                end
                out = '';
                for j=1:numel(p)
                    if isstruct(p) || isnumeric(p) || islogical(p)
                        dot = p(j);
                    elseif iscell(p)
                        dot = p{j};
                    elseif ischar(p)
                        % don't loop over chars
                    else
                        error('Cannot apply range on %s.',class(p));
                    end
                    out = [out parse_template({[boom{i+1:en}]},dot,tplist)];
                end
                [boom{i+1:en}] = deal('');
                boom{i} = out;
            end
            dot = dot_backup;
        case 'template'
            % {{template "name"}}
            % {{template "name" pipeline}}
            [n,p] = strtok(r);
            tplname = n(2:end-1);
            if ~isempty(p)
                dot = evaluate(p,dot);
            end
            if isfield(tplist,tplname)
                % fileread might not be utf-8 ready
                boom{i} = parse_template(tplist.(tplname),dot,tplist);
            else
                warning('Template %s not found.',tplname);
            end
        case 'block'
            % {{block "name" pipeline}} T1 {{end}}
            %   {{define "name"}} T1 {{end}}
            %   {{template "name" pipeline}}
            [n,p] = strtok(r);
            tplname = n(2:end-1);
            [el,en] = findelseend(boom,bt,i);
            tplist.(tplname) = {[boom{i+1:en-1}]};
            [boom{i+1:en}] = deal('');
            if ~isempty(p)
                dot = evaluate(p,dot);
            end
            boom{i} = parse_template(tplist.(tplname),dot,tplist);
        case 'with'
            % {{with pipeline}} T1 {{end}}
            % {{with pipeline}} T1 {{else}} T0 {{end}}
            if isfield(data,r(2:end)) && ~isempty(data.(r(2:end)))
                dot = data.(r(2:end));
            else
                %dot = '';
                %fprintf('ignore %s\n',r(2:end));
                [el,en] = findelseend(boom,bt,i);
                [boom{i+1:en}] = deal(''); % until {{else}} or {{end}}
            end
        case 'else'
        case 'end'
            dot = data; % change namespace
        case 'define'
            % {{define "name"}} T1 {{end}}
            tplname = r(2:end-1);
            [el,en] = findelseend(boom,bt,i);
            tplist.(tplname) = {[boom{i+1:en-1}]};
            [boom{i+1:en}] = deal('');
        otherwise
            % {{pipeline}}
            boom{i} = evaluate(act,dot);
    end
end

%-Return executed template
%--------------------------------------------------------------------------
tpl = [boom{:}];

while ~isempty(strfind(tpl,'{{'))
    % recursively execute templates as long as they contain actions
    tpl = parse_template({tpl},dot,tplist);
end


%==========================================================================
function [el,en] = findelseend(boom,bt,i)
counten = 0;
el = 0;
bt(1:i-1) = false;
for j=find(bt)
    if strcmp(strrep(strrep(boom{j},' ',''),'-',''),'{{end}}') % TODO: improve
        if ~counten
            en = j;
            return;
        else
            counten = counten - 1;
        end
    elseif strcmp(strrep(strrep(boom{j},' ',''),'-',''),'{{else}}') % TODO: improve
        if ~counten
            el = j;
        end
    %elseif any(strcmp(boom{j},{'{{ if }}','{{ range }}','{{ with }}'})) % TODO they can have a pipeline
    elseif ~isempty(regexp(boom{j},'\{\{\s*(if|range|with)[^}]*\}\}'))
        counten = counten + 1;
    end
end

%==========================================================================
function tf = evaluate(pipeline,dot)
% A pipeline is a sequence of commands separated with pipeline character '|'
%   Argument
%     boolean, chararacter array, numeric
%     '.'
%     nil
%     variable name $
%     name of a field/key of the data .Field1.Key1.Field2.Key2 or $x.Field1.Field2
%     name of a niladic method of the data .Method or $x.Method1.Field
%     name of a niladic function, eg fun
%     parenthesized instance of one the above, for grouping, eg print (.F1 arg1) (.F2 arg2)
%   .Method [Argument...]
%   functionName [Argument...]

pipeline = strsplit(pipeline,'|');
if numel(pipeline) > 1
    for i=1:numel(pipeline)
        % TODO: evaluate the chained pipeline
        tf = evaluate(pipeline,dot);
    end
else
    pipeline = pipeline{1};
end

pipeline = strtrim(pipeline);
if isequal(pipeline,'.')
    tf = dot;
elseif pipeline(1) == '.'
    var = 'dot';
    if strmatch('.Site',pipeline)
        var = 'cfg';
        pipeline = pipeline(6:end);
    end
    try
        tf = eval([var pipeline]); % sanitize...
    catch
        fprintf('Ignored %s\n',pipeline);
        tf = '';
    end
else
    % eval a function
    s = strsplit(pipeline);
    switch s{1}
        case 'eq'
            s = s(2:end);
            for i=1:numel(s)
                if ~isempty(s{i})
                    if s{i}(1) == '.'
                        try
                            s{i} = eval(['dot' s{i}]); % sanitize...
                        catch
                            s{i} = NaN; % something that it cannot match to?
                        end
                    elseif s{i}(1) == '"'
                        s{i}([1 end]) = '';
                    end
                end
            end
            tf = isequal(s{:});
        otherwise
            tf = true; % no! unknown functions
    end
end
% only if logical output is expected
% if isstruct(tf) || iscell(tf) || ischar(tf)
%     tf = isempty(tf);
% elseif isnumeric(tf) || islogical(tf)
%     tf = all(logical(tf(:)));
% end
