function [H,S] = structsect( H, S )
% Fields of intersection in two structures
%   [H,S] = Fields( H, S )
% Input Arguments:
%   H = Structure
%   S = Structure to seperate
% Ouput Arguments:
%   H = Structure of fields in S also in H
%   S = Structure of fields in S not in H
% Notes:
%   Intersections are not case sensitive. The intersecting structure fields are set
% to the case in H and the value in S.
%
% Example: Rename fields in structure S to match reference structure case
%    Ref = struct('FieldA',[],'FieldB',[])
%    S   = struct('fielda',1,'FIELDB',2,'FieldC',3)
%    S   = structcat(S,structsect(struct(Ref),S))
%
% See also STRUCTCAT
% ------------------------------------------------------------------------------------
if isempty(H), H = struct([]); return; end

Hnames = fieldnames(H); % Hvals = struct2cell(H);
Snames = fieldnames(S);   Svals = struct2cell(S);

[I,i,j] = intersect(upper(Snames),upper(Hnames));
if isempty(i), H = struct([]); else H = cell2struct(Svals(i),Hnames(j),1); end

[I,i] = setdiff(upper(Snames),upper(Hnames));
if isempty(i), S = struct([]); else S = cell2struct(Svals(i),Snames(i),1); end


