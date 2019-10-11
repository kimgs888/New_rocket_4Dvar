function S = structcat( varargin )
% Concatinate structures and param,value pairs
%    S = structcat( varargin )
% Input Arguments:
%    Structures, cells of param,value pairs or param,value pairs
% Output Arguments:
%    S = Structure containing a concatination of all fields, later
%        field definitions override prior
% Notes:
%    Structure names are case insensitive. The case used is set to that
% of the last occuring field names.
% Example:
%    S = struct('ParamA',11,'ParamB',12,'ParamD',13);
%    structcat('Parama',1,'Paramb',2,'Paramc',3,S)
%
% See also STRUCTSECT
% ----------------------------------------------------------------------

S = Sexpand(varargin,struct([]));


% ------------------------------------------------------------------------------------
function S = Sexpand( V, S)
% Expand cells recursively and append structure fields
   i = 1; while i <= length(V)
      if     iscell(V{i}),    S = Sexpand(V{i},S);
      elseif isstruct(V{i}),  S = Scat(S,V{i});
      elseif ischar(V{i})
         if i+1<=length(V),   S = Scat(S,struct(V{i},V(i+1))); i=i+1;
         else error(sprintf('Param field ''%s'' without matching Value',V{i})); end
      elseif isempty(V{i});
      else error('Arguments must be of type cell,struct or char'); end
      i = i+1;
   end

% ------------------------------------------------------------------------------------
function C = Scat( A, B )
% Concatinate two structures
%   C = Scat( A, B )
% Input Arguments:
%   A = Structure
%   B = Structure to concatinate
% Ouput Arguments:
%   C = Concatinated structure of fields in A and B
% Notes:
%   Fields in B override those in A if duplicated
% ------------------------------------------------------------------------------------
   if     isempty(A), C = B;
   elseif isempty(B), C = A;
   else
      [D,E] = structsect(B,A);
      if     isempty(E), C = B;
      elseif isempty(B), C = E;
      else
         C = cell2struct( [struct2cell(E) ; struct2cell(B)],...
                          [fieldnames(E)  ; fieldnames(B) ], 1);
      end
   end

