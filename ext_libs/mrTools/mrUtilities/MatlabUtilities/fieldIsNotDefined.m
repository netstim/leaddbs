function notDefined = fieldIsNotDefined(structure,fieldName)
%
% notDefined = fieldIsNotDefined(struct,fieldName)
%
%      Author: Julien Besle
%        Date: 18/12/2010
%         $Id$
%     Purpose: Determine if a field is defined in structure structure
%               A field is defined if (a) it exists and (b) it is not empty.
%
%
%  This routine should replace the many calls of the form
%    if ~isfield(structure,'fieldname') || isempty(structure.fieldname)
%  with the call
%    if fieldIsNotDefined(structure,'fieldname')

if (~ischar(fieldName)), error('Field name must be a string'); end
if (~isstruct(structure)), error('First argument must be a structure'); end

notDefined = ~isfield(structure,fieldName) || isempty(structure.(fieldName));
