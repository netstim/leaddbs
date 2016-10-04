function score=ea_exportscore()

% for now static sample export:

score.studyOrigin='Earlystim'; % Tag label study name
score.studyDoi='10.1056/NEJMoa1205158'; % Tag label study name
score.metric='UPDRS-III'; % clinical score label
score.dbsOff.vals=[1,2,0,0,0,0,0,0,0,0,0,0,0,0,2,3,2,2,3,3,1,1,1,2,2,2,2];
score.dbsOff.mean=(mean(score.dbsOff.vals));
score.dbsOn.vals=[0,2,0,0,0,0,0,2,2,0,0,0,0,0,3,3,1,2,1,3,0,2,1,2,2,2,2];
score.dbsOn.mean=(mean(score.dbsOn.vals));
score.patient.DOB='05-10-1951'; % date of birth
score.patient.DOS='03-02-2014'; % date of surgery
score.patient.disease='Parkinson''''s Disease';
score.patient.diseaseSubtype='Tremor Dominant';
score.version=1.0;

ea_savejson('',score,ea_setjsopts('sample.json'));
















function opt=ea_setjsopts(fname)


   opt.FileName=fname;
%        opt.FloatFormat ['%.10g'|string]: format to show each numeric element
%                         of a 1D/2D array;
%        opt.ArrayIndent [1|0]: if 1, output explicit data array with
%                         precedent indentation; if 0, no indentation
%        opt.ArrayToStruct[0|1]: when set to 0, savejson outputs 1D/2D
%                         array in JSON array format; if sets to 1, an
%                         array will be shown as a struct with fields
%                         "_ArrayType_", "_ArraySize_" and "_ArrayData_"; for
%                         sparse arrays, the non-zero elements will be
%                         saved to _ArrayData_ field in triplet-format i.e.
%                         (ix,iy,val) and "_ArrayIsSparse_" will be added
%                         with a value of 1; for a complex array, the 
%                         _ArrayData_ array will include two columns 
%                         (4 for sparse) to record the real and imaginary 
%                         parts, and also "_ArrayIsComplex_":1 is added. 
%        opt.ParseLogical [0|1]: if this is set to 1, logical array elem
%                         will use true/false rather than 1/0.
%        opt.NoRowBracket [1|0]: if this is set to 1, arrays with a single
%                         numerical element will be shown without a square
%                         bracket, unless it is the root object; if 0, square
%                         brackets are forced for any numerical arrays.
%        opt.ForceRootName [0|1]: when set to 1 and rootname is empty, savejson
%                         will use the name of the passed obj variable as the 
%                         root object name; if obj is an expression and 
%                         does not have a name, 'root' will be used; if this 
%                         is set to 0 and rootname is empty, the root level 
%                         will be merged down to the lower level.
%        opt.Inf ['"$1_Inf_"'|string]: a customized regular expression pattern
%                         to represent +/-Inf. The matched pattern is '([-+]*)Inf'
%                         and $1 represents the sign. For those who want to use
%                         1e999 to represent Inf, they can set opt.Inf to '$11e999'
%        opt.NaN ['"_NaN_"'|string]: a customized regular expression pattern
%                         to represent NaN
%        opt.JSONP [''|string]: to generate a JSONP output (JSON with padding),
%                         for example, if opt.JSON='foo', the JSON data is
%                         wrapped inside a function call as 'foo(...);'
%        opt.UnpackHex [1|0]: conver the 0x[hex code] output by loadjson 
%                         back to the string form