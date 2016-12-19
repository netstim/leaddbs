% COMMON
%
% Files
%   b_img_nr               - b_img_nr.m 	reads parts of a sequence of Bruker images
%   barweb                 - Usage: handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend)
%   betacf                 - BETACF(a,b,x) is a continued fraction representation
%   betai                  - Incomplete Beta function.
%   calc_i                 - function result= calc_i(D, bTensor)
%   calculate_dti          - calculate_dti: calculate DTDstruct from DTI data file (Susanne Schnell)
%   createBmatrix          - createBmatrix
%   DE_iso                 - DE_schema12_b0_1000
%   DE_schema12_b0_1000    - DE_schema12_b0_1000 (Kamil)
%   dedirs_ge_15           - DE_schema15ge (Susanne Schnell)
%   dedirs_ge_21           - DE_schema21ge (Susanne Schnell)
%   dedirs_ge_25           - function DEdirs_GE_25 gives DE_scheme with b0s if amount of b0s is given (Susanne Schnell)
%   dedirs_ge_27           - function DEdirs_GE_27 gives DE_scheme with b0s if amount of b0s is given (Susanne Schnell)
%   dedirs_ge_45           - function DEdirs_GE_45 gives DE_scheme with b0s if amount of b0s is given (Susanne Schnell)
%   dedirs_ge_55           - function DEdirs_GE_55 gives DE_scheme with b0s if amount of b0s is given (Susanne Schnell)
%   dedirs_ge_6            - DE_schema55ge (Susanne Schnell)
%   dedirs_ge_61           - function DEdirs_GE_61 gives DE_scheme with b0s if amount of b0s is given (Susanne Schnell)
%   DicomReadSerieV2       - DicomReadSerieV2 - reads serie of dicom images and saves in mrstruct (Kamil)
%   dir2ten                - function B_tensor= dir2ten(b_val, dirVc, measureNo, fName)
%   dtdstruct_calculate    - function [res, errStr, oArg1]= dtdstruct_calculate(commandStr [, arg1[, arg2[... [, argN]]]])
%   dtdstruct_export       - function [res, errStr, oArg]= dtdstruct_export(dtdStruct, comStr, [, op1[, op2[... [, opN]]]])
%   dtdstruct_import       - function [dtdStrcut, errStr, oArg1, oArg2]= dtdstruct_import(typeStr, [, op1[, op2[... [, opN]]]])
%   dtdstruct_init         - function [res, errStr]= dtdstruct_init(typeStr, mrStruct, [mrData1, nameStr1, [mrData2, nameStr2, ...]])
%   dtdstruct_istype       - function [res, typeNameStr, errStr]dtdstruct_istype(arg)
%   dtdstruct_modify       - function [res, errStr, oArg1]= dtdstruct_modify(dtdStruct, command, parm1[, parm2[, ...[, parmN]]])
%   dtdstruct_query        - function [res, errStr]= dtdstruct_query(dtdStruct, command[, op1[, op2[... [, opN]]]])
%   dtdstruct_read         - function [res, errStr, fName]= dtdstruct_read(fName)
%   dtdstruct_write        - function [res, errStr]= dtdstruct_write(dtdStruct, fileName)
%   dw_data_admin          - function [res, errStr, res2]= dw_data_admin(commandStr, [arg1[, arg2[, ...])
%   DWI_Dicom2ADCv3        - DWI_Dicom2ADCv2   (Kamil)
%   DWI_Dicom2DTDv4        - DWI_Dicom2DTD (Kamil)
%   f_do_dti_calculation   - Kamil
%   fdist                  - FDIST( F, v1, v2) returns Q(F|v1,v2), the probability
%   fiberstruct_istype     - function [res, errStr]= fiberstruct_istype(fiberStruct)
%   fiberstruct_modify     - function [res, errStr]= fiberstruct_query(fiberStruct, command, parm1[, parm2[, ...[, parmN]]])
%   fiberstruct_read       - function [res, errStr, fName]= fiberstruct_read(fName)
%   fiberstruct_write      - function [res, errStr]= fiberstruct_write(fiberStruct, fileName)
%   files_lowercase_rec    - Converts all .m and .mat files in the current directory to lowercase filenames
%   ftrstruct_init         - function [res, errStr]= ftrstruct_init(typeStr)
%   ftrstruct_istype       - function [res, errStr]= ftrstruct_istype(dataIn, fiberType)
%   ftrstruct_modify       - function [res, errStr, res2]= ftrstruct_modify(ftrStruct, command, parm1[, parm2[, ...[, parmN]]])
%   ftrstruct_query        - function [res, errStr, res2]= ftrstruct_query(ftrStruct, commandStr[, op1[, op2[... [, opN]]]])
%   ftrstruct_read         - function [res, errStr, fName]= ftrstruct_read(fName)
%   ftrstruct_write        - function [res, errStr]= ftrstruct_write(ftrStruct, fileName)
%   gammln                 - Natural log of the complete Gamma function.
%   get_unique_str         - function [nameStr, errStr]= get_unique_str(entriesCell, defaultStr, nameGenFlag, strucFlag, forceDialogFlag)
%   hm_analyse             - function [res, errStr]= hm_analyse(inAy)
%   hm_create              - function [res, errStr]= hm_create(inAy)
%   hm_invert              - function [res, errStr]= hm_invert(M)
%   hm_mult                - function res= hm_mult(M, vec)
%   hm_rot                 - function [res, errStr]= hm_rot(alpha[rad], axis, M, flag)
%   hm_scale               - function [res, errStr]= hm_scale(scale, axis, M)
%   hm_trans               - function [res, errStr]= hm_trans(t_vect, M)
%   marching_cube          - function [verteces, triAy, errStr]= marching_cube(dataAy, verbHd)
%   maskstruct_init        - function [res, errStr]= maskstruct_init(mrStruct)
%   maskstruct_istype      - function [res, errStr]= maskstruct_istype(dataIn, maskType)
%   maskstruct_modify      - function [res, errStr, res2, res3]= maskstruct_modify(maskStruct, command, parm1[, parm2[, ...[, parmN]]])
%   maskstruct_query       - function [res, errStr, res2]= maskstruct_query(maskStruct, commandStr[, op1[, op2[... [, opN]]]])
%   maskstruct_read        - function [res, errStr, fName]= maskstruct_read(fName)
%   maskstruct_write       - function [res, errStr]= ftrstruct_write(maskStruct, fName)
%   morph_data             - function [res, errStr]= morph_data(mask, kernelAy, kernelPos, OPStr, noVoxInKernel, forceKPos)
%   morph_data_ext         - function [res, errStr]= morph_data_ext(mask, commandStr, par1, par2, par3)
%   mrstruct_export        - function [dtdStrcut, errStr, oArg1, oArg2]= mrstruct_export(typeStr, mrType, [, op1[, op2[... [, opN]]]])
%   mrstruct_import        - function [dtdStrcut, errStr, oArg1, oArg2]= mrstruct_import(typeStr, mrType, [, op1[, op2[... [, opN]]]])
%   mrstruct_transform     - function [res, errStr]= mrtsruct_transfortm(mrStruct, commandStr[, op1[, op2[... [, opN]]]])
%   my_image_zoom          - function erg= my_image_zoom(small, high, axesHd, viewStr)
%   my_mod                 - overloded function:
%   my_mrstruct_segmentate - queryStr   string which contains the name of the segmentation algorithmus
%   myboxplot              - ' myboxplot(x, y, dx, Col, Width, noWhisker)' Produces a single box plot.
%   myboxplot2             - ' myboxplot(x, y, dx, Col, Width, noWhisker)' Produces a single box plot.
%   myfactorial            - 
%   myFigPortrait          - 'h = myFigPortrait()' 
%   myFigure               - 'h = myFigure()' 
%   norm_numbers           - function [numStr, mant, expo]= norm_nummbers(num, expOffset, mantLength)
%   open_brooker           - function M= open_brooker(fName, res, gradNo, vox)
%   read_dti_dicom         - read_dti_dicom: Read all dicom DTI data indepently from manufacturer
%   ReadDicomFilesDTI_MZv5 - ReadDicomFilesDTIMZv5 for FR-Siemens sequences (Kamil)
%   RotateDeGradients      - corrects DE grads direction for the obligue slice (Kamil)
%   save_brooker           - function erg = save_brooker(fName, mrStruct)
%   spherical_harm         - function [erg, errStr, oArg1, oArg2, oArg3]= spherical_harm(commandStr, arg1, arg2, ....)
%   split_filename         - function [fName, path, errStr]= split_filename(longFileName)
%   str2filename           - function [fileStr, errStr]= str2filename(nameStr, structFlag)
%   ten2dir                - function [b_val, dirVc, b0Idx, bValIdx]= ten2dir(bTensor)
%   probstruct_init        - creates a new probStruct object
