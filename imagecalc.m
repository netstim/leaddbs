%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
%%
matlabbatch{1}.spm.util.imcalc.input = {
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/DettmerWl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/DettmerWr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/DreissigDl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/DreissigDr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/FleischmannNr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/FrommerTonil.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/FrommerTonir.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/GadMl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/GadMr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/GenschmerGerr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/HartmannDietl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/HartmannDietr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/HolwasJr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/HornHl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/HornHr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/KeilMarcol.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/KeilMarcor.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/Kerberl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/KieslichHl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/KieslichHr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/KnausOlgal.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/KnausOlgar.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/KohrtWl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/LachmannGl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/LachmannGr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/LangeWalterl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/LehmannBl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/LehmannBr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/MildeJoachiml.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/OttmarManfrel.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/OttmarManfrer.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/PrahstMirjaml.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/PrahstMirjamr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/PuchertBrigil.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/PuchertBrigir.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/RegenbergIngl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/RegenbergIngr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/ReiskiMarianl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/RoloffBirgitl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/SchabackerLul.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/SchabackerLur.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/SchatzMartinl.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/SchatzMartinr.nii,2'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/SchiekeIl.nii,1'
                                        '/Volumes/MacHD/Dropbox/MATLAB/eAuto/qm/2mm/SchiekeIr.nii,1'
                                        };
%%
matlabbatch{1}.spm.util.imcalc.output = 'mean_2mm.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16+i17+i18+i19+i20+i21+i22+i23+i24+i25+i26+i27+i28+i29+i30+i31+i32+i33+i34+i35+i36+i37+i38+i39+i40+i41+i42+i43+i44+i45)/45';
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
