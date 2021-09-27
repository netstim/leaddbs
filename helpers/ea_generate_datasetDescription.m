function ea_generate_datasetDescription(dest_filepath)
    [parent_dir,parent_file,ext] = fileparts(dest_filepath);
    dataset_description.Name = parent_file;
    dataset_description.BIDSVersion = '1.6.0';
    dataset_description.DatasetType = 'raw'; %for backwards compatibility, as suggested by BIDS
    output_file = fullfile(parent_dir,parent_file,'dataset_description.json');
    json_fid = fopen(output_file,'w');
    encodejson = jsonencode(dataset_description);
    fprintf(json_fid,encodejson);
    output_file_ignore = fullfile(parent_dir,parent_file,'.bidsignore');
    opfile_ignore_fid = fopen(output_file_ignore,'w');
  
    fprintf(opfile_ignore_fid,'/derivatives\n');
    fprintf(opfile_ignore_fid,'/rawdata/**/**/**/*_CT.nii.gz');
    fprintf(opfile_ignore_fid,'*.xslx');
end