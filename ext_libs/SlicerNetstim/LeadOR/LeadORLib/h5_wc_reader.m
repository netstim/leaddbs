classdef h5_wc_reader < handle
	properties
        sr
        min_index
        max_segments
        raw_filename
        sr_infile
        t0_segments
        max_index
        segmentLength
        spikes_file
    end 
	methods 
        function obj = h5_wc_reader(par, raw_filename)
            obj.sr = [];
            obj.max_segments = [];
            finfo = h5info(raw_filename);
            obj.raw_filename = raw_filename;
            obj.spikes_file = false;
            
            if ismember('sr',{finfo.Datasets.Name})  %if is possible, load sr from file; else from set_parameters
                sr = h5read(raw_filename,'/sr'); 
                obj.sr_infile = true;
                obj.sr = double(sr);
            else
                obj.sr_infile = false;
                if isfield(par,'sr')
                    obj.sr = par.sr;
                end
            end
            
            if ismember('spikes',{finfo.Datasets.Name}) && ismember('index',{finfo.Datasets.Name})
               obj.spikes_file = true;
            end
            if isfield(par,'tmax') && ismember('data',{finfo.Datasets.Name})
                
                data_info =  h5info(raw_filename,'/data');
                obj.min_index = floor(par.tmin * obj.sr);
                if strcmp(par.tmax,'all')
                    obj.max_index = max(data_info.Dataspace.Size);
                else
                    obj.max_index = ceil(par.tmax * obj.sr);
                end
                n = obj.max_index -  obj.min_index;

                %Segments the data in par.segments pieces
                obj.max_segments = ceil(n/ obj.sr / ...
                        (par.segments_length * 60));         %number of segments in which data is cutted


                obj.segmentLength = floor (n/obj.max_segments);

                obj.t0_segments = ones(1,obj.max_segments);
                obj.t0_segments(1) = (obj.min_index-1)/obj.sr*1000;
                for i = 2:obj.max_segments
                    obj.t0_segments(i) = obj.t0_segments(i-1) + obj.segmentLength/obj.sr*1000;
                end
            end
                
                        
        end

        function [spikes, index] = load_spikes(obj)
             spikes = h5read(obj.raw_filename, '/spikes'); 
             index = h5read(obj.raw_filename, '/index'); 
        end
        
        function [sr,max_segments,with_raw,with_spikes] = get_info(obj)
            if isempty(obj.max_segments)
                with_raw = false;
                max_segments = 0;
            else
                with_raw = true;
                max_segments = obj.max_segments;
            end
            
            with_spikes = obj.spikes_file;
            if obj.sr_infile
                sr = obj.sr;
            else
                sr = [];
            end
        end
        
        function index_ts = index2ts(obj,index,i)
            index_ts = (index)/obj.sr*1000 + obj.t0_segments(i);
        end
      
        function x = get_segment(obj,i)
            
            data = h5read(obj.raw_filename,'/data');
            
            if i ~= obj.max_segments
                x = data(obj.min_index+obj.segmentLength*(i-1)+1 : obj.min_index+obj.segmentLength*i);
            else
                x = data(obj.min_index+obj.segmentLength*(i-1)+1 : obj.max_index);
            end
        end
    end
end
