lead dbs

mass_convert_config.lr = 'right';
mass_convert_config.ap = 'anterior';
mass_convert_config.ab = 'above';

uipatdir={'./BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/01_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/02_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/03_MDST_re'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/04_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/05_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/06_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/07_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/08_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/09_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/10_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/11_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/12_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/13_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/14_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/15_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/16_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/17_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/18_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/19_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/20_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/21_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/22_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/23_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/24_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/25_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/26_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/27_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/28_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/29_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/30_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/31_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/32_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/33_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/34_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/35_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/36_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/37_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/38_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/39_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/40_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/41_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/42_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/43_MDST'
    './BER_MDST_Reprocessed_VTAs_Jan2021/All_MDST/44_MDST'
    };

output_save_filepath = 'parallelized_v2_output_17.mat';
ea_mass_convert_acpc2mni(mass_convert_config, "sample_mass_convert_coordinates.xlsx", uipatdir, output_save_filepath)