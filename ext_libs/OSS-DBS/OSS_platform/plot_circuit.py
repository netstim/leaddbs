from GUI_BG_Th_Cort_circuit import drawResults
import pickle

"""
    Usage:
        drawResults(results, electrodeLocation)

    The order or case of the names do not matter.
"""

dict_connect_short = {      
'HDP_MC_STN'                    :  0.0 ,
'HDP_STN_GPi'                    :  0.0 ,
'HDP_STN_SN'                    :  0.0 ,

'Direct_MC_Str'                    :  0.0 ,
'Direct_Str_GPi'                    :  0.0 ,
'Direct_Str_SN'                    :  0.0 ,

'Indirect_MC_Str'                    :  0.0 ,
'Indirect_Str_GPe'                    :  0.0 ,
'Indirect_GPe_GPi'                    :  0.0,
'Indirect_GPe_STN'                    :  0.0 ,
'Indirect_GPe_SN'                    :  0.0 ,
'Indirect_STN_GPi'                    :  0.0 ,
'Indirect_STN_SN'                    :  0.0 ,

'Excitatory_STN_GPe'                  :  0.0 ,
'Inhibitory_GPi_Th'                  :  0.0 ,
'Inhibitory_SN_Th'                  :  0.0 ,
}

with open('connections_status.pkl', "rb") as f:
            dict_connect = pickle.load(f)
            
#print(dict_connect)
#scale the connections with multiple branches
i=0
for key in dict_connect:
    if dict_connect[key][0]!=0:
        dict_connect_short[key]=dict_connect[key][0]/dict_connect[key][1]
        #print(dict_connect_short[key])
    i+=1

print(dict_connect_short)

# Example: electrode is located at the thalamus. 'th' or 'thalamus' could used. this is also case insensitive
drawResults(dict_connect_short, 'STN')
