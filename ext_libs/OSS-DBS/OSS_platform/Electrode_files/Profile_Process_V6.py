# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

########## function read a text line and compare and detect a given string 
def words_detect (find_word , my_string):
    w_len = len (find_word); 
    str_len = len (my_string);
    idx=0;
    str_index = []; # [word_found, indx1,index2.....indexn ]
    i=0;            # for string index word_found position
    word_found = False;
    for idx in range (0,str_len):
        if ((str_len -idx) >= w_len):
            if (find_word == my_string [idx:idx+w_len]):
                word_found = True;
                if (i == 0): str_index.append (word_found);i=1; # write the status to first of string index
                str_index.append (idx); # record string index found            
            else:  word_found = False;         
                
                        
        else:
            if( i==0):str_index.append (word_found);
            break;
    return str_index
###### test
##mystring = 'my papa say hello to my mama hello hello hello sasa '
#result = words_detect ( word, mystring )
#print result
def find_component (pattern,text):
# find the components name t=in the profile return the strings name

    myx = words_detect (pattern, text) ;
    if myx[0] : 
        test =text [myx[1]+len(pattern):]    
        name = test.split()[-2] ;
        #print name
        name = name[1:len(name)-1];
        return name
        
def find_line (pattern,line, temp_file):
    my_line = words_detect (pattern, line) ;
    if my_line[0]:        
        temp_file.write(line);
        
        
###### test 
#my_file = open('temp_file.txt','w+') #            
#text = "geompy.addToStudy( encap_body, 'encap_body' )"
#print find_component("geompy.addToStudy(",text)
#text=" geompy.addToStudyInFather( Partition_encap, Group_volume_encap, 'Group_volume_encap' )" 
#print find_component("geompy.addToStudyInFather(",text) 
#find_line("addToStudy",text,my_file) 
#my_file.close();       