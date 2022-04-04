# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 16:14:19 2018

@author: Trieu, modified by Konstantin

These dummy functions just extract indices of mesh entities from Meshes/Mesh_unref.msh
Should be substituted with a more unified approach

"""

from Electrode_files.Profile_Process_V6 import words_detect
import os

def read_mesh_indicies(dictionary):
    f1=open(os.environ['PATIENTDIR']+'/Meshes/Mesh_unref.msh','r')
    rst = 0;
    ot = 0;
    co = 0;
    en = 0;
    en_con = 0;
    roi = 0;
    flt=-1;

    list_of_contacts=[]
    for index,line in enumerate(f1):
        C1_1 = words_detect("C1_1",line)
        C1_2 = words_detect("C1_2",line)
        C1_3 = words_detect("C1_3",line)
        C1_4 = words_detect("C1_4",line)
        C1_5 = words_detect("C1_5",line)
        C1_6 = words_detect("C1_6",line)
        C1_7 = words_detect("C1_7",line)
        C1_8 = words_detect("C1_8",line)

        en_rest = words_detect("Encap_rest",line)
        en_contact = words_detect("Encap_contact",line)
        ROI = words_detect("RegOfInt",line)
        Rst = words_detect("Rst",line)
        Flt_cnt = words_detect("Flt_cnt",line)

        if (C1_1[0]):
            c1 = int(line[C1_1[1]-4:C1_1[1]-2]);
            list_of_contacts.insert(0,c1)
        if (C1_2[0]):
            c2 = int(line[C1_2[1]-4:C1_2[1]-2]);
            list_of_contacts.insert(1,c2)
        if (C1_3[0]):
            c3 = int(line[C1_3[1]-4:C1_3[1]-2]);
            list_of_contacts.insert(2,c3)
        if (C1_4[0]):
            c4 = int(line[C1_4[1]-4:C1_4[1]-2]);
            list_of_contacts.insert(3,c4)
        if (C1_5[0]):
            c5 = int(line[C1_5[1]-4:C1_5[1]-2]);
            list_of_contacts.insert(4,c5)
        if (C1_6[0]):
            c6 = int(line[C1_6[1]-4:C1_6[1]-2]);
            list_of_contacts.insert(5,c6)
        if (C1_7[0]):
            c7 = int(line[C1_7[1]-4:C1_7[1]-2]);
            list_of_contacts.insert(6,c7)
        if (C1_8[0]):
            c8 = int(line[C1_8[1]-4:C1_8[1]-2]);
            list_of_contacts.insert(7,c8)


        if (en_rest [0]):
            en = int(line[en_rest[1]-4:en_rest[1]-2]);
#            print 'Encap_Rest:{}\n'.format(en)
        if (en_contact [0]):
            en_con = int(line[en_contact[1]-4:en_contact[1]-2]);
 #           print 'Encap_Rest:{}\n'.format(en_con)
        if (ROI[0]):
            roi = int (line[ROI[1]-4:ROI[1]-2]);
 #           print 'ROI:{}\n'.format(roi)
        if (Rst[0]):
            rst = int(line[Rst[1]-4:Rst[1]-2]);
        if (Flt_cnt[0]):
            flt = int(line[Flt_cnt[1]-4:Flt_cnt[1]-2]);
 #           print 'Rst:{}\n'.format(rst)

    dictionary['Tis_ind']     = [rst,en];
    #print("Tis ind: ", dictionary['Tis_ind'])
    dictionary['ROI_ind']     = roi;
    dictionary['Contact_ind'] = en_con;
    dictionary['Rest_ind']    = rst;
    dictionary['Flt_cnt']     = flt;
    dictionary['Encup_ind']   = [en,en_con];
    dictionary['Contacts']   = list_of_contacts;
    f1.close();

    #print(dictionary)


def read_mesh_indicies_extended(dictionary):
    f1=open(os.environ['PATIENTDIR']+'/Meshes/Mesh_unref.msh','r')
    rst = 0;
    ot = 0;
    co = 0;
    en = 0;
    en_con = 0;
    roi = 0;
    flt=-1;

    list_of_contacts=[]                     #this contains indices of contacts in gmsh
    list_of_floating=[]                     #this contains indices of floatin in gmsh
    N_contacts_on_lead=[]                    #contains numbers of the contacts on the lead that are active   (starting from 1)
    N_floats_on_lead=[]                    #contains numbers of the floats on the lead  (starting from 1)
    for index,line in enumerate(f1):
        C1_1 = words_detect("C1_1",line)
        C1_2 = words_detect("C1_2",line)
        C1_3 = words_detect("C1_3",line)
        C1_4 = words_detect("C1_4",line)
        C1_5 = words_detect("C1_5",line)
        C1_6 = words_detect("C1_6",line)
        C1_7 = words_detect("C1_7",line)
        C1_8 = words_detect("C1_8",line)
        en_rest = words_detect("Encap_rest",line)
        en_contact = words_detect("Encap_contact",line)
        ROI = words_detect("RegOfInt",line)
        Rst = words_detect("Rst",line)
        Flt_cnt1 = words_detect("Flt_cnt1",line)
        Flt_cnt2 = words_detect("Flt_cnt2",line)
        Flt_cnt3 = words_detect("Flt_cnt3",line)
        Flt_cnt4 = words_detect("Flt_cnt4",line)
        Flt_cnt5 = words_detect("Flt_cnt5",line)
        Flt_cnt6 = words_detect("Flt_cnt6",line)
        Flt_cnt7 = words_detect("Flt_cnt7",line)
        Flt_cnt8 = words_detect("Flt_cnt8",line)

        if (C1_1[0]):
            c1 = int(line[C1_1[1]-4:C1_1[1]-2]);
            list_of_contacts.insert(0,c1)
            N_contacts_on_lead.append(1)
        if (C1_2[0]):
            c2 = int(line[C1_2[1]-4:C1_2[1]-2]);
            list_of_contacts.insert(1,c2)
            N_contacts_on_lead.append(2)
        if (C1_3[0]):
            c3 = int(line[C1_3[1]-4:C1_3[1]-2]);
            list_of_contacts.insert(2,c3)
            N_contacts_on_lead.append(3)
        if (C1_4[0]):
            c4 = int(line[C1_4[1]-4:C1_4[1]-2]);
            list_of_contacts.insert(3,c4)
            N_contacts_on_lead.append(4)
        if (C1_5[0]):
            c5 = int(line[C1_5[1]-4:C1_5[1]-2]);
            list_of_contacts.insert(4,c5)
            N_contacts_on_lead.append(5)
        if (C1_6[0]):
            c6 = int(line[C1_6[1]-4:C1_6[1]-2]);
            list_of_contacts.insert(5,c6)
            N_contacts_on_lead.append(6)
        if (C1_7[0]):
            c7 = int(line[C1_7[1]-4:C1_7[1]-2]);
            list_of_contacts.insert(6,c7)
            N_contacts_on_lead.append(7)
        if (C1_8[0]):
            c8 = int(line[C1_8[1]-4:C1_8[1]-2]);
            list_of_contacts.insert(7,c8)
            N_contacts_on_lead.append(8)
        if (en_rest [0]):
            en = int(line[en_rest[1]-4:en_rest[1]-2]);
#            print 'Encap_Rest:{}\n'.format(en)
        if (en_contact [0]):
            en_con = int(line[en_contact[1]-4:en_contact[1]-2]);
 #           print 'Encap_Rest:{}\n'.format(en_con)
        if (ROI[0]):
            roi = int (line[ROI[1]-4:ROI[1]-2]);
 #           print 'ROI:{}\n'.format(roi)
        if (Rst[0]):
            rst = int(line[Rst[1]-4:Rst[1]-2]);

        if (Flt_cnt1[0]):
            flt1 = int(line[Flt_cnt1[1]-4:Flt_cnt1[1]-2]);
            list_of_floating.insert(0,flt1)
            N_floats_on_lead.append(1)
        if (Flt_cnt2[0]):
            flt2 = int(line[Flt_cnt2[1]-4:Flt_cnt2[1]-2]);
            list_of_floating.insert(1,flt2)
            N_floats_on_lead.append(2)
        if (Flt_cnt3[0]):
            flt3 = int(line[Flt_cnt3[1]-4:Flt_cnt3[1]-2]);
            list_of_floating.insert(2,flt3)
            N_floats_on_lead.append(3)
        if (Flt_cnt4[0]):
            flt4 = int(line[Flt_cnt4[1]-4:Flt_cnt4[1]-2]);
            list_of_floating.insert(3,flt4)
            N_floats_on_lead.append(4)
        if (Flt_cnt5[0]):
            flt5 = int(line[Flt_cnt5[1]-4:Flt_cnt5[1]-2]);
            list_of_floating.insert(4,flt5)
            N_floats_on_lead.append(5)
        if (Flt_cnt6[0]):
            flt6 = int(line[Flt_cnt6[1]-4:Flt_cnt6[1]-2]);
            list_of_floating.insert(5,flt6)
            N_floats_on_lead.append(6)
        if (Flt_cnt7[0]):
            flt7 = int(line[Flt_cnt7[1]-4:Flt_cnt7[1]-2]);
            list_of_floating.insert(6,flt7)
            N_floats_on_lead.append(7)
        if (Flt_cnt8[0]):
            flt8 = int(line[Flt_cnt8[1]-4:Flt_cnt8[1]-2]);
            list_of_floating.insert(7,flt8)
            N_floats_on_lead.append(8)



 #           print 'Rst:{}\n'.format(rst)

    dictionary['Tis_ind']     = [rst,en];
    #print("Tis ind: ", dictionary['Tis_ind'])
    dictionary['ROI_ind']     = roi;
    dictionary['Contact_ind'] = en_con;
    dictionary['Rest_ind']    = rst;
    dictionary['Flt_contacts']     = list_of_floating;
    dictionary['Encup_ind']   = [en,en_con];
    dictionary['Contacts']   = list_of_contacts;
    dictionary['Active_contacts_on_lead']   = N_contacts_on_lead;
    dictionary['Float_contacts_on_lead']   = N_floats_on_lead;
    f1.close();
