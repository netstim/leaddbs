#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 11:03:09 2020

@author: butenko
"""

# this sript prerefines cells where the external grounding will be used
from dolfin import *
import numpy as np

def simple_mesh_refiner(mesh_old,boundaries,subdomains_assigned,cell_markers):
    parameters['linear_algebra_backend']='PETSc'
    parameters["refinement_algorithm"] = "plaza_with_parent_facets"
    parameters["allow_extrapolation"] = True;
    
    # refine marked cells on the old mesh and adapt other entities
    mesh_new = refine(mesh_old, cell_markers)
    subdomains_assigned_new=adapt(subdomains_assigned,mesh_new)
    boundaries_new = adapt(boundaries,mesh_new) # put function space

    return (mesh_new,boundaries_new,subdomains_assigned_new)

def prerefine_ground(mesh,subdomains_assigned,boundaries,Domains):
    
    tdim = mesh.topology().dim()
    mesh.init(tdim-1, tdim)

    zmin = mesh.coordinates()[:, 2].min()   #assuming that z is dorso-ventral axis
    #print("zmin: ",zmin)
    ground_height=1000.0
  
    cells_to_ref= MeshFunction('bool',mesh,3)
    cells_to_ref.set_all(False)

    subdomains_just_here= MeshFunction('size_t',mesh,3)
    subdomains_just_here.set_all(0)

    if Domains.Float_contacts!=-1:
        if isinstance(Domains.Float_contacts,int):
            subdomains_just_here.array()[subdomains_assigned.array()==Domains.Float_contacts]=5
        else:
            for i in range(len(Domains.Float_contacts)):
                subdomains_just_here.array()[subdomains_assigned.array()==Domains.Float_contacts[i]]=5          #5 is index for float contact (very high cond and perm)            
         
    for i in range(len(Domains.Encup_index)):
        subdomains_just_here.array()[subdomains_assigned.array()==Domains.Encup_index[i]]=4          #4 is index of encap    

    facets_bc=MeshFunction('size_t',mesh, 2)
    facets_bc.set_all(0)

    for cell in cells(mesh): 
        z_coord=cell.midpoint().z()
        if z_coord<zmin+ground_height and subdomains_just_here[cell]!=4 and subdomains_just_here[cell]!=5:
            for facet in facets(cell):
                if facet.exterior():
                    facets_bc[facet]=1
                    cells_to_ref[cell]=True
                    break
                    #facets_bc[facet] = 1    

    dss=Measure("ds",domain=mesh,subdomain_data=facets_bc)  
    Ground_surface=assemble(1.0*dss(1))
    print("Ground_surface: ",Ground_surface)

    file=File('/opt/Patient/Results_adaptive/grounding_cells.pvd')
    file<<cells_to_ref  
                    
    mesh_refined,boundaries_refined,subdomains_assigned_refined=simple_mesh_refiner(mesh,boundaries,subdomains_assigned,cells_to_ref)
    
    return mesh_refined,subdomains_assigned_refined,boundaries_refined,Ground_surface

def refine_external_ground(Domains):        #if we want to skip adaptive mesh refinement
 
    mesh = Mesh("/opt/Patient/Meshes/Mesh_unref.xml")
    boundaries = MeshFunction('size_t',mesh,'/opt/Patient/Meshes/Mesh_unref_facet_region.xml')
    subdomains_assigned=MeshFunction('size_t',mesh,"/opt/Patient/Meshes/Mesh_unref_physical_region.xml")

    print("Before ground ref: ",mesh.num_cells())
    
    # grounding_surface_div=100.0   #100%  
    #mesh,subdomains_assigned,boundaries,ground_surface_old=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    #mesh,subdomains_assigned,boundaries,ground_surface_old=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    #mesh,subdomains_assigned,boundaries,ground_surface_old=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    #mesh,subdomains_assigned,boundaries,ground_surface_old=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    # #only after two prerefinement iterations the ground is properly placed, but need to check further
    # while grounding_surface_div>5.0:
    
    #     mesh,subdomains_assigned,boundaries,ground_surface=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    #     grounding_surface_div=abs(ground_surface-ground_surface_old)*100.0/ground_surface
    #     ground_surface_old=ground_surface

    print("After ground ref: ",mesh.num_cells())
        
    mesh_file=File('/opt/Patient/Results_adaptive/mesh_adapt.xml.gz')
    boundaries_file = File('/opt/Patient/Results_adaptive/boundaries_adapt.xml')
    subdomains_assigned_file=File('/opt/Patient/Results_adaptive/subdomains_assigned_adapt.xml')
    
    mesh_file<<mesh
    boundaries_file<<boundaries
    subdomains_assigned_file<<subdomains_assigned    

    return True
