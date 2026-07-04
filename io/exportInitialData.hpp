#ifndef EXPORT_INITIAL_DATA_HPP
#define EXPORT_INITIAL_DATA_HPP

#include "../Declarations.hpp"
#include "cgnslib.h"
#include <iostream>

inline void ExportInitialData(const std::string& filename,const Mesh& mesh,const std::vector<ConservedVars>& U)
{
    int file_id=0;
    if(cg_open(filename.c_str(),CG_MODE_WRITE,&file_id)!=CG_OK)return;

    int base_id=1;
    cg_base_write(file_id,"Base",3,3,&base_id);

    cgsize_t zone_sizes[3]={static_cast<cgsize_t>(mesh.nodes.size()),static_cast<cgsize_t>(mesh.cells.size()),0};

    int zone_id=1;
    cg_zone_write(file_id,base_id,"Zone",zone_sizes,CGNS_ENUMV(Unstructured),&zone_id);

    std::vector<double>x(mesh.nodes.size()),y(mesh.nodes.size()),z(mesh.nodes.size());

    for(size_t i=0;i<mesh.nodes.size();i++)
    {
        x[i]=mesh.nodes[i].coord.x;
        y[i]=mesh.nodes[i].coord.y;
        z[i]=mesh.nodes[i].coord.z;
    }

    int coord_id=0;
    cg_coord_write(file_id,base_id,zone_id,RealDouble,"CoordinateX",x.data(),&coord_id);
    cg_coord_write(file_id,base_id,zone_id,RealDouble,"CoordinateY",y.data(),&coord_id);
    cg_coord_write(file_id,base_id,zone_id,RealDouble,"CoordinateZ",z.data(),&coord_id);


    int flow_id=0;
    cg_sol_write(file_id,base_id,zone_id, "InitialField", CellCenter, &flow_id);

    std::vector<double>field(mesh.cells.size());
    int field_id=0;

    for(size_t i=0;i<mesh.cells.size();i++)field[i]=U[i].rho;
    cg_field_write(file_id,base_id,zone_id,flow_id,RealDouble,"Density",field.data(),&field_id);

    for(size_t i=0;i<mesh.cells.size();i++)field[i]=U[i].rhou;
    cg_field_write(file_id,base_id,zone_id,flow_id,RealDouble,"MomentumX",field.data(),&field_id);

    for(size_t i=0;i<mesh.cells.size();i++)field[i]=U[i].rhov;
    cg_field_write(file_id,base_id,zone_id,flow_id,RealDouble,"MomentumY",field.data(),&field_id);

    for(size_t i=0;i<mesh.cells.size();i++)field[i]=U[i].rhow;
    cg_field_write(file_id,base_id,zone_id,flow_id,RealDouble,"MomentumZ",field.data(),&field_id);

    for(size_t i=0;i<mesh.cells.size();i++)field[i]=U[i].E;
    cg_field_write(file_id,base_id,zone_id,flow_id,RealDouble,"Energy",field.data(),&field_id);


    int wall=0,inlet=0,outlet=0,farfield=0,symmetry=0,unknown=0;

    for(const auto& face:mesh.faces)
    {
        if(face.right_cell==-1)
        {
            if(face.boco_marker==1)inlet++;
            else if(face.boco_marker==2)wall++;
            else if(face.boco_marker==3)farfield++;
            else if(face.boco_marker==4)outlet++;
            else if(face.boco_marker==5)symmetry++;
            else unknown++;
        }
    }


    cg_close(file_id);

    std::cout<<"\nBoundary Summary\n";
    std::cout<<"----------------------------\n";
    std::cout<<"Wall        : "<<wall<<" faces\n";
    std::cout<<"Inlet       : "<<inlet<<" faces\n";
    std::cout<<"Outlet      : "<<outlet<<" faces\n";
    std::cout<<"Farfield    : "<<farfield<<" faces\n";
    std::cout<<"Symmetry    : "<<symmetry<<" faces\n";

    if(unknown)std::cout<<"Unknown     : "<<unknown<<" faces\n";
    std::cout<<"----------------------------\n";
    std::cout<<"Initial CGNS written : "<<filename<<"\n";
}

#endif