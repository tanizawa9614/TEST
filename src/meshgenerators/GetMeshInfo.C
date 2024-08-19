//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GetMeshInfo.h"
#include "MooseMesh.h"
#include "Conversion.h"
#include "MooseMeshUtils.h"
#include "CastUniquePointer.h"
#include "InputParameters.h"
#include "MooseTypes.h"

#include "libmesh/distributed_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/parallel_elem.h"
#include "libmesh/parallel_node.h"
#include "libmesh/compare_elems_by_level.h"
#include "libmesh/mesh_communication.h"

#include "timpi/parallel_sync.h"


#include <set>
#include <typeinfo>


#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/node.h"
// #include <iostream>
#include <fstream>


registerMooseObject("MooseApp", GetMeshInfo);

InputParameters
GetMeshInfo::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addRequiredParam<FileName>("output_file_name", "The Output file of Mesh Info");
  params.addClassDescription("Creates infomation file of mesh.");

  return params;
}


GetMeshInfo::GetMeshInfo(const InputParameters & parameters)
  : MeshGenerator(parameters), 
    _input(getMesh("input")), 
    _file_name(this->template getParam<FileName>("output_file_name"))
{
}

std::unique_ptr<MeshBase>
GetMeshInfo::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);
  // Open a file to write the output
  std::ofstream outfile(_file_name);

  // Write Node IDs
  outfile << "Node IDs:" << std::endl;
  for (const auto &elem : mesh->element_ptr_range()){
    for (const auto & node : elem->node_ref_range())
    {
        outfile << node.id() << " ";
    }
    outfile << "\n";
  }
  // Write Element IDs and connectivity
  // outfile << "\nElement IDs and Connectivity:" << std::endl;
  // for (const auto & elem : mesh.element_ptr_range())
  // {
  //     outfile << "Element ID: " << elem->id() << " Connectivity: ";
  //     for (unsigned int i = 0; i < elem->n_nodes(); ++i)
  //     {
  //         outfile << elem->node(i) << " ";
  //     }
  //     outfile << std::endl;
  // }

  // Close the output file
  outfile.close();

  return mesh;
}
