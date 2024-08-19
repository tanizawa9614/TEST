//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshGenerator.h"

/*
 * Mesh generator to create a new node set and a new boundary with the nodes the user provides
 */
class BlockFromNodesGeneratorFromFile : public MeshGenerator
{
public:
  static InputParameters validParams();

  BlockFromNodesGeneratorFromFile(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

  std::unique_ptr<MeshBase> 
  GenerateNodeSet(
          std::unique_ptr<MeshBase> mesh,
          const std::vector<unsigned int> &node_ids,
          const std::vector<BoundaryName> &boundary_name);
  
  std::unique_ptr<MeshBase> 
  GenerateSideSet(
            std::unique_ptr<MeshBase> mesh);
  
  std::unique_ptr<MeshBase>
  GenerateBlock(
          std::unique_ptr<MeshBase> mesh,
          const std::vector<BoundaryName> &boundary_name,
          const SubdomainName &block_name);

protected:
  /// mesh to modify
  std::unique_ptr<MeshBase> & _input;
  const FileName & _file_name;
  std::vector<std::string> _block_names;
  std::vector<unsigned int> _node_ids;
  std::vector<unsigned int> _node_offsets;
};
