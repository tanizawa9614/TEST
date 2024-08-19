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

/**
 * Creates lower-dimensional elements on the specified sidesets
 */
class BlocksFromSideSetsGeneratorFromFile : public MeshGenerator
{
public:
  static InputParameters validParams();

  BlocksFromSideSetsGeneratorFromFile(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;
  std::unique_ptr<MeshBase> generate2(std::unique_ptr<MeshBase> mesh, const std::vector<BoundaryName>& sideset_name, const SubdomainName& block_name);

protected:
  std::unique_ptr<MeshBase> & _input;
  /// a vector of the names of the sidesets to add the lower-D elements to
  std::vector<SubdomainID> _sideset_ids;
  std::vector<BoundaryName> _sideset_names;
  std::vector<SubdomainName> _block_names;
  const FileName& _file_name;
  const unsigned int& _col_index;

};
