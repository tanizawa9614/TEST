//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowPermeabilityBase.h"
#include <string>

/**
 * Material to provide permeability taken from a variable. This material
 * is primarily designed for use with heterogeneous reservoir models
 * where the components of the permeability tensor are provided by an
 * elemental aux variables that do not change.
 * The three diagonal entries corresponding to the x, y, and z directions
 * must be given. Optionally, the off-diagonal components of the full
 * permeability tensor can be given. If they are not provided, they will be
 * initialised to zero.
 */
template <bool is_ad>
class PorousFlowPermeabilityAllBlocksTempl : public PorousFlowPermeabilityBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  PorousFlowPermeabilityAllBlocksTempl(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  const FileName& _file_name;
  const unsigned int& _col_index;
  std::vector<Real> _permeability_data; // ブロックIDとインデックスの辞書


  usingPorousFlowPermeabilityBaseMembers;
};

typedef PorousFlowPermeabilityAllBlocksTempl<false> PorousFlowPermeabilityAllBlocks;
typedef PorousFlowPermeabilityAllBlocksTempl<true> ADPorousFlowPermeabilityAllBlocks;
