//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowPorosityBase.h" 

/**
 * Material to provide a constant value of porosity. This can be specified
 * by either a constant value in the input file, or taken from an aux variable.
 * Note: this material assumes that the porosity remains constant throughout a
 * simulation, so the coupled aux variable porosity must also remain constant.
 */
template <bool is_ad>
class PorousFlowPorosityMultiBlocksTempl : public PorousFlowPorosityBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  PorousFlowPorosityMultiBlocksTempl(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Constant porosity (Real constant Monomial variable only so no AD version)
  // const VariableValue & _input_porosity;
  const std::vector<Real> & _input_porosity;
  // const Function & _func;

  std::vector<SubdomainID> _block_ids;
  std::map<SubdomainID, size_t> _block_to_index; // ブロックIDとインデックスの辞書

  usingPorousFlowPorosityBaseMembers;
};

typedef PorousFlowPorosityMultiBlocksTempl<false> PorousFlowPorosityMultiBlocks;
typedef PorousFlowPorosityMultiBlocksTempl<true> ADPorousFlowPorosityMultiBlocks;
