//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPorosityMultiBlocks.h"
#include "libmesh/elem.h"
// #include "libmesh/subdomain.h"

#include <iostream>

registerMooseObject("PorousFlowApp", PorousFlowPorosityMultiBlocks);
// registerMooseObject("PorousFlowApp", ADPorousFlowPorosityMultiBlocks);

template <bool is_ad>
InputParameters
PorousFlowPorosityMultiBlocksTempl<is_ad>::validParams()
{
  InputParameters params = PorousFlowPorosityBaseTempl<is_ad>::validParams();
  params.addRequiredParam<std::vector<Real>>(
      "porosity",
      "The porosity (assumed indepenent of porepressure, temperature, "
      "strain, etc, for this material).  This should be a real number, or "
      "a constant monomial variable (not a linear lagrange or other kind of variable).");
  params.addRequiredParam<std::vector<SubdomainID>>("block_id", "Block id vector");
  // params.addRequiredCoupledVar(
  //     "porosity",
  //     "The porosity (assumed indepenent of porepressure, temperature, "
  //     "strain, etc, for this material).  This should be a real number, or "
  //     "a constant monomial variable (not a linear lagrange or other kind of variable).");
  params.addClassDescription("This Material calculates the porosity assuming it is constant");
  return params;
}

template <bool is_ad>
PorousFlowPorosityMultiBlocksTempl<is_ad>::PorousFlowPorosityMultiBlocksTempl(
    const InputParameters & parameters)
  : PorousFlowPorosityBaseTempl<is_ad>(parameters), 
    // _input_porosity(coupledValue("porosity"))
    _input_porosity(this->template getParam<std::vector<Real>>("porosity")),
    _block_ids(this->template getParam<std::vector<SubdomainID>>("block_id"))
{
  // ブロックIDとインデックスの辞書を作成
  for (size_t i = 0; i < _block_ids.size(); ++i)
  {
    _block_to_index[_block_ids[i]] = i;
  }
}

template <bool is_ad>
void
PorousFlowPorosityMultiBlocksTempl<is_ad>::initQpStatefulProperties()
{
  // ガウス点の要素を取得
  const Elem * elem = this->_current_elem;
  // const Elem * elem = this->currentElem();
  
  // 要素が属するブロックのIDを取得
  SubdomainID block_id = elem->subdomain_id();

  size_t file_idx = _block_to_index[block_id];
  _porosity[_qp] = _input_porosity[file_idx];
}

template <bool is_ad>
void
PorousFlowPorosityMultiBlocksTempl<is_ad>::computeQpProperties()
{
  initQpStatefulProperties();

  if (!is_ad)
  {
    // The derivatives are zero for all time
    (*_dporosity_dvar)[_qp].assign(_num_var, 0.0);
    (*_dporosity_dgradvar)[_qp].assign(_num_var, RealGradient());
  }
}

template class PorousFlowPorosityMultiBlocksTempl<false>;
template class PorousFlowPorosityMultiBlocksTempl<true>;
