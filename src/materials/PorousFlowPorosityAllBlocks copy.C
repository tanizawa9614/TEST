//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPorosityAllBlocks.h"
#include "libmesh/elem.h"
// #include "libmesh/subdomain.h"

#include <iostream>

registerMooseObject("PorousFlowApp", PorousFlowPorosityAllBlocks);
// registerMooseObject("PorousFlowApp", ADPorousFlowPorosityAllBlocks);

template <bool is_ad>
InputParameters
PorousFlowPorosityAllBlocksTempl<is_ad>::validParams()
{
  InputParameters params = PorousFlowPorosityBaseTempl<is_ad>::validParams();
  params.addRequiredParam<std::string>("file_name");
  params.addRequiredParam<unsigned int>("col_index");
  params.addClassDescription("This Material calculates the porosity assuming it is constant");
  return params;
}

template <bool is_ad>
PorousFlowPorosityAllBlocksTempl<is_ad>::PorousFlowPorosityAllBlocksTempl(
    const InputParameters & parameters)
  : PorousFlowPorosityBaseTempl<is_ad>(parameters), 
    _file_name(this->template getParam<std::string>("file_name")),
    _col_index(this->template getParam<unsigned int>("col_index")),
{
  std::ifstream inFile(_file_name);
  std::string line;
  unsigned int row_num;

  while (std::getline(file, line))
  {
      std::stringstream ss(line);
      SubdomainID subdomain_id;
      Real porosity;
      char comma;

      ss >> subdomain_id >> comma >> porosity;
      if (ss.fail() || comma != ',')
        mooseError("Invalid CSV format in file ", filename);
      _porosity_data[subdomain_id] = porosity;
  }
}

template <bool is_ad>
void
PorousFlowPorosityAllBlocksTempl<is_ad>::initQpStatefulProperties()
{
  // ガウス点の要素を取得
  const Elem * elem = this->_current_elem;
  // const Elem * elem = this->currentElem();
  
  // 要素が属するブロックのIDを取得
  SubdomainID block_id = elem->subdomain_id();

  _porosity[_qp] = _porosity_data[block_id];
}

template <bool is_ad>
void
PorousFlowPorosityAllBlocksTempl<is_ad>::computeQpProperties()
{
  initQpStatefulProperties();

  if (!is_ad)
  {
    // The derivatives are zero for all time
    (*_dporosity_dvar)[_qp].assign(_num_var, 0.0);
    (*_dporosity_dgradvar)[_qp].assign(_num_var, RealGradient());
  }
}

template class PorousFlowPorosityAllBlocksTempl<false>;
template class PorousFlowPorosityAllBlocksTempl<true>;
