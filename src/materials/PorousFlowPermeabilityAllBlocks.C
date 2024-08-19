//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPermeabilityAllBlocks.h"
#include "libmesh/elem.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <ctime>


registerMooseObject("PorousFlowApp", PorousFlowPermeabilityAllBlocks);
registerMooseObject("PorousFlowApp", ADPorousFlowPermeabilityAllBlocks);

template <bool is_ad>
InputParameters
PorousFlowPermeabilityAllBlocksTempl<is_ad>::validParams()
{
  InputParameters params = PorousFlowPermeabilityBase::validParams();
  params.addRequiredParam<FileName>("file_name", "The CSV file containing the porosity values for each block");
  params.addRequiredParam<unsigned int>("col_index", "The index of the column to read the porosity value from");
  params.addClassDescription(
      "This Material calculates the permeability tensor given by the input variables");
  return params;
}

// BOMを検出してスキップする関数
namespace PermeabilityAllBlocks{
void skipBOM(std::ifstream& inFile) {
    char c1, c2, c3;
    inFile.get(c1);
    inFile.get(c2);
    inFile.get(c3);
    if (!(c1 == char(0xEF) && c2 == char(0xBB) && c3 == char(0xBF))) {
        inFile.seekg(0); // BOMがない場合はファイルの先頭に戻る
    }
}
}

template <bool is_ad>
PorousFlowPermeabilityAllBlocksTempl<is_ad>::PorousFlowPermeabilityAllBlocksTempl(
    const InputParameters & parameters)
  : PorousFlowPermeabilityBaseTempl<is_ad>(parameters),
    _file_name(this->template getParam<FileName>("file_name")),
    _col_index(this->template getParam<unsigned int>("col_index"))
{
  std::ifstream inFile(_file_name);
  std::string line;

  if (!inFile.is_open()) {mooseError("Unable to open file: ", _file_name);}

  PermeabilityAllBlocks::skipBOM(inFile); // BOMをスキップ
  std::ofstream outFile("perm.txt", std::ios::out | std::ios::app);

  while (std::getline(inFile, line))
  {
      std::stringstream ss(line);
      std::string cell;
      std::vector<std::string> cells;

      unsigned int cnt = -1;
      while (std::getline(ss, cell, ',')){
        cnt++;
        if (cnt == _col_index) break;
      }

      if (cnt != _col_index) mooseError("Invalid CSV format or col_index out of range in file ", _file_name);

      _permeability_data.push_back(std::stod(cell));

      std::time_t now = std::time(nullptr);
      std::tm* local_time = std::localtime(&now);
      outFile << _permeability_data.size() << " " << std::put_time(local_time, "%Y-%m-%d %H:%M:%S") << "\n";
  }
}

template <bool is_ad>
void
PorousFlowPermeabilityAllBlocksTempl<is_ad>::computeQpProperties()
{
  const Elem * elem = this->_current_elem;
  SubdomainID block_id = elem->subdomain_id();

  Real permeability;
  if (block_id - 1 >= 0 || (_permeability_data.size() <= block_id)){
    permeability = _permeability_data[block_id - 1];
  }else{
    mooseError("Block ID ", block_id, " not found in CSV data");
  }
  const RealTensorValue permeability_tensor(permeability,
                                     0,
                                     0,
                                     0,
                                     permeability,
                                     0,
                                     0,
                                     0,
                                     permeability);

  _permeability_qp[_qp] = permeability_tensor;

  if (!is_ad)
  {
    (*_dpermeability_qp_dvar)[_qp].resize(_num_var, RealTensorValue());
    (*_dpermeability_qp_dgradvar)[_qp].resize(LIBMESH_DIM);

    for (const auto i : make_range(Moose::dim))
      (*_dpermeability_qp_dgradvar)[_qp][i].resize(_num_var, RealTensorValue());
  }
}

template class PorousFlowPermeabilityAllBlocksTempl<false>;
template class PorousFlowPermeabilityAllBlocksTempl<true>;
