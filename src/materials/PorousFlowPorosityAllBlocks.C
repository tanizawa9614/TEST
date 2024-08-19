#include "PorousFlowPorosityAllBlocks.h"
#include "libmesh/elem.h"
#include <fstream>
#include <sstream>
#include <iostream>

registerMooseObject("PorousFlowApp", PorousFlowPorosityAllBlocks);

template <bool is_ad>
InputParameters
PorousFlowPorosityAllBlocksTempl<is_ad>::validParams()
{
  InputParameters params = PorousFlowPorosityBaseTempl<is_ad>::validParams();
  params.addRequiredParam<FileName>("file_name", "The CSV file containing the porosity values for each block");
  params.addRequiredParam<unsigned int>("col_index", "The index of the column to read the porosity value from");
  params.addClassDescription("This Material calculates the porosity assuming it is constant");
  return params;
}

// BOMを検出してスキップする関数
namespace PorosityAllBlocks{
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
PorousFlowPorosityAllBlocksTempl<is_ad>::PorousFlowPorosityAllBlocksTempl(
    const InputParameters & parameters)
  : PorousFlowPorosityBaseTempl<is_ad>(parameters), 
    _file_name(this->template getParam<FileName>("file_name")),
    _col_index(this->template getParam<unsigned int>("col_index"))
{
  std::ifstream inFile(_file_name);
  std::string line;

  if (!inFile.is_open()) mooseError("Unable to open file: ", _file_name);

  PorosityAllBlocks::skipBOM(inFile); // BOMをスキップ

  while (std::getline(inFile, line)){
      std::stringstream ss(line);
      std::string cell;
      std::vector<std::string> cells;

      unsigned int cnt = -1;
      while (std::getline(ss, cell, ',')){
        cnt++;
        if (cnt == _col_index) break;
      }

      if (cnt != _col_index) mooseError("Invalid CSV format or col_index out of range in file ", _file_name);

      _porosity_data.push_back(std::stod(cell));
  }
}

template <bool is_ad>
void
PorousFlowPorosityAllBlocksTempl<is_ad>::initQpStatefulProperties()
{
  const Elem * elem = this->_current_elem;
  SubdomainID block_id = elem->subdomain_id();

  if (block_id - 1 >= 0 || (_porosity_data.size() <= block_id))
  {
    this->_porosity[_qp] = _porosity_data[block_id-1];
  }
  else
  {
    mooseError("Block ID ", block_id, " not found in CSV data");
  }
}

template <bool is_ad>
void
PorousFlowPorosityAllBlocksTempl<is_ad>::computeQpProperties()
{
  initQpStatefulProperties();

  if (!is_ad)
  {
    // The derivatives are zero for all time
    (*_dporosity_dvar)[_qp].assign(this->_num_var, 0.0);
    (*_dporosity_dgradvar)[_qp].assign(this->_num_var, RealGradient());
  }
}

template class PorousFlowPorosityAllBlocksTempl<false>;
template class PorousFlowPorosityAllBlocksTempl<true>;
