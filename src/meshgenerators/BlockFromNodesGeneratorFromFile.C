//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BlockFromNodesGeneratorFromFile.h"
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


registerMooseObject("MooseApp", BlockFromNodesGeneratorFromFile);

InputParameters
BlockFromNodesGeneratorFromFile::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addRequiredParam<FileName>("file_name", "The CSV file containing the porosity values for each block");
  params.addRequiredParam<unsigned int>("block_name_column_index", "The index of the column to read the porosity value from");
  params.addRequiredParam<unsigned int>("node_offset_column_index", "The index of the column to read the porosity value from");
  params.addClassDescription("Creates new node sets, new boundaries and new blocks made with the nodes the user provides.");

  return params;
}

// BOMを検出してスキップする関数
namespace BlockFromNodesGeneratorFromFileBOM{
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

// Used to temporarily store information about which lower-dimensional
// sides to add and what subdomain id to use for the added sides.
namespace BlockFromNodesFromFile{
struct ElemSideDouble
{
  ElemSideDouble(Elem * elem_in, unsigned short int side_in) : elem(elem_in), side(side_in) {}

  Elem * elem;
  unsigned short int side;
};
}

BlockFromNodesGeneratorFromFile::BlockFromNodesGeneratorFromFile(const InputParameters & parameters)
  : MeshGenerator(parameters), 
    _input(getMesh("input")), 
    _file_name(this->template getParam<FileName>("file_name"))
{
  const unsigned int block_col_index= 
              this->template getParam<unsigned int>("block_name_column_index");
  const unsigned int node_offset_col_index = 
              this->template getParam<unsigned int>("node_offset_column_index");

  std::ifstream inFile(_file_name);
  std::string line;

  if (!inFile.is_open()) mooseError("Unable to open file: ", _file_name);

  BlockFromNodesGeneratorFromFileBOM::skipBOM(inFile); // BOMをスキップ

  while (std::getline(inFile, line)){
    std::stringstream ss(line);
    std::string cell;

    unsigned int cnt = -1;
    while (std::getline(ss, cell, ',')){
      cnt++;
      cell.erase(cell.find_last_not_of(" \n\r\t") + 1);
      if (cnt == 1 && cell == "0") break;
      if (cnt == block_col_index)
        _block_names.push_back(cell);
      if (cnt == node_offset_col_index) 
      {
        const unsigned int node_offset = std::stoi(cell);
        _node_offsets.push_back(node_offset);
        for (unsigned int i = 0; i < node_offset; ++i)
        {
          std::getline(ss, cell, ',');
          _node_ids.push_back(std::stoi(cell));
        }
        break;
      }
    }
  }
}

std::unique_ptr<MeshBase>
BlockFromNodesGeneratorFromFile::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);
  unsigned int idx = 0;
  for (unsigned int i = 0; i < _block_names.size(); i++)
  {
    const unsigned int node_offset = _node_offsets[i];
    const std::vector<unsigned int> node_ids(_node_ids.begin()+idx, _node_ids.begin()+idx+node_offset);
    idx += node_offset;
    const std::vector<BoundaryName> boundary_name(_block_names.begin()+i, _block_names.begin()+i+1);
    const SubdomainName &block_name = _block_names[i]; 
    mesh = BlockFromNodesGeneratorFromFile::GenerateNodeSet(std::move(mesh), node_ids, boundary_name);
  }

  mesh = BlockFromNodesGeneratorFromFile::GenerateSideSet(std::move(mesh));

  for (unsigned int i = 0; i < _block_names.size(); i++)
  {
    const std::vector<BoundaryName> boundary_name(_block_names.begin()+i, _block_names.begin()+i+1);
    const SubdomainName &block_name = _block_names[i]; 
    mesh = BlockFromNodesGeneratorFromFile::GenerateBlock(std::move(mesh), boundary_name, block_name);
  }
  return mesh;
}

std::unique_ptr<MeshBase> 
BlockFromNodesGeneratorFromFile::GenerateNodeSet(
          std::unique_ptr<MeshBase> mesh,
          const std::vector<unsigned int> &node_ids,
          const std::vector<BoundaryName> &boundary_name)
{
  // ***********************************************************************************
  // NodeSetsGenerator***********************************************************

  // Get the BoundaryIDs from the mesh
  std::vector<boundary_id_type> boundary_id =
      MooseMeshUtils::getBoundaryIDs(*mesh, boundary_name, true);

  // Get a reference to our BoundaryInfo object
  BoundaryInfo & boundary_info = mesh->get_boundary_info();

  // add nodes with their ids
  for (const auto & node_id : node_ids)
  {
    // Our mesh may be distributed and this node may not exist on this process
    if (!mesh->query_node_ptr(node_id))
      continue;

    boundary_info.add_node(node_id, boundary_id[0]);
  }

  boundary_info.nodeset_name(boundary_id[0]) = boundary_name[0];

  // This is a terrible hack that we'll want to remove once BMBBG isn't terrible
  if (!_app.getMeshGeneratorSystem().hasBreakMeshByBlockGenerator())
    mesh->set_isnt_prepared();
  return mesh;
}

std::unique_ptr<MeshBase> 
BlockFromNodesGeneratorFromFile::GenerateSideSet(
            std::unique_ptr<MeshBase> mesh)
{
  // ***********************************************************************************
  // SideSetsGenerator***********************************************************
  mesh->get_boundary_info().build_side_list_from_node_list();
  return mesh;
}


std::unique_ptr<MeshBase>
BlockFromNodesGeneratorFromFile::GenerateBlock(
          std::unique_ptr<MeshBase> mesh,
          const std::vector<BoundaryName> &boundary_name,
          const SubdomainName &block_name)
{
  // ***********************************************************************************
  // BlockGenerator***********************************************************
  // Generate a new block id 
  SubdomainID new_block_id = MooseMeshUtils::getNextFreeSubdomainID(*mesh);

  // Make sure our boundary info and parallel counts are setup
  if (!mesh->is_prepared())
  {
    const bool allow_remote_element_removal = mesh->allow_remote_element_removal();
    // We want all of our boundary elements available, so avoid removing them if they haven't
    // already been so
    mesh->allow_remote_element_removal(false);
    mesh->prepare_for_use();
    mesh->allow_remote_element_removal(allow_remote_element_removal);
  }

  auto boundary_id = MooseMeshUtils::getBoundaryIDs(*mesh, boundary_name, true);
  std::set<boundary_id_type> sidesets(boundary_id.begin(), boundary_id.end());

  auto side_list = mesh->get_boundary_info().build_side_list();
  if (!mesh->is_serial() && mesh->comm().size() > 1)
  {
    std::vector<Elem *> elements_to_send;
    unsigned short i_need_boundary_elems = 0;
    for (const auto & [elem_id, side, bc_id] : side_list)
    {
      libmesh_ignore(side);
      if (sidesets.count(bc_id))
      {
        // Whether we have this boundary information through our locally owned element or a ghosted
        // element, we'll need the boundary elements for parallel consistent addition
        i_need_boundary_elems = 1;
        auto * elem = mesh->elem_ptr(elem_id);
        if (elem->processor_id() == mesh->processor_id())
          elements_to_send.push_back(elem);
      }
    }

    std::set<const Elem *, CompareElemIdsByLevel> connected_elements(elements_to_send.begin(),
                                                                     elements_to_send.end());
    std::set<const Node *> connected_nodes;
    reconnect_nodes(connected_elements, connected_nodes);
    std::set<dof_id_type> connected_node_ids;
    for (auto * nd : connected_nodes)
      connected_node_ids.insert(nd->id());

    std::vector<unsigned short> need_boundary_elems(mesh->comm().size());
    mesh->comm().allgather(i_need_boundary_elems, need_boundary_elems);
    std::unordered_map<processor_id_type, decltype(elements_to_send)> push_element_data;
    std::unordered_map<processor_id_type, decltype(connected_nodes)> push_node_data;

    for (const auto pid : index_range(mesh->comm()))
      // Don't need to send to self
      if (pid != mesh->processor_id() && need_boundary_elems[pid])
      {
        if (elements_to_send.size())
          push_element_data[pid] = elements_to_send;
        if (connected_nodes.size())
          push_node_data[pid] = connected_nodes;
      }

    auto node_action_functor = [](processor_id_type, const auto &)
    {
      // Node packing specialization already has unpacked node into mesh, so nothing to do
    };
    Parallel::push_parallel_packed_range(
        mesh->comm(), push_node_data, mesh.get(), node_action_functor);
    auto elem_action_functor = [](processor_id_type, const auto &)
    {
      // Elem packing specialization already has unpacked elem into mesh, so nothing to do
    };
    TIMPI::push_parallel_packed_range(
        mesh->comm(), push_element_data, mesh.get(), elem_action_functor);

    // now that we've gathered everything, we need to rebuild the side list
    side_list = mesh->get_boundary_info().build_side_list();
  }

  std::vector<std::pair<dof_id_type, BlockFromNodesFromFile::ElemSideDouble>> element_sides_on_boundary;
  dof_id_type counter = 0;
  for (const auto & triple : side_list)
    if (sidesets.count(std::get<2>(triple)))
    {
      if (auto elem = mesh->query_elem_ptr(std::get<0>(triple)))
      {
        if (!elem->active())
          mooseError(
              "Only active, level 0 elements can be made interior parents of new level 0 lower-d "
              "elements. Make sure that ",
              type(),
              "s are run before any refinement generators");
        element_sides_on_boundary.push_back(
            std::make_pair(counter, BlockFromNodesFromFile::ElemSideDouble(elem, std::get<1>(triple))));
      }
      ++counter;
    }

  dof_id_type max_elem_id = mesh->max_elem_id();
  unique_id_type max_unique_id = mesh->parallel_max_unique_id();

  // Making an important assumption that at least our boundary elements are the same on all
  // processes even in distributed mesh mode (this is reliant on the correct ghosting functors
  // existing on the mesh)
  for (auto & [i, elem_side] : element_sides_on_boundary)
  {
    Elem * elem = elem_side.elem;

    const auto side = elem_side.side;

    // Build a non-proxy element from this side.
    std::unique_ptr<Elem> side_elem(elem->build_side_ptr(side, /*proxy=*/false));

    // The side will be added with the same processor id as the parent.
    side_elem->processor_id() = elem->processor_id();

    // Add subdomain ID
    side_elem->subdomain_id() = new_block_id;

    // Also assign the side's interior parent, so it is always
    // easy to figure out the Elem we came from.
    side_elem->set_interior_parent(elem);

    // Add id
    side_elem->set_id(max_elem_id + i);
    side_elem->set_unique_id(max_unique_id + i);

    // Finally, add the lower-dimensional element to the Mesh.
    mesh->add_elem(side_elem.release());
  };

  // Assign block name
  mesh->subdomain_name(new_block_id) = block_name;

  const bool skip_partitioning_old = mesh->skip_partitioning();
  mesh->skip_partitioning(true);
  mesh->prepare_for_use();
  mesh->skip_partitioning(skip_partitioning_old);

  return dynamic_pointer_cast<MeshBase>(mesh);
}
