#ifndef dealii_dof_containers_h
#define dealii_dof_containers_h

#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
// And this is the file in which the functions are declared that create grids:
#include <deal.II/grid/grid_generator.h>

// The next three files contain classes which are needed for loops over all
// cells and to get the information from the cell objects. The first two have
// been used before to get geometric information from cells; the last one is
// new and provides information about the degrees of freedom local to a cell:
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

// And this file is needed for the creation of sparsity patterns of sparse
// matrices, as shown in previous examples:
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

template <int dofs_per_cell, int dim, int spacedim=dim>
class CellDoFContainer
{
public:

  CellDoFContainer() = default;	
  void initialize(const DoFHandler<dim,spacedim> &dof_handler);
  std::vector<unsigned int> get_cell_dofs();
  void fill_coordinates(std::vector<double> &coordinates, const DoFHandler<dim,spacedim> &dof_handler, const FiniteElement<dim,spacedim> &fe);
private:

  std::vector<unsigned int> cell_dof_vector;
};

template <int dofs_per_cell, int dim, int spacedim=dim>
class FaceDoFContainer
{
public:

  FaceDoFContainer() = default;
  void initialize(const DoFHandler<dim,spacedim> &dof_handler);
  std::vector<std::array<unsigned int,4>> get_face_dofs();
  std::vector<std::array<unsigned int,2>> get_boundary_face_dofs();
private:

  std::vector<std::array<unsigned int, 4>> interior_face_dof_array;
  // additional array for boundary faces
  std::vector<std::array<unsigned int, 2>> boundary_face_dof_array;
};

template <int dofs_per_cell, int dim, int spacedim>
void CellDoFContainer<dofs_per_cell, dim, spacedim>::initialize(const DoFHandler<dim,spacedim> &dof_handler)
{
  if(dof_handler.get_fe().dofs_per_cell != dofs_per_cell)
    std::cout << "Error: Inconsistent #dof for CellDoFContainer" << std::endl;
  else{	
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    for (const auto& cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_indices(local_dof_indices);
      cell_dof_vector.push_back(local_dof_indices[0]);
    }
  }  
}

template<int dofs_per_cell, int dim, int spacedim>
void CellDoFContainer<dofs_per_cell, dim, spacedim>::fill_coordinates(std::vector<double> &coordinates, const DoFHandler<dim,spacedim> &dof_handler, const FiniteElement<dim,spacedim> &fe)
{
  FESystem<dim,spacedim> fe_sys(fe,2);
  int sys_dofs_per_cell = fe_sys.dofs_per_cell;
  for (const auto& cell : dof_handler.active_cell_iterators()){
    for(int i = 0; i < sys_dofs_per_cell / dim; ++i){
      Point<dim> vertex_coords = cell->vertex(i);
      for(int j = 0; j < dim; ++j){
        coordinates.push_back(vertex_coords(j));	
      }
    }	
  }
}

template <int dofs_per_cell, int dim, int spacedim>
std::vector<unsigned int> CellDoFContainer<dofs_per_cell, dim, spacedim>::get_cell_dofs()
{
  return cell_dof_vector;
}

// operator()
// size
// num_dofs gegen dofs_per_cell checken


template <int dofs_per_cell, int dim, int spacedim>
void FaceDoFContainer<dofs_per_cell, dim, spacedim>::initialize(const DoFHandler<dim,spacedim> &dof_handler)
{
  if(dof_handler.get_fe().dofs_per_cell != dofs_per_cell)
    std::cout << "Error: Inconsistent #dof for FaceDoFContainer" << std::endl;
  else{
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    for (const auto& cell : dof_handler.active_cell_iterators()){
      cell->get_dof_indices(local_dof_indices);
      std::array<unsigned int, 4> local_array;
      local_array[0] = local_dof_indices[0];
      for(int i = 0; i < 2*dim; ++i){
        if(cell->at_boundary(i)){
          std::array<unsigned int, 2> local_boundary_array;
  	  local_boundary_array[0] = local_array[0];
	  local_boundary_array[1] = i;
	  boundary_face_dof_array.push_back(local_boundary_array);
        }
        auto neighbor = cell->neighbor(i); 
        neighbor->get_dof_indices(local_dof_indices);
        if(local_dof_indices[0] > local_array[0]){
          local_array[1] = local_dof_indices[0];
          local_array[2] = i;
          local_array[3] = cell->neighbor_face_no(i);
          face_dof_array.push_back(local_array);
        } 
      }
    }
  }  
}

template <int dofs_per_cell, int dim, int spacedim>
std::vector<std::array<unsigned int, 4>> FaceDoFContainer<dofs_per_cell, dim, spacedim>::get_face_dofs()
{
  return face_dof_array;
}

template <int dofs_per_cell, int dim, int spacedim>
std::vector<std::array<unsigned int, 2>> FaceDoFContainer<dofs_per_cell, dim, spacedim>::get_boundary_face_dofs()
{
  return boundary_face_dof_array;
}


DEAL_II_NAMESPACE_CLOSE
#endif
