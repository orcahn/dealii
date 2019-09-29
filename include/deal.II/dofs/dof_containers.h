#ifndef dealii_dof_containers_h
#define dealii_dof_containers_h

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

// This file contains the description of the Lagrange interpolation finite
// element:
#include <deal.II/fe/fe_dgq.h>

#include <fstream>
#include <iostream>
using namespace dealii;

template <int dim, int spacedim>
class CellDoFContainer
{
public:
  
  CellDoFContainer(DoFHandler<dim> handler, int degree);

  void get_dofs(std::vector<int> &out);
private:

  DoFHandler<dim> dof_handler;
  FE_DGQ<dim> fe; 
};
#endif
