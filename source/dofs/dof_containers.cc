#include <deal.II/dofs/dof_containers.h>

template <int dim, int spacedim>
CellDoFContainer<dim, spacedim>::CellDoFContainer(const DoFHandler<dim>* handler, const FiniteElement<dim>* fe_)
: dof_handler(handler)
, fe(fe_)
{
}

template <int dim, int spacedim>
void CellDoFContainer<dim, spacedim>::get_dofs(std::vector<int> &out)
{
  const unsigned int dofs_per_cell = (*fe).dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto& cell : (*dof_handler).active_cell_iterators())
  {
    cell->get_dof_indices(local_dof_indices);
    out.push_back(local_dof_indices[0]); 
  } 
}

template class CellDoFContainer<1,1>;
template class CellDoFContainer<2,2>;
template class CellDoFContainer<3,3>;

