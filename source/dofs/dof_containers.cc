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

template <int dim, int spacedim>
FaceDoFContainer<dim, spacedim>::FaceDoFContainer(const DoFHandler<dim>* handler, const FiniteElement<dim>* fe_)
: dof_handler(handler)
, fe(fe_)
{
}

bool compare_dofs(std::array<unsigned int, 4 > arr1, std::array<unsigned int, 4 > arr2){
  return (arr1[0] == arr2[1] && arr1[1] == arr2[0]);
}	


template <int dim, int spacedim>
void FaceDoFContainer<dim, spacedim>::get_dofs_and_indices(std::vector<std::array<unsigned int, 4>> &out)
{
  const unsigned int dofs_per_cell = (*fe).dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto& cell : (*dof_handler).active_cell_iterators()){
    cell->get_dof_indices(local_dof_indices);
    std::array<unsigned int, 4> local_array;
    for(int i = 0; i < 2*dim; ++i){
      local_array[0] = local_dof_indices[0];
      auto neighbor = cell->neighbor(i);
      neighbor->get_dof_indices(local_dof_indices);
      local_array[1] = local_dof_indices[0];
      local_array[2] = i;
      local_array[3] = cell->neighbor_face_no(i);
    }
    out.push_back(local_array); 
  }
  out.erase(std::unique(out.begin(), out.end(), compare_dofs));
}

// Explicit template instantiation
template class CellDoFContainer<1,1>;
template class CellDoFContainer<2,2>;
template class CellDoFContainer<3,3>;

template class FaceDoFContainer<1,1>;
template class FaceDoFContainer<2,2>;
template class FaceDoFContainer<3,3>;
