// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/table.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_MPI
#  include <deal.II/lac/la_parallel_vector.h>
#endif

#include <deal.II/base/std_cxx14/memory.h>

#include <algorithm>
#include <array>
#include <memory>
#include <numeric>

DEAL_II_NAMESPACE_OPEN



namespace DoFTools
{
  namespace internal
  {
    namespace
    {
      inline bool
      check_master_dof_list(
        const FullMatrix<double> &                  face_interpolation_matrix,
        const std::vector<types::global_dof_index> &master_dof_list)
      {
        const unsigned int N = master_dof_list.size();

        FullMatrix<double> tmp(N, N);
        for (unsigned int i = 0; i < N; ++i)
          for (unsigned int j = 0; j < N; ++j)
            tmp(i, j) = face_interpolation_matrix(master_dof_list[i], j);

        // then use the algorithm from FullMatrix::gauss_jordan on this matrix
        // to find out whether it is singular. the algorithm there does pivoting
        // and at the end swaps rows back into their proper order -- we omit
        // this step here, since we don't care about the inverse matrix, all we
        // care about is whether the matrix is regular or singular

        // first get an estimate of the size of the elements of this matrix, for
        // later checks whether the pivot element is large enough, or whether we
        // have to fear that the matrix is not regular
        double diagonal_sum = 0;
        for (unsigned int i = 0; i < N; ++i)
          diagonal_sum += std::fabs(tmp(i, i));
        const double typical_diagonal_element = diagonal_sum / N;

        // initialize the array that holds the permutations that we find during
        // pivot search
        std::vector<unsigned int> p(N);
        for (unsigned int i = 0; i < N; ++i)
          p[i] = i;

        for (unsigned int j = 0; j < N; ++j)
          {
            // pivot search: search that part of the line on and right of the
            // diagonal for the largest element
            double       max = std::fabs(tmp(j, j));
            unsigned int r   = j;
            for (unsigned int i = j + 1; i < N; ++i)
              {
                if (std::fabs(tmp(i, j)) > max)
                  {
                    max = std::fabs(tmp(i, j));
                    r   = i;
                  }
              }
            // check whether the pivot is too small. if that is the case, then
            // the matrix is singular and we shouldn't use this set of master
            // dofs
            if (max < 1.e-12 * typical_diagonal_element)
              return false;

            // row interchange
            if (r > j)
              {
                for (unsigned int k = 0; k < N; ++k)
                  std::swap(tmp(j, k), tmp(r, k));

                std::swap(p[j], p[r]);
              }

            // transformation
            const double hr = 1. / tmp(j, j);
            tmp(j, j)       = hr;
            for (unsigned int k = 0; k < N; ++k)
              {
                if (k == j)
                  continue;
                for (unsigned int i = 0; i < N; ++i)
                  {
                    if (i == j)
                      continue;
                    tmp(i, k) -= tmp(i, j) * tmp(j, k) * hr;
                  }
              }
            for (unsigned int i = 0; i < N; ++i)
              {
                tmp(i, j) *= hr;
                tmp(j, i) *= -hr;
              }
            tmp(j, j) = hr;
          }

        // everything went fine, so we can accept this set of master dofs (at
        // least as far as they have already been collected)
        return true;
      }



      /**
       * When restricting, on a face, the degrees of freedom of fe1 to the space
       * described by fe2 (for example for the complex case described
       * in the @ref hp_paper "hp paper"), we have to select fe2.dofs_per_face
       * out of the fe1.dofs_per_face face DoFs as the
       * master DoFs, and the rest become slave dofs. This function selects
       * which ones will be masters, and which ones will be slaves.
       *
       * The function assumes that master_dofs already has size
       * fe1.dofs_per_face. After the function, exactly fe2.dofs_per_face
       * entries will be true.
       *
       * The function is a bit complicated since it has to figure out a set a
       * DoFs so that the corresponding rows in the face interpolation matrix
       * are all linearly independent. we have a good heuristic (see the
       * function body) for selecting these rows, but there are cases where this
       * fails and we have to pick them differently. what we do is to run the
       * heuristic and then go back to determine whether we have a set of rows
       * with full row rank. if this isn't the case, go back and select dofs
       * differently
       */
      template <int dim, int spacedim>
      void
      select_master_dofs_for_face_restriction(
        const FiniteElement<dim, spacedim> &fe1,
        const FiniteElement<dim, spacedim> &fe2,
        const FullMatrix<double> &          face_interpolation_matrix,
        std::vector<bool> &                 master_dof_mask)
      {
        Assert(fe1.dofs_per_face >= fe2.dofs_per_face, ExcInternalError());
        AssertDimension(master_dof_mask.size(), fe1.dofs_per_face);

        Assert(fe2.dofs_per_vertex <= fe1.dofs_per_vertex, ExcInternalError());
        Assert(fe2.dofs_per_line <= fe1.dofs_per_line, ExcInternalError());
        Assert((dim < 3) || (fe2.dofs_per_quad <= fe1.dofs_per_quad),
               ExcInternalError());

        // the idea here is to designate as many DoFs in fe1 per object (vertex,
        // line, quad) as master as there are such dofs in fe2 (indices are int,
        // because we want to avoid the 'unsigned int < 0 is always false
        // warning for the cases at the bottom in 1d and 2d)
        //
        // as mentioned in the paper, it is not always easy to find a set of
        // master dofs that produces an invertible matrix. to this end, we check
        // in each step whether the matrix is still invertible and simply
        // discard this dof if the matrix is not invertible anymore.
        //
        // the cases where we did have trouble in the past were with adding more
        // quad dofs when Q3 and Q4 elements meet at a refined face in 3d (see
        // the hp/crash_12 test that tests that we can do exactly this, and
        // failed before we had code to compensate for this case). the other
        // case are system elements: if we have say a Q1Q2 vs a Q2Q3 element,
        // then we can't just take all master dofs on a line from a single base
        // element, since the shape functions of that base element are
        // independent of that of the other one. this latter case shows up when
        // running hp/hp_constraints_q_system_06

        std::vector<types::global_dof_index> master_dof_list;
        unsigned int                         index = 0;
        for (int v = 0;
             v < static_cast<signed int>(GeometryInfo<dim>::vertices_per_face);
             ++v)
          {
            unsigned int dofs_added = 0;
            unsigned int i          = 0;
            while (dofs_added < fe2.dofs_per_vertex)
              {
                // make sure that we were able to find a set of master dofs and
                // that the code down below didn't just reject all our efforts
                Assert(i < fe1.dofs_per_vertex, ExcInternalError());

                // tentatively push this vertex dof
                master_dof_list.push_back(index + i);

                // then see what happens. if it succeeds, fine
                if (check_master_dof_list(face_interpolation_matrix,
                                          master_dof_list) == true)
                  ++dofs_added;
                else
                  // well, it didn't. simply pop that dof from the list again
                  // and try with the next dof
                  master_dof_list.pop_back();

                // forward counter by one
                ++i;
              }
            index += fe1.dofs_per_vertex;
          }

        for (int l = 0;
             l < static_cast<signed int>(GeometryInfo<dim>::lines_per_face);
             ++l)
          {
            // same algorithm as above
            unsigned int dofs_added = 0;
            unsigned int i          = 0;
            while (dofs_added < fe2.dofs_per_line)
              {
                Assert(i < fe1.dofs_per_line, ExcInternalError());

                master_dof_list.push_back(index + i);
                if (check_master_dof_list(face_interpolation_matrix,
                                          master_dof_list) == true)
                  ++dofs_added;
                else
                  master_dof_list.pop_back();

                ++i;
              }
            index += fe1.dofs_per_line;
          }

        for (int q = 0;
             q < static_cast<signed int>(GeometryInfo<dim>::quads_per_face);
             ++q)
          {
            // same algorithm as above
            unsigned int dofs_added = 0;
            unsigned int i          = 0;
            while (dofs_added < fe2.dofs_per_quad)
              {
                Assert(i < fe1.dofs_per_quad, ExcInternalError());

                master_dof_list.push_back(index + i);
                if (check_master_dof_list(face_interpolation_matrix,
                                          master_dof_list) == true)
                  ++dofs_added;
                else
                  master_dof_list.pop_back();

                ++i;
              }
            index += fe1.dofs_per_quad;
          }

        AssertDimension(index, fe1.dofs_per_face);
        AssertDimension(master_dof_list.size(), fe2.dofs_per_face);

        // finally copy the list into the mask
        std::fill(master_dof_mask.begin(), master_dof_mask.end(), false);
        for (const auto dof : master_dof_list)
          master_dof_mask[dof] = true;
      }



      /**
       * Make sure that the mask exists that determines which dofs will be the
       * masters on refined faces where an fe1 and a fe2 meet.
       */
      template <int dim, int spacedim>
      void
      ensure_existence_of_master_dof_mask(
        const FiniteElement<dim, spacedim> &fe1,
        const FiniteElement<dim, spacedim> &fe2,
        const FullMatrix<double> &          face_interpolation_matrix,
        std::unique_ptr<std::vector<bool>> &master_dof_mask)
      {
        if (master_dof_mask == nullptr)
          {
            master_dof_mask =
              std_cxx14::make_unique<std::vector<bool>>(fe1.dofs_per_face);
            select_master_dofs_for_face_restriction(fe1,
                                                    fe2,
                                                    face_interpolation_matrix,
                                                    *master_dof_mask);
          }
      }



      /**
       * Make sure that the given @p face_interpolation_matrix pointer points
       * to a valid matrix. If the pointer is zero beforehand, create an entry
       * with the correct data. If it is nonzero, don't touch it.
       */
      template <int dim, int spacedim>
      void
      ensure_existence_of_face_matrix(
        const FiniteElement<dim, spacedim> & fe1,
        const FiniteElement<dim, spacedim> & fe2,
        std::unique_ptr<FullMatrix<double>> &matrix)
      {
        if (matrix == nullptr)
          {
            matrix =
              std_cxx14::make_unique<FullMatrix<double>>(fe2.dofs_per_face,
                                                         fe1.dofs_per_face);
            fe1.get_face_interpolation_matrix(fe2, *matrix);
          }
      }



      /**
       * Same, but for subface interpolation matrices.
       */
      template <int dim, int spacedim>
      void
      ensure_existence_of_subface_matrix(
        const FiniteElement<dim, spacedim> & fe1,
        const FiniteElement<dim, spacedim> & fe2,
        const unsigned int                   subface,
        std::unique_ptr<FullMatrix<double>> &matrix)
      {
        if (matrix == nullptr)
          {
            matrix =
              std_cxx14::make_unique<FullMatrix<double>>(fe2.dofs_per_face,
                                                         fe1.dofs_per_face);
            fe1.get_subface_interpolation_matrix(fe2, subface, *matrix);
          }
      }



      /**
       * Given the face interpolation matrix between two elements, split it into
       * its master and slave parts and invert the master part as
       * explained in the @ref hp_paper "hp paper".
       */
      void
      ensure_existence_of_split_face_matrix(
        const FullMatrix<double> &face_interpolation_matrix,
        const std::vector<bool> & master_dof_mask,
        std::unique_ptr<std::pair<FullMatrix<double>, FullMatrix<double>>>
          &split_matrix)
      {
        AssertDimension(master_dof_mask.size(), face_interpolation_matrix.m());
        Assert(std::count(master_dof_mask.begin(),
                          master_dof_mask.end(),
                          true) ==
                 static_cast<signed int>(face_interpolation_matrix.n()),
               ExcInternalError());

        if (split_matrix == nullptr)
          {
            split_matrix = std_cxx14::make_unique<
              std::pair<FullMatrix<double>, FullMatrix<double>>>();

            const unsigned int n_master_dofs = face_interpolation_matrix.n();
            const unsigned int n_dofs        = face_interpolation_matrix.m();

            Assert(n_master_dofs <= n_dofs, ExcInternalError());

            // copy and invert the master component, copy the slave component
            split_matrix->first.reinit(n_master_dofs, n_master_dofs);
            split_matrix->second.reinit(n_dofs - n_master_dofs, n_master_dofs);

            unsigned int nth_master_dof = 0, nth_slave_dof = 0;

            for (unsigned int i = 0; i < n_dofs; ++i)
              if (master_dof_mask[i] == true)
                {
                  for (unsigned int j = 0; j < n_master_dofs; ++j)
                    split_matrix->first(nth_master_dof, j) =
                      face_interpolation_matrix(i, j);
                  ++nth_master_dof;
                }
              else
                {
                  for (unsigned int j = 0; j < n_master_dofs; ++j)
                    split_matrix->second(nth_slave_dof, j) =
                      face_interpolation_matrix(i, j);
                  ++nth_slave_dof;
                }

            AssertDimension(nth_master_dof, n_master_dofs);
            AssertDimension(nth_slave_dof, n_dofs - n_master_dofs);

            // TODO[WB]: We should make sure very small entries are removed
            // after inversion
            split_matrix->first.gauss_jordan();
          }
      }


      // a template that can determine statically whether a given DoFHandler
      // class supports different finite element elements
      template <typename>
      struct DoFHandlerSupportsDifferentFEs
      {
        static const bool value = true;
      };


      template <int dim, int spacedim>
      struct DoFHandlerSupportsDifferentFEs<dealii::DoFHandler<dim, spacedim>>
      {
        static const bool value = false;
      };


      /**
       * A function that returns how many different finite elements a dof
       * handler uses. This is one for non-hp DoFHandlers and
       * dof_handler.get_fe().size() for the hp-versions.
       */
      template <int dim, int spacedim>
      unsigned int
      n_finite_elements(
        const dealii::hp::DoFHandler<dim, spacedim> &dof_handler)
      {
        return dof_handler.get_fe_collection().size();
      }


      template <typename DoFHandlerType>
      unsigned int
      n_finite_elements(const DoFHandlerType &)
      {
        return 1;
      }



      /**
       * Copy constraints into an AffineConstraints object.
       *
       * This function removes zero constraints and those, which constrain a DoF
       * which was already eliminated in one of the previous steps of the hp
       * hanging node procedure.
       *
       * It also suppresses very small entries in the AffineConstraints object
       * to avoid making the sparsity pattern fuller than necessary.
       */
      template <typename number1, typename number2>
      void
      filter_constraints(
        const std::vector<types::global_dof_index> &master_dofs,
        const std::vector<types::global_dof_index> &slave_dofs,
        const FullMatrix<number1> &                 face_constraints,
        AffineConstraints<number2> &                constraints)
      {
        Assert(face_constraints.n() == master_dofs.size(),
               ExcDimensionMismatch(master_dofs.size(), face_constraints.n()));
        Assert(face_constraints.m() == slave_dofs.size(),
               ExcDimensionMismatch(slave_dofs.size(), face_constraints.m()));

        const unsigned int n_master_dofs = master_dofs.size();
        const unsigned int n_slave_dofs  = slave_dofs.size();

        // check for a couple conditions that happened in parallel distributed
        // mode
        for (unsigned int row = 0; row != n_slave_dofs; ++row)
          Assert(slave_dofs[row] != numbers::invalid_dof_index,
                 ExcInternalError());
        for (unsigned int col = 0; col != n_master_dofs; ++col)
          Assert(master_dofs[col] != numbers::invalid_dof_index,
                 ExcInternalError());


        for (unsigned int row = 0; row != n_slave_dofs; ++row)
          if (constraints.is_constrained(slave_dofs[row]) == false)
            {
              bool constraint_already_satisfied = false;

              // Check if we have an identity constraint, which is already
              // satisfied by unification of the corresponding global dof
              // indices
              for (unsigned int i = 0; i < n_master_dofs; ++i)
                if (face_constraints(row, i) == 1.0)
                  if (master_dofs[i] == slave_dofs[row])
                    {
                      constraint_already_satisfied = true;
                      break;
                    }

              if (constraint_already_satisfied == false)
                {
                  // add up the absolute values of all constraints in this line
                  // to get a measure of their absolute size
                  number1 abs_sum = 0;
                  for (unsigned int i = 0; i < n_master_dofs; ++i)
                    abs_sum += std::abs(face_constraints(row, i));

                  // then enter those constraints that are larger than
                  // 1e-14*abs_sum. everything else probably originated from
                  // inexact inversion of matrices and similar effects. having
                  // those constraints in here will only lead to problems
                  // because it makes sparsity patterns fuller than necessary
                  // without producing any significant effect
                  constraints.add_line(slave_dofs[row]);
                  for (unsigned int i = 0; i < n_master_dofs; ++i)
                    if ((face_constraints(row, i) != 0) &&
                        (std::fabs(face_constraints(row, i)) >=
                         1e-14 * abs_sum))
                      constraints.add_entry(slave_dofs[row],
                                            master_dofs[i],
                                            face_constraints(row, i));
                  constraints.set_inhomogeneity(slave_dofs[row], 0.);
                }
            }
      }

    } // namespace


    template <typename number>
    void
    make_hp_hanging_node_constraints(const dealii::DoFHandler<1> &,
                                     AffineConstraints<number> &)
    {
      // nothing to do for regular dof handlers in 1d
    }


    template <typename number>
    void
    make_oldstyle_hanging_node_constraints(const dealii::DoFHandler<1> &,
                                           AffineConstraints<number> &,
                                           std::integral_constant<int, 1>)
    {
      // nothing to do for regular dof handlers in 1d
    }


    template <typename number>
    void
    make_hp_hanging_node_constraints(
      const dealii::hp::DoFHandler<1> & /*dof_handler*/,
      AffineConstraints<number> & /*constraints*/)
    {
      // we may have to compute constraints for vertices. gotta think about that
      // a bit more

      // TODO[WB]: think about what to do here...
    }


    template <typename number>
    void
    make_oldstyle_hanging_node_constraints(
      const dealii::hp::DoFHandler<1> & /*dof_handler*/,
      AffineConstraints<number> & /*constraints*/,
      std::integral_constant<int, 1>)
    {
      // we may have to compute constraints for vertices. gotta think about that
      // a bit more

      // TODO[WB]: think about what to do here...
    }


    template <typename number>
    void
    make_hp_hanging_node_constraints(const dealii::DoFHandler<1, 2> &,
                                     AffineConstraints<number> &)
    {
      // nothing to do for regular dof handlers in 1d
    }


    template <typename number>
    void
    make_oldstyle_hanging_node_constraints(const dealii::DoFHandler<1, 2> &,
                                           AffineConstraints<number> &,
                                           std::integral_constant<int, 1>)
    {
      // nothing to do for regular dof handlers in 1d
    }


    template <typename number>
    void
    make_hp_hanging_node_constraints(const dealii::DoFHandler<1, 3> &,
                                     AffineConstraints<number> &)
    {
      // nothing to do for regular dof handlers in 1d
    }


    template <typename number>
    void
    make_oldstyle_hanging_node_constraints(const dealii::DoFHandler<1, 3> &,
                                           AffineConstraints<number> &,
                                           std::integral_constant<int, 1>)
    {
      // nothing to do for regular dof handlers in 1d
    }



    template <typename DoFHandlerType, typename number>
    void
    make_oldstyle_hanging_node_constraints(
      const DoFHandlerType &     dof_handler,
      AffineConstraints<number> &constraints,
      std::integral_constant<int, 2>)
    {
      const unsigned int dim = 2;

      const unsigned int spacedim = DoFHandlerType::space_dimension;

      std::vector<types::global_dof_index> dofs_on_mother;
      std::vector<types::global_dof_index> dofs_on_children;

      // loop over all lines; only on lines there can be constraints. We do so
      // by looping over all active cells and checking whether any of the faces
      // are refined which can only be from the neighboring cell because this
      // one is active. In that case, the face is subject to constraints
      //
      // note that even though we may visit a face twice if the neighboring
      // cells are equally refined, we can only visit each face with hanging
      // nodes once
      typename DoFHandlerType::active_cell_iterator cell = dof_handler
                                                             .begin_active(),
                                                    endc = dof_handler.end();
      for (; cell != endc; ++cell)
        {
          // artificial cells can at best neighbor ghost cells, but we're not
          // interested in these interfaces
          if (cell->is_artificial())
            continue;

          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               ++face)
            if (cell->face(face)->has_children())
              {
                // in any case, faces can have at most two active fe indices,
                // but here the face can have only one (namely the same as that
                // from the cell we're sitting on), and each of the children can
                // have only one as well. check this
                Assert(cell->face(face)->n_active_fe_indices() == 1,
                       ExcInternalError());
                Assert(cell->face(face)->fe_index_is_active(
                         cell->active_fe_index()) == true,
                       ExcInternalError());
                for (unsigned int c = 0; c < cell->face(face)->n_children();
                     ++c)
                  if (!cell->neighbor_child_on_subface(face, c)
                         ->is_artificial())
                    Assert(cell->face(face)->child(c)->n_active_fe_indices() ==
                             1,
                           ExcInternalError());

                // right now, all that is implemented is the case that both
                // sides use the same fe
                for (unsigned int c = 0; c < cell->face(face)->n_children();
                     ++c)
                  if (!cell->neighbor_child_on_subface(face, c)
                         ->is_artificial())
                    Assert(cell->face(face)->child(c)->fe_index_is_active(
                             cell->active_fe_index()) == true,
                           ExcNotImplemented());

                // ok, start up the work
                const FiniteElement<dim, spacedim> &fe = cell->get_fe();
                const unsigned int fe_index = cell->active_fe_index();

                const unsigned int n_dofs_on_mother =
                                     2 * fe.dofs_per_vertex + fe.dofs_per_line,
                                   n_dofs_on_children =
                                     fe.dofs_per_vertex + 2 * fe.dofs_per_line;

                dofs_on_mother.resize(n_dofs_on_mother);
                // we might not use all of those in case of artificial cells, so
                // do not resize(), but reserve() and use push_back later.
                dofs_on_children.clear();
                dofs_on_children.reserve(n_dofs_on_children);

                Assert(n_dofs_on_mother == fe.constraints().n(),
                       ExcDimensionMismatch(n_dofs_on_mother,
                                            fe.constraints().n()));
                Assert(n_dofs_on_children == fe.constraints().m(),
                       ExcDimensionMismatch(n_dofs_on_children,
                                            fe.constraints().m()));

                const typename DoFHandlerType::line_iterator this_face =
                  cell->face(face);

                // fill the dofs indices. Use same enumeration scheme as in
                // @p{FiniteElement::constraints()}
                unsigned int next_index = 0;
                for (unsigned int vertex = 0; vertex < 2; ++vertex)
                  for (unsigned int dof = 0; dof != fe.dofs_per_vertex; ++dof)
                    dofs_on_mother[next_index++] =
                      this_face->vertex_dof_index(vertex, dof, fe_index);
                for (unsigned int dof = 0; dof != fe.dofs_per_line; ++dof)
                  dofs_on_mother[next_index++] =
                    this_face->dof_index(dof, fe_index);
                AssertDimension(next_index, dofs_on_mother.size());

                for (unsigned int dof = 0; dof != fe.dofs_per_vertex; ++dof)
                  dofs_on_children.push_back(
                    this_face->child(0)->vertex_dof_index(1, dof, fe_index));
                for (unsigned int child = 0; child < 2; ++child)
                  {
                    // skip artificial cells
                    if (cell->neighbor_child_on_subface(face, child)
                          ->is_artificial())
                      continue;
                    for (unsigned int dof = 0; dof != fe.dofs_per_line; ++dof)
                      dofs_on_children.push_back(
                        this_face->child(child)->dof_index(dof, fe_index));
                  }
                // note: can get fewer DoFs when we have artificial cells
                Assert(dofs_on_children.size() <= n_dofs_on_children,
                       ExcInternalError());

                // for each row in the AffineConstraints object for this line:
                for (unsigned int row = 0; row != dofs_on_children.size();
                     ++row)
                  {
                    constraints.add_line(dofs_on_children[row]);
                    for (unsigned int i = 0; i != dofs_on_mother.size(); ++i)
                      constraints.add_entry(dofs_on_children[row],
                                            dofs_on_mother[i],
                                            fe.constraints()(row, i));

                    constraints.set_inhomogeneity(dofs_on_children[row], 0.);
                  }
              }
            else
              {
                // this face has no children, but it could still be that it is
                // shared by two cells that use a different fe index. check a
                // couple of things, but ignore the case that the neighbor is an
                // artificial cell
                if (!cell->at_boundary(face) &&
                    !cell->neighbor(face)->is_artificial())
                  {
                    Assert(cell->face(face)->n_active_fe_indices() == 1,
                           ExcNotImplemented());
                    Assert(cell->face(face)->fe_index_is_active(
                             cell->active_fe_index()) == true,
                           ExcInternalError());
                  }
              }
        }
    }



    template <typename DoFHandlerType, typename number>
    void
    make_oldstyle_hanging_node_constraints(
      const DoFHandlerType &     dof_handler,
      AffineConstraints<number> &constraints,
      std::integral_constant<int, 3>)
    {
      const unsigned int dim = 3;

      std::vector<types::global_dof_index> dofs_on_mother;
      std::vector<types::global_dof_index> dofs_on_children;

      // loop over all quads; only on quads there can be constraints. We do so
      // by looping over all active cells and checking whether any of the faces
      // are refined which can only be from the neighboring cell because this
      // one is active. In that case, the face is subject to constraints
      //
      // note that even though we may visit a face twice if the neighboring
      // cells are equally refined, we can only visit each face with hanging
      // nodes once
      typename DoFHandlerType::active_cell_iterator cell = dof_handler
                                                             .begin_active(),
                                                    endc = dof_handler.end();
      for (; cell != endc; ++cell)
        {
          // artificial cells can at best neighbor ghost cells, but we're not
          // interested in these interfaces
          if (cell->is_artificial())
            continue;

          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               ++face)
            if (cell->face(face)->has_children())
              {
                // first of all, make sure that we treat a case which is
                // possible, i.e. either no dofs on the face at all or no
                // anisotropic refinement
                if (cell->get_fe().dofs_per_face == 0)
                  continue;

                Assert(cell->face(face)->refinement_case() ==
                         RefinementCase<dim - 1>::isotropic_refinement,
                       ExcNotImplemented());

                // in any case, faces can have at most two active fe indices,
                // but here the face can have only one (namely the same as that
                // from the cell we're sitting on), and each of the children can
                // have only one as well. check this
                AssertDimension(cell->face(face)->n_active_fe_indices(), 1);
                Assert(cell->face(face)->fe_index_is_active(
                         cell->active_fe_index()) == true,
                       ExcInternalError());
                for (unsigned int c = 0; c < cell->face(face)->n_children();
                     ++c)
                  if (!cell->neighbor_child_on_subface(face, c)
                         ->is_artificial())
                    AssertDimension(
                      cell->face(face)->child(c)->n_active_fe_indices(), 1);

                // right now, all that is implemented is the case that both
                // sides use the same fe, and not only that but also that all
                // lines bounding this face and the children have the same fe
                for (unsigned int c = 0; c < cell->face(face)->n_children();
                     ++c)
                  if (!cell->neighbor_child_on_subface(face, c)
                         ->is_artificial())
                    {
                      Assert(cell->face(face)->child(c)->fe_index_is_active(
                               cell->active_fe_index()) == true,
                             ExcNotImplemented());
                      for (unsigned int e = 0; e < 4; ++e)
                        {
                          Assert(cell->face(face)
                                     ->child(c)
                                     ->line(e)
                                     ->n_active_fe_indices() == 1,
                                 ExcNotImplemented());
                          Assert(cell->face(face)
                                     ->child(c)
                                     ->line(e)
                                     ->fe_index_is_active(
                                       cell->active_fe_index()) == true,
                                 ExcNotImplemented());
                        }
                    }
                for (unsigned int e = 0; e < 4; ++e)
                  {
                    Assert(cell->face(face)->line(e)->n_active_fe_indices() ==
                             1,
                           ExcNotImplemented());
                    Assert(cell->face(face)->line(e)->fe_index_is_active(
                             cell->active_fe_index()) == true,
                           ExcNotImplemented());
                  }

                // ok, start up the work
                const FiniteElement<dim> &fe       = cell->get_fe();
                const unsigned int        fe_index = cell->active_fe_index();

                const unsigned int n_dofs_on_mother = fe.dofs_per_face;
                const unsigned int n_dofs_on_children =
                  (5 * fe.dofs_per_vertex + 12 * fe.dofs_per_line +
                   4 * fe.dofs_per_quad);

                // TODO[TL]: think about this and the following in case of
                // anisotropic refinement

                dofs_on_mother.resize(n_dofs_on_mother);
                // we might not use all of those in case of artificial cells, so
                // do not resize(), but reserve() and use push_back later.
                dofs_on_children.clear();
                dofs_on_children.reserve(n_dofs_on_children);

                Assert(n_dofs_on_mother == fe.constraints().n(),
                       ExcDimensionMismatch(n_dofs_on_mother,
                                            fe.constraints().n()));
                Assert(n_dofs_on_children == fe.constraints().m(),
                       ExcDimensionMismatch(n_dofs_on_children,
                                            fe.constraints().m()));

                const typename DoFHandlerType::face_iterator this_face =
                  cell->face(face);

                // fill the dofs indices. Use same enumeration scheme as in
                // @p{FiniteElement::constraints()}
                unsigned int next_index = 0;
                for (unsigned int vertex = 0; vertex < 4; ++vertex)
                  for (unsigned int dof = 0; dof != fe.dofs_per_vertex; ++dof)
                    dofs_on_mother[next_index++] =
                      this_face->vertex_dof_index(vertex, dof, fe_index);
                for (unsigned int line = 0; line < 4; ++line)
                  for (unsigned int dof = 0; dof != fe.dofs_per_line; ++dof)
                    dofs_on_mother[next_index++] =
                      this_face->line(line)->dof_index(dof, fe_index);
                for (unsigned int dof = 0; dof != fe.dofs_per_quad; ++dof)
                  dofs_on_mother[next_index++] =
                    this_face->dof_index(dof, fe_index);
                AssertDimension(next_index, dofs_on_mother.size());

                // TODO: assert some consistency assumptions

                // TODO[TL]: think about this in case of anisotropic refinement

                Assert(dof_handler.get_triangulation()
                           .get_anisotropic_refinement_flag() ||
                         ((this_face->child(0)->vertex_index(3) ==
                           this_face->child(1)->vertex_index(2)) &&
                          (this_face->child(0)->vertex_index(3) ==
                           this_face->child(2)->vertex_index(1)) &&
                          (this_face->child(0)->vertex_index(3) ==
                           this_face->child(3)->vertex_index(0))),
                       ExcInternalError());

                for (unsigned int dof = 0; dof != fe.dofs_per_vertex; ++dof)
                  dofs_on_children.push_back(
                    this_face->child(0)->vertex_dof_index(3, dof));

                // dof numbers on the centers of the lines bounding this face
                for (unsigned int line = 0; line < 4; ++line)
                  for (unsigned int dof = 0; dof != fe.dofs_per_vertex; ++dof)
                    dofs_on_children.push_back(
                      this_face->line(line)->child(0)->vertex_dof_index(
                        1, dof, fe_index));

                // next the dofs on the lines interior to the face; the order of
                // these lines is laid down in the FiniteElement class
                // documentation
                for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
                  dofs_on_children.push_back(
                    this_face->child(0)->line(1)->dof_index(dof, fe_index));
                for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
                  dofs_on_children.push_back(
                    this_face->child(2)->line(1)->dof_index(dof, fe_index));
                for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
                  dofs_on_children.push_back(
                    this_face->child(0)->line(3)->dof_index(dof, fe_index));
                for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
                  dofs_on_children.push_back(
                    this_face->child(1)->line(3)->dof_index(dof, fe_index));

                // dofs on the bordering lines
                for (unsigned int line = 0; line < 4; ++line)
                  for (unsigned int child = 0; child < 2; ++child)
                    {
                      for (unsigned int dof = 0; dof != fe.dofs_per_line; ++dof)
                        dofs_on_children.push_back(
                          this_face->line(line)->child(child)->dof_index(
                            dof, fe_index));
                    }

                // finally, for the dofs interior to the four child faces
                for (unsigned int child = 0; child < 4; ++child)
                  {
                    // skip artificial cells
                    if (cell->neighbor_child_on_subface(face, child)
                          ->is_artificial())
                      continue;
                    for (unsigned int dof = 0; dof != fe.dofs_per_quad; ++dof)
                      dofs_on_children.push_back(
                        this_face->child(child)->dof_index(dof, fe_index));
                  }

                // note: can get fewer DoFs when we have artificial cells:
                Assert(dofs_on_children.size() <= n_dofs_on_children,
                       ExcInternalError());

                // for each row in the AffineConstraints object for this line:
                for (unsigned int row = 0; row != dofs_on_children.size();
                     ++row)
                  {
                    constraints.add_line(dofs_on_children[row]);
                    for (unsigned int i = 0; i != dofs_on_mother.size(); ++i)
                      constraints.add_entry(dofs_on_children[row],
                                            dofs_on_mother[i],
                                            fe.constraints()(row, i));

                    constraints.set_inhomogeneity(dofs_on_children[row], 0.);
                  }
              }
            else
              {
                // this face has no children, but it could still be that it is
                // shared by two cells that use a different fe index. check a
                // couple of things, but ignore the case that the neighbor is an
                // artificial cell
                if (!cell->at_boundary(face) &&
                    !cell->neighbor(face)->is_artificial())
                  {
                    Assert(cell->face(face)->n_active_fe_indices() == 1,
                           ExcNotImplemented());
                    Assert(cell->face(face)->fe_index_is_active(
                             cell->active_fe_index()) == true,
                           ExcInternalError());
                  }
              }
        }
    }



    template <typename DoFHandlerType, typename number>
    void
    make_hp_hanging_node_constraints(const DoFHandlerType &     dof_handler,
                                     AffineConstraints<number> &constraints)
    {
      // note: this function is going to be hard to understand if you haven't
      // read the hp paper. however, we try to follow the notation laid out
      // there, so go read the paper before you try to understand what is going
      // on here

      const unsigned int dim = DoFHandlerType::dimension;

      const unsigned int spacedim = DoFHandlerType::space_dimension;


      // a matrix to be used for constraints below. declared here and simply
      // resized down below to avoid permanent re-allocation of memory
      FullMatrix<double> constraint_matrix;

      // similarly have arrays that will hold master and slave dof numbers, as
      // well as a scratch array needed for the complicated case below
      std::vector<types::global_dof_index> master_dofs;
      std::vector<types::global_dof_index> slave_dofs;
      std::vector<types::global_dof_index> scratch_dofs;

      // caches for the face and subface interpolation matrices between
      // different (or the same) finite elements. we compute them only once,
      // namely the first time they are needed, and then just reuse them
      Table<2, std::unique_ptr<FullMatrix<double>>> face_interpolation_matrices(
        n_finite_elements(dof_handler), n_finite_elements(dof_handler));
      Table<3, std::unique_ptr<FullMatrix<double>>>
        subface_interpolation_matrices(
          n_finite_elements(dof_handler),
          n_finite_elements(dof_handler),
          GeometryInfo<dim>::max_children_per_face);

      // similarly have a cache for the matrices that are split into their
      // master and slave parts, and for which the master part is inverted.
      // these two matrices are derived from the face interpolation matrix
      // as described in the @ref hp_paper "hp paper"
      Table<2,
            std::unique_ptr<std::pair<FullMatrix<double>, FullMatrix<double>>>>
        split_face_interpolation_matrices(n_finite_elements(dof_handler),
                                          n_finite_elements(dof_handler));

      // finally, for each pair of finite elements, have a mask that states
      // which of the degrees of freedom on the coarse side of a refined face
      // will act as master dofs.
      Table<2, std::unique_ptr<std::vector<bool>>> master_dof_masks(
        n_finite_elements(dof_handler), n_finite_elements(dof_handler));

      // loop over all faces
      //
      // note that even though we may visit a face twice if the neighboring
      // cells are equally refined, we can only visit each face with hanging
      // nodes once
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          // artificial cells can at best neighbor ghost cells, but we're not
          // interested in these interfaces
          if (cell->is_artificial())
            continue;

          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               ++face)
            if (cell->face(face)->has_children())
              {
                // first of all, make sure that we treat a case which is
                // possible, i.e. either no dofs on the face at all or no
                // anisotropic refinement
                if (cell->get_fe().dofs_per_face == 0)
                  continue;

                Assert(cell->face(face)->refinement_case() ==
                         RefinementCase<dim - 1>::isotropic_refinement,
                       ExcNotImplemented());

                // so now we've found a face of an active cell that has
                // children. that means that there are hanging nodes here.

                // in any case, faces can have at most two sets of active fe
                // indices, but here the face can have only one (namely the same
                // as that from the cell we're sitting on), and each of the
                // children can have only one as well. check this
                Assert(cell->face(face)->n_active_fe_indices() == 1,
                       ExcInternalError());
                Assert(cell->face(face)->fe_index_is_active(
                         cell->active_fe_index()) == true,
                       ExcInternalError());
                for (unsigned int c = 0; c < cell->face(face)->n_children();
                     ++c)
                  if (!cell->neighbor_child_on_subface(face, c)
                         ->is_artificial())
                    Assert(cell->face(face)->child(c)->n_active_fe_indices() ==
                             1,
                           ExcInternalError());

                // first find out whether we can constrain each of the subfaces
                // to the mother face. in the lingo of the hp paper, this would
                // be the simple case. note that we can short-circuit this
                // decision if the dof_handler doesn't support hp at all
                //
                // ignore all interfaces with artificial cells
                FiniteElementDomination::Domination mother_face_dominates =
                  FiniteElementDomination::either_element_can_dominate;

                // auxiliary variable which holds FE indices of the mother face
                // and its subfaces. This knowledge will be needed in hp-case
                // with neither_element_dominates.
                std::set<unsigned int> fe_ind_face_subface;
                fe_ind_face_subface.insert(cell->active_fe_index());

                if (DoFHandlerSupportsDifferentFEs<DoFHandlerType>::value ==
                    true)
                  for (unsigned int c = 0;
                       c < cell->face(face)->number_of_children();
                       ++c)
                    {
                      const auto subcell =
                        cell->neighbor_child_on_subface(face, c);
                      if (!subcell->is_artificial())
                        {
                          mother_face_dominates =
                            mother_face_dominates &
                            (cell->get_fe().compare_for_domination(
                              subcell->get_fe(), /*codim=*/1));
                          fe_ind_face_subface.insert(
                            subcell->active_fe_index());
                        }
                    }

                switch (mother_face_dominates)
                  {
                    case FiniteElementDomination::this_element_dominates:
                    case FiniteElementDomination::either_element_can_dominate:
                      {
                        // Case 1 (the simple case and the only case that can
                        // happen for non-hp DoFHandlers): The coarse element
                        // dominates the elements on the subfaces (or they are
                        // all the same)
                        //
                        // so we are going to constrain the DoFs on the face
                        // children against the DoFs on the face itself
                        master_dofs.resize(cell->get_fe().dofs_per_face);

                        cell->face(face)->get_dof_indices(
                          master_dofs, cell->active_fe_index());

                        // Now create constraints for the subfaces and
                        // assemble it. ignore all interfaces with artificial
                        // cells because we can only get to such interfaces if
                        // the current cell is a ghost cell
                        for (unsigned int c = 0;
                             c < cell->face(face)->n_children();
                             ++c)
                          {
                            if (cell->neighbor_child_on_subface(face, c)
                                  ->is_artificial())
                              continue;

                            const typename DoFHandlerType::active_face_iterator
                              subface = cell->face(face)->child(c);

                            Assert(subface->n_active_fe_indices() == 1,
                                   ExcInternalError());

                            const unsigned int subface_fe_index =
                              subface->nth_active_fe_index(0);

                            // we sometime run into the situation where for
                            // example on one big cell we have a FE_Q(1) and on
                            // the subfaces we have a mixture of FE_Q(1) and
                            // FE_Nothing. In that case, the face domination is
                            // either_element_can_dominate for the whole
                            // collection of subfaces, but on the particular
                            // subface between FE_Q(1) and FE_Nothing, there are
                            // no constraints that we need to take care of. in
                            // that case, just continue
                            if (cell->get_fe().compare_for_domination(
                                  subface->get_fe(subface_fe_index),
                                  /*codim=*/1) ==
                                FiniteElementDomination::no_requirements)
                              continue;

                            // Same procedure as for the mother cell. Extract
                            // the face DoFs from the cell DoFs.
                            slave_dofs.resize(
                              subface->get_fe(subface_fe_index).dofs_per_face);
                            subface->get_dof_indices(slave_dofs,
                                                     subface_fe_index);

                            for (const types::global_dof_index slave_dof :
                                 slave_dofs)
                              {
                                (void)slave_dof;
                                Assert(slave_dof != numbers::invalid_dof_index,
                                       ExcInternalError());
                              }

                            // Now create the element constraint for this
                            // subface.
                            //
                            // As a side remark, one may wonder the following:
                            // neighbor_child is clearly computed correctly,
                            // i.e. taking into account face_orientation (just
                            // look at the implementation of that function).
                            // however, we don't care about this here, when we
                            // ask for subface_interpolation on subface c. the
                            // question rather is: do we have to translate 'c'
                            // here as well?
                            //
                            // the answer is in fact 'no'. if one does that,
                            // results are wrong: constraints are added twice
                            // for the same pair of nodes but with differing
                            // weights. in addition, one can look at the
                            // deal.II/project_*_03 tests that look at exactly
                            // this case: there, we have a mesh with at least
                            // one face_orientation==false and hanging nodes,
                            // and the results of those tests show that the
                            // result of projection verifies the approximation
                            // properties of a finite element onto that mesh
                            ensure_existence_of_subface_matrix(
                              cell->get_fe(),
                              subface->get_fe(subface_fe_index),
                              c,
                              subface_interpolation_matrices
                                [cell->active_fe_index()][subface_fe_index][c]);

                            // Add constraints to global AffineConstraints
                            // object.
                            filter_constraints(master_dofs,
                                               slave_dofs,
                                               *(subface_interpolation_matrices
                                                   [cell->active_fe_index()]
                                                   [subface_fe_index][c]),
                                               constraints);
                          } // loop over subfaces

                        break;
                      } // Case 1

                    case FiniteElementDomination::other_element_dominates:
                    case FiniteElementDomination::neither_element_dominates:
                      {
                        // Case 2 (the "complex" case): at least one (the
                        // neither_... case) of the finer elements or all of
                        // them (the other_... case) is dominating. See the hp
                        // paper for a way how to deal with this situation
                        //
                        // since this is something that can only happen for hp
                        // dof handlers, add a check here...
                        Assert(DoFHandlerSupportsDifferentFEs<
                                 DoFHandlerType>::value == true,
                               ExcInternalError());

                        const dealii::hp::FECollection<dim, spacedim>
                          &fe_collection = dof_handler.get_fe_collection();
                        // we first have to find the finite element that is able
                        // to generate a space that all the other ones can be
                        // constrained to. At this point we potentially have
                        // different scenarios:
                        //
                        // 1) sub-faces dominate mother face and there is a
                        // dominating FE among sub faces. We could loop over sub
                        // faces to find the needed FE index. However, this will
                        // not work in the case when ...
                        //
                        // 2) there is no dominating FE among sub faces (e.g.
                        // Q1xQ2 vs Q2xQ1), but subfaces still dominate mother
                        // face (e.g. Q2xQ2). To cover this case we would have
                        // to find the least dominating element amongst all
                        // finite elements on sub faces.
                        //
                        // 3) Finally, it could happen that we got here because
                        // neither_element_dominates (e.g. Q1xQ1xQ2 and Q1xQ2xQ1
                        // for subfaces and Q2xQ1xQ1 for mother face). This
                        // requires finding the least dominating element amongst
                        // all finite elements on sub faces and the mother face.
                        //
                        // Note that the last solution covers the first two
                        // scenarios, thus we stick with it assuming that we
                        // won't lose much time/efficiency.
                        const unsigned int dominating_fe_index =
                          fe_collection.find_dominating_fe_extended(
                            fe_ind_face_subface, /*codim=*/1);

                        AssertThrow(
                          dominating_fe_index != numbers::invalid_unsigned_int,
                          ExcMessage(
                            "Could not find a least face dominating FE."));

                        const FiniteElement<dim, spacedim> &dominating_fe =
                          dof_handler.get_fe(dominating_fe_index);

                        // first get the interpolation matrix from the mother to
                        // the virtual dofs
                        Assert(dominating_fe.dofs_per_face <=
                                 cell->get_fe().dofs_per_face,
                               ExcInternalError());

                        ensure_existence_of_face_matrix(
                          dominating_fe,
                          cell->get_fe(),
                          face_interpolation_matrices[dominating_fe_index]
                                                     [cell->active_fe_index()]);

                        // split this matrix into master and slave components.
                        // invert the master component
                        ensure_existence_of_master_dof_mask(
                          cell->get_fe(),
                          dominating_fe,
                          (*face_interpolation_matrices
                             [dominating_fe_index][cell->active_fe_index()]),
                          master_dof_masks[dominating_fe_index]
                                          [cell->active_fe_index()]);

                        ensure_existence_of_split_face_matrix(
                          *face_interpolation_matrices[dominating_fe_index]
                                                      [cell->active_fe_index()],
                          (*master_dof_masks[dominating_fe_index]
                                            [cell->active_fe_index()]),
                          split_face_interpolation_matrices
                            [dominating_fe_index][cell->active_fe_index()]);

                        const FullMatrix<double>
                          &restrict_mother_to_virtual_master_inv =
                            (split_face_interpolation_matrices
                               [dominating_fe_index][cell->active_fe_index()]
                                 ->first);

                        const FullMatrix<double>
                          &restrict_mother_to_virtual_slave =
                            (split_face_interpolation_matrices
                               [dominating_fe_index][cell->active_fe_index()]
                                 ->second);

                        // now compute the constraint matrix as the product
                        // between the inverse matrix and the slave part
                        constraint_matrix.reinit(cell->get_fe().dofs_per_face -
                                                   dominating_fe.dofs_per_face,
                                                 dominating_fe.dofs_per_face);
                        restrict_mother_to_virtual_slave.mmult(
                          constraint_matrix,
                          restrict_mother_to_virtual_master_inv);

                        // then figure out the global numbers of master and
                        // slave dofs and apply constraints
                        scratch_dofs.resize(cell->get_fe().dofs_per_face);
                        cell->face(face)->get_dof_indices(
                          scratch_dofs, cell->active_fe_index());

                        // split dofs into master and slave components
                        master_dofs.clear();
                        slave_dofs.clear();
                        for (unsigned int i = 0;
                             i < cell->get_fe().dofs_per_face;
                             ++i)
                          if ((*master_dof_masks[dominating_fe_index]
                                                [cell->active_fe_index()])[i] ==
                              true)
                            master_dofs.push_back(scratch_dofs[i]);
                          else
                            slave_dofs.push_back(scratch_dofs[i]);

                        AssertDimension(master_dofs.size(),
                                        dominating_fe.dofs_per_face);
                        AssertDimension(slave_dofs.size(),
                                        cell->get_fe().dofs_per_face -
                                          dominating_fe.dofs_per_face);

                        filter_constraints(master_dofs,
                                           slave_dofs,
                                           constraint_matrix,
                                           constraints);



                        // next we have to deal with the subfaces. do as
                        // discussed in the hp paper
                        for (unsigned int sf = 0;
                             sf < cell->face(face)->n_children();
                             ++sf)
                          {
                            // ignore interfaces with artificial cells as well
                            // as interfaces between ghost cells in 2d
                            if (cell->neighbor_child_on_subface(face, sf)
                                  ->is_artificial() ||
                                (dim == 2 && cell->is_ghost() &&
                                 cell->neighbor_child_on_subface(face, sf)
                                   ->is_ghost()))
                              continue;

                            Assert(cell->face(face)
                                       ->child(sf)
                                       ->n_active_fe_indices() == 1,
                                   ExcInternalError());

                            const unsigned int subface_fe_index =
                              cell->face(face)->child(sf)->nth_active_fe_index(
                                0);
                            const FiniteElement<dim, spacedim> &subface_fe =
                              dof_handler.get_fe(subface_fe_index);

                            // first get the interpolation matrix from the
                            // subface to the virtual dofs
                            Assert(dominating_fe.dofs_per_face <=
                                     subface_fe.dofs_per_face,
                                   ExcInternalError());
                            ensure_existence_of_subface_matrix(
                              dominating_fe,
                              subface_fe,
                              sf,
                              subface_interpolation_matrices
                                [dominating_fe_index][subface_fe_index][sf]);

                            const FullMatrix<double>
                              &restrict_subface_to_virtual = *(
                                subface_interpolation_matrices
                                  [dominating_fe_index][subface_fe_index][sf]);

                            constraint_matrix.reinit(
                              subface_fe.dofs_per_face,
                              dominating_fe.dofs_per_face);

                            restrict_subface_to_virtual.mmult(
                              constraint_matrix,
                              restrict_mother_to_virtual_master_inv);

                            slave_dofs.resize(subface_fe.dofs_per_face);
                            cell->face(face)->child(sf)->get_dof_indices(
                              slave_dofs, subface_fe_index);

                            filter_constraints(master_dofs,
                                               slave_dofs,
                                               constraint_matrix,
                                               constraints);
                          } // loop over subfaces

                        break;
                      } // Case 2

                    case FiniteElementDomination::no_requirements:
                      // there are no continuity requirements between the two
                      // elements. record no constraints
                      break;

                    default:
                      // we shouldn't get here
                      Assert(false, ExcInternalError());
                  }
              }
            else
              {
                // this face has no children, but it could still be that it is
                // shared by two cells that use a different fe index
                Assert(cell->face(face)->fe_index_is_active(
                         cell->active_fe_index()) == true,
                       ExcInternalError());

                // see if there is a neighbor that is an artificial cell. in
                // that case, we're not interested in this interface. we test
                // this case first since artificial cells may not have an
                // active_fe_index set, etc
                if (!cell->at_boundary(face) &&
                    cell->neighbor(face)->is_artificial())
                  continue;

                // Only if there is a neighbor with a different active_fe_index
                // and the same h-level, some action has to be taken.
                if ((DoFHandlerSupportsDifferentFEs<DoFHandlerType>::value ==
                     true) &&
                    !cell->face(face)->at_boundary() &&
                    (cell->neighbor(face)->active_fe_index() !=
                     cell->active_fe_index()) &&
                    (!cell->face(face)->has_children() &&
                     !cell->neighbor_is_coarser(face)))
                  {
                    const typename DoFHandlerType::level_cell_iterator
                      neighbor = cell->neighbor(face);

                    // see which side of the face we have to constrain
                    switch (
                      cell->get_fe().compare_for_domination(neighbor->get_fe(),
                                                            /*codim=*/1))
                      {
                        case FiniteElementDomination::this_element_dominates:
                          {
                            // Get DoFs on dominating and dominated side of the
                            // face
                            master_dofs.resize(cell->get_fe().dofs_per_face);
                            cell->face(face)->get_dof_indices(
                              master_dofs, cell->active_fe_index());

                            // break if the n_master_dofs == 0, because we are
                            // attempting to constrain to an element that has no
                            // face dofs
                            if (master_dofs.size() == 0)
                              break;

                            slave_dofs.resize(neighbor->get_fe().dofs_per_face);
                            cell->face(face)->get_dof_indices(
                              slave_dofs, neighbor->active_fe_index());

                            // make sure the element constraints for this face
                            // are available
                            ensure_existence_of_face_matrix(
                              cell->get_fe(),
                              neighbor->get_fe(),
                              face_interpolation_matrices
                                [cell->active_fe_index()]
                                [neighbor->active_fe_index()]);

                            // Add constraints to global constraint matrix.
                            filter_constraints(
                              master_dofs,
                              slave_dofs,
                              *(face_interpolation_matrices
                                  [cell->active_fe_index()]
                                  [neighbor->active_fe_index()]),
                              constraints);

                            break;
                          }

                        case FiniteElementDomination::other_element_dominates:
                          {
                            // we don't do anything here since we will come back
                            // to this face from the other cell, at which time
                            // we will fall into the first case clause above
                            break;
                          }

                        case FiniteElementDomination::
                          either_element_can_dominate:
                          {
                            // it appears as if neither element has any
                            // constraints on its neighbor. this may be because
                            // neither element has any DoFs on faces at all. or
                            // that the two elements are actually the same,
                            // although they happen to run under different
                            // fe_indices (this is what happens in
                            // hp/hp_hanging_nodes_01 for example).
                            //
                            // another possibility is what happens in crash_13.
                            // there, we have FESystem(FE_Q(1),FE_DGQ(0)) vs.
                            // FESystem(FE_Q(1),FE_DGQ(1)). neither of them
                            // dominates the other.
                            //
                            // a final possibility is that we have something
                            // like FESystem(FE_Q(1),FE_Q(1)) vs
                            // FESystem(FE_Q(1),FE_Nothing()), see
                            // hp/fe_nothing_18/19.
                            //
                            // in any case, the point is that it doesn't matter.
                            // there is nothing to do here.
                            break;
                          }

                        case FiniteElementDomination::neither_element_dominates:
                          {
                            // make sure we don't get here twice from each cell
                            if (cell < neighbor)
                              break;

                            // our best bet is to find the common space among
                            // other FEs in FECollection and then constrain both
                            // FEs to that one. More precisely, we follow the
                            // strategy outlined on page 17 of the hp paper:
                            // First we find the dominant FE space S. Then we
                            // divide our dofs in master and slave such that
                            // I^{face,master}_{S^{face}->S} is invertible. And
                            // finally constrain slave dofs to master dofs based
                            // on the interpolation matrix.

                            const unsigned int this_fe_index =
                              cell->active_fe_index();
                            const unsigned int neighbor_fe_index =
                              neighbor->active_fe_index();
                            std::set<unsigned int> fes;
                            fes.insert(this_fe_index);
                            fes.insert(neighbor_fe_index);
                            const dealii::hp::FECollection<dim, spacedim>
                              &fe_collection = dof_handler.get_fe_collection();

                            const unsigned int dominating_fe_index =
                              fe_collection.find_dominating_fe_extended(
                                fes, /*codim=*/1);

                            AssertThrow(
                              dominating_fe_index !=
                                numbers::invalid_unsigned_int,
                              ExcMessage(
                                "Could not find the dominating FE for " +
                                cell->get_fe().get_name() + " and " +
                                neighbor->get_fe().get_name() +
                                " inside FECollection."));

                            const FiniteElement<dim, spacedim> &dominating_fe =
                              fe_collection[dominating_fe_index];

                            // TODO: until we hit the second face, the code is a
                            // copy-paste from h-refinement case...

                            // first get the interpolation matrix from main FE
                            // to the virtual dofs
                            Assert(dominating_fe.dofs_per_face <=
                                     cell->get_fe().dofs_per_face,
                                   ExcInternalError());

                            ensure_existence_of_face_matrix(
                              dominating_fe,
                              cell->get_fe(),
                              face_interpolation_matrices
                                [dominating_fe_index][cell->active_fe_index()]);

                            // split this matrix into master and slave
                            // components. invert the master component
                            ensure_existence_of_master_dof_mask(
                              cell->get_fe(),
                              dominating_fe,
                              (*face_interpolation_matrices
                                 [dominating_fe_index]
                                 [cell->active_fe_index()]),
                              master_dof_masks[dominating_fe_index]
                                              [cell->active_fe_index()]);

                            ensure_existence_of_split_face_matrix(
                              *face_interpolation_matrices
                                [dominating_fe_index][cell->active_fe_index()],
                              (*master_dof_masks[dominating_fe_index]
                                                [cell->active_fe_index()]),
                              split_face_interpolation_matrices
                                [dominating_fe_index][cell->active_fe_index()]);

                            const FullMatrix<
                              double> &restrict_mother_to_virtual_master_inv =
                              (split_face_interpolation_matrices
                                 [dominating_fe_index][cell->active_fe_index()]
                                   ->first);

                            const FullMatrix<
                              double> &restrict_mother_to_virtual_slave =
                              (split_face_interpolation_matrices
                                 [dominating_fe_index][cell->active_fe_index()]
                                   ->second);

                            // now compute the constraint matrix as the product
                            // between the inverse matrix and the slave part
                            constraint_matrix.reinit(
                              cell->get_fe().dofs_per_face -
                                dominating_fe.dofs_per_face,
                              dominating_fe.dofs_per_face);
                            restrict_mother_to_virtual_slave.mmult(
                              constraint_matrix,
                              restrict_mother_to_virtual_master_inv);

                            // then figure out the global numbers of master and
                            // slave dofs and apply constraints
                            scratch_dofs.resize(cell->get_fe().dofs_per_face);
                            cell->face(face)->get_dof_indices(
                              scratch_dofs, cell->active_fe_index());

                            // split dofs into master and slave components
                            master_dofs.clear();
                            slave_dofs.clear();
                            for (unsigned int i = 0;
                                 i < cell->get_fe().dofs_per_face;
                                 ++i)
                              if ((*master_dof_masks[dominating_fe_index]
                                                    [cell->active_fe_index()])
                                    [i] == true)
                                master_dofs.push_back(scratch_dofs[i]);
                              else
                                slave_dofs.push_back(scratch_dofs[i]);

                            AssertDimension(master_dofs.size(),
                                            dominating_fe.dofs_per_face);
                            AssertDimension(slave_dofs.size(),
                                            cell->get_fe().dofs_per_face -
                                              dominating_fe.dofs_per_face);

                            filter_constraints(master_dofs,
                                               slave_dofs,
                                               constraint_matrix,
                                               constraints);

                            // now do the same for another FE this is pretty
                            // much the same we do above to resolve h-refinement
                            // constraints
                            Assert(dominating_fe.dofs_per_face <=
                                     neighbor->get_fe().dofs_per_face,
                                   ExcInternalError());

                            ensure_existence_of_face_matrix(
                              dominating_fe,
                              neighbor->get_fe(),
                              face_interpolation_matrices
                                [dominating_fe_index]
                                [neighbor->active_fe_index()]);

                            const FullMatrix<double>
                              &restrict_secondface_to_virtual =
                                *(face_interpolation_matrices
                                    [dominating_fe_index]
                                    [neighbor->active_fe_index()]);

                            constraint_matrix.reinit(
                              neighbor->get_fe().dofs_per_face,
                              dominating_fe.dofs_per_face);

                            restrict_secondface_to_virtual.mmult(
                              constraint_matrix,
                              restrict_mother_to_virtual_master_inv);

                            slave_dofs.resize(neighbor->get_fe().dofs_per_face);
                            cell->face(face)->get_dof_indices(
                              slave_dofs, neighbor->active_fe_index());

                            filter_constraints(master_dofs,
                                               slave_dofs,
                                               constraint_matrix,
                                               constraints);

                            break;
                          }

                        case FiniteElementDomination::no_requirements:
                          {
                            // nothing to do here
                            break;
                          }

                        default:
                          // we shouldn't get here
                          Assert(false, ExcInternalError());
                      }
                  }
              }
        }
    }
  } // namespace internal



  template <typename DoFHandlerType, typename number>
  void
  make_hanging_node_constraints(const DoFHandlerType &     dof_handler,
                                AffineConstraints<number> &constraints)
  {
    // Decide whether to use the new or old make_hanging_node_constraints
    // function. If all the FiniteElement or all elements in a FECollection
    // support the new face constraint matrix, the new code will be used.
    // Otherwise, the old implementation is used for the moment.
    if (dof_handler.get_fe_collection().hp_constraints_are_implemented())
      internal::make_hp_hanging_node_constraints(dof_handler, constraints);
    else
      internal::make_oldstyle_hanging_node_constraints(
        dof_handler,
        constraints,
        std::integral_constant<int, DoFHandlerType::dimension>());
  }



  namespace
  {
    /**
     * @internal
     *
     * Internally used in make_periodicity_constraints.
     *
     * enter constraints for periodicity into the given AffineConstraints
     * object. this function is called when at least one of the two face
     * iterators corresponds to an active object without further children
     *
     * @param transformation A matrix that maps degrees of freedom from one face
     * to another. If the DoFs on the two faces are supposed to match exactly,
     * then the matrix so provided will be the identity matrix. if face 2 is
     * once refined from face 1, then the matrix needs to be the interpolation
     * matrix from a face to this particular child
     *
     * @precondition: face_1 is supposed to be active
     *
     * @note We have to be careful not to accidentally create constraint
     * cycles when adding periodic constraints: For example, as the
     * corresponding testcase bits/periodicity_05 demonstrates, we can
     * occasionally get into trouble if we already have the constraint
     * x1 == x2 and want to insert x2 == x1. We avoid this by skipping
     * such "identity constraints" if the opposite constraint already
     * exists.
     */
    template <typename FaceIterator, typename number>
    void
    set_periodicity_constraints(
      const FaceIterator &                         face_1,
      const typename identity<FaceIterator>::type &face_2,
      const FullMatrix<double> &                   transformation,
      AffineConstraints<number> &                  affine_constraints,
      const ComponentMask &                        component_mask,
      const bool                                   face_orientation,
      const bool                                   face_flip,
      const bool                                   face_rotation)
    {
      static const int dim      = FaceIterator::AccessorType::dimension;
      static const int spacedim = FaceIterator::AccessorType::space_dimension;

      // we should be in the case where face_1 is active, i.e. has no children:
      Assert(!face_1->has_children(), ExcInternalError());

      Assert(face_1->n_active_fe_indices() == 1, ExcInternalError());

      // If face_2 does have children, then we need to iterate over these
      // children and set periodic constraints in the inverse direction:

      if (face_2->has_children())
        {
          Assert(face_2->n_children() ==
                   GeometryInfo<dim>::max_children_per_face,
                 ExcNotImplemented());

          const unsigned int dofs_per_face =
            face_1->get_fe(face_1->nth_active_fe_index(0)).dofs_per_face;
          FullMatrix<double> child_transformation(dofs_per_face, dofs_per_face);
          FullMatrix<double> subface_interpolation(dofs_per_face,
                                                   dofs_per_face);

          for (unsigned int c = 0; c < face_2->n_children(); ++c)
            {
              // get the interpolation matrix recursively from the one that
              // interpolated from face_1 to face_2 by multiplying from the left
              // with the one that interpolates from face_2 to its child
              const auto &fe = face_1->get_fe(face_1->nth_active_fe_index(0));
              fe.get_subface_interpolation_matrix(fe, c, subface_interpolation);
              subface_interpolation.mmult(child_transformation, transformation);

              set_periodicity_constraints(face_1,
                                          face_2->child(c),
                                          child_transformation,
                                          affine_constraints,
                                          component_mask,
                                          face_orientation,
                                          face_flip,
                                          face_rotation);
            }
          return;
        }

      //
      // If we reached this point then both faces are active. Now all
      // that is left is to match the corresponding DoFs of both faces.
      //

      const unsigned int face_1_index = face_1->nth_active_fe_index(0);
      const unsigned int face_2_index = face_2->nth_active_fe_index(0);
      Assert(face_1->get_fe(face_1_index) == face_2->get_fe(face_2_index),
             ExcMessage(
               "Matching periodic cells need to use the same finite element"));

      const FiniteElement<dim, spacedim> &fe = face_1->get_fe(face_1_index);

      Assert(component_mask.represents_n_components(fe.n_components()),
             ExcMessage(
               "The number of components in the mask has to be either "
               "zero or equal to the number of components in the finite "
               "element."));

      const unsigned int dofs_per_face = fe.dofs_per_face;

      std::vector<types::global_dof_index> dofs_1(dofs_per_face);
      std::vector<types::global_dof_index> dofs_2(dofs_per_face);

      face_1->get_dof_indices(dofs_1, face_1_index);
      face_2->get_dof_indices(dofs_2, face_2_index);

      // If either of the two faces has an invalid dof index, stop. This is
      // so that there is no attempt to match artificial cells of parallel
      // distributed triangulations.
      //
      // While it seems like we ought to be able to avoid even calling
      // set_periodicity_constraints for artificial faces, this situation
      // can arise when a face that is being made periodic is only
      // partially touched by the local subdomain.
      // make_periodicity_constraints will be called recursively even for
      // the section of the face that is not touched by the local
      // subdomain.
      //
      // Until there is a better way to determine if the cells that
      // neighbor a face are artificial, we simply test to see if the face
      // does not have a valid dof initialization.

      for (unsigned int i = 0; i < dofs_per_face; i++)
        if (dofs_1[i] == numbers::invalid_dof_index ||
            dofs_2[i] == numbers::invalid_dof_index)
          {
            return;
          }

      // Well, this is a hack:
      //
      // There is no
      //   face_to_face_index(face_index,
      //                      face_orientation,
      //                      face_flip,
      //                      face_rotation)
      // function in FiniteElementData, so we have to use
      //   face_to_cell_index(face_index, face
      //                      face_orientation,
      //                      face_flip,
      //                      face_rotation)
      // But this will give us an index on a cell - something we cannot work
      // with directly. But luckily we can match them back :-]

      std::map<unsigned int, unsigned int> cell_to_rotated_face_index;

      // Build up a cell to face index for face_2:
      for (unsigned int i = 0; i < dofs_per_face; ++i)
        {
          const unsigned int cell_index =
            fe.face_to_cell_index(i,
                                  0, /* It doesn't really matter, just
                                      * assume we're on the first face...
                                      */
                                  true,
                                  false,
                                  false // default orientation
            );
          cell_to_rotated_face_index[cell_index] = i;
        }

      //
      // Loop over all dofs on face 2 and constrain them against all
      // matching dofs on face 1:
      //

      for (unsigned int i = 0; i < dofs_per_face; ++i)
        {
          // Obey the component mask
          if ((component_mask.n_selected_components(fe.n_components()) !=
               fe.n_components()) &&
              !component_mask[fe.face_system_to_component_index(i).first])
            continue;

          // We have to be careful to treat so called "identity
          // constraints" special. These are constraints of the form
          // x1 == factor * x_2. In this case, if the constraint
          // x2 == 1./factor * x1 already exists we are in trouble.
          //
          // Consequently, we have to check that we have indeed such an
          // "identity constraint". We do this by looping over all entries
          // of the row of the transformation matrix and check whether we
          // find exactly one nonzero entry. If this is the case, set
          // "is_identity_constrained" to true and record the corresponding
          // index and factor.

          bool         is_identity_constrained = false;
          unsigned int target                  = numbers::invalid_unsigned_int;
          double       factor                  = 1.;

          constexpr double eps = 1.e-13;
          for (unsigned int jj = 0; jj < dofs_per_face; ++jj)
            {
              const auto entry = transformation(i, jj);
              if (std::abs(entry) > eps)
                {
                  if (is_identity_constrained)
                    {
                      // We did encounter more than one nonzero entry, so
                      // the dof is not identity constrained. Set the
                      // boolean to false and break out of the for loop.
                      is_identity_constrained = false;
                      break;
                    }
                  is_identity_constrained = true;
                  target                  = jj;
                  factor                  = entry;
                }
            }

          // Next, we work on all constraints that are not identity
          // constraints, i.e., constraints that involve an interpolation
          // step that constrains the current dof (on face 2) to more than
          // one dof on face 1.

          if (!is_identity_constrained)
            {
              // The current dof is already constrained. There is nothing
              // left to do.
              if (affine_constraints.is_constrained(dofs_2[i]))
                continue;

              // Enter the constraint piece by piece:
              affine_constraints.add_line(dofs_2[i]);

              for (unsigned int jj = 0; jj < dofs_per_face; ++jj)
                {
                  // Get the correct dof index on face_1 respecting the
                  // given orientation:
                  const unsigned int j =
                    cell_to_rotated_face_index[fe.face_to_cell_index(
                      jj, 0, face_orientation, face_flip, face_rotation)];

                  if (std::abs(transformation(i, jj)) > eps)
                    affine_constraints.add_entry(dofs_2[i],
                                                 dofs_1[j],
                                                 transformation(i, jj));
                }

              // Continue with next dof.
              continue;
            }

          // We are left with an "identity constraint".

          // Get the correct dof index on face_1 respecting the given
          // orientation:
          const unsigned int j =
            cell_to_rotated_face_index[fe.face_to_cell_index(
              target, 0, face_orientation, face_flip, face_rotation)];

          auto dof_left  = dofs_1[j];
          auto dof_right = dofs_2[i];

          // If dof_left is already constrained, or dof_left < dof_right we
          // flip the order to ensure that dofs are constrained in a stable
          // manner on different MPI processes.
          if (affine_constraints.is_constrained(dof_left) ||
              (dof_left < dof_right &&
               !affine_constraints.is_constrained(dof_right)))
            {
              std::swap(dof_left, dof_right);
              factor = 1. / factor;
            }

          // Next, we try to enter the constraint
          //   dof_left = factor * dof_right;

          // If both degrees of freedom are constrained, there is nothing we
          // can do. Simply continue with the next dof.
          if (affine_constraints.is_constrained(dof_left) &&
              affine_constraints.is_constrained(dof_right))
            continue;

          // We have to be careful that adding the current identity
          // constraint does not create a constraint cycle. Thus, check for
          // a dependency cycle:

          bool   constraints_are_cyclic = true;
          number cycle_factor           = factor;

          for (auto test_dof = dof_right; test_dof != dof_left;)
            {
              if (!affine_constraints.is_constrained(test_dof))
                {
                  constraints_are_cyclic = false;
                  break;
                }

              const auto &constraint_entries =
                *affine_constraints.get_constraint_entries(test_dof);
              if (constraint_entries.size() == 1)
                {
                  test_dof = constraint_entries[0].first;
                  cycle_factor *= constraint_entries[0].second;
                }
              else
                {
                  constraints_are_cyclic = false;
                  break;
                }
            }

          // In case of a dependency cycle we, either
          //  - do nothing if cycle_factor == 1. In this case all degrees
          //    of freedom are already periodically constrained,
          //  - otherwise, force all dofs to zero (by setting dof_left to
          //    zero). The reasoning behind this is the fact that
          //    cycle_factor != 1 occurs in situations such as
          //    x1 == x2 and x2 == -1. * x1. This system is only solved by
          //    x_1 = x_2 = 0.

          if (constraints_are_cyclic)
            {
              if (std::abs(cycle_factor - 1.) > eps)
                affine_constraints.add_line(dof_left);
            }
          else
            {
              affine_constraints.add_line(dof_left);
              affine_constraints.add_entry(dof_left, dof_right, factor);
            }
        } /* for dofs_per_face */
    }



    // Internally used in make_periodicity_constraints.
    //
    // Build up a (possibly rotated) interpolation matrix that is used in
    // set_periodicity_constraints with the help of user supplied matrix and
    // first_vector_components.
    template <int dim, int spacedim>
    FullMatrix<double>
    compute_transformation(
      const FiniteElement<dim, spacedim> &fe,
      const FullMatrix<double> &          matrix,
      const std::vector<unsigned int> &   first_vector_components)
    {
      Assert(matrix.m() == matrix.n(), ExcInternalError());

      const unsigned int n_dofs_per_face = fe.dofs_per_face;

      if (matrix.m() == n_dofs_per_face)
        {
          // In case of m == n == n_dofs_per_face the supplied matrix is already
          // an interpolation matrix, so we use it directly:
          return matrix;
        }

      if (first_vector_components.empty() && matrix.m() == 0)
        {
          // Just the identity matrix in case no rotation is specified:
          return IdentityMatrix(n_dofs_per_face);
        }

      // The matrix describes a rotation and we have to build a transformation
      // matrix, we assume that for a 0* rotation we would have to build the
      // identity matrix

      Assert(matrix.m() == spacedim, ExcInternalError())

        Quadrature<dim - 1>
          quadrature(fe.get_unit_face_support_points());

      // have an array that stores the location of each vector-dof tuple we want
      // to rotate.
      using DoFTuple = std::array<unsigned int, spacedim>;

      // start with a pristine interpolation matrix...
      FullMatrix<double> transformation = IdentityMatrix(n_dofs_per_face);

      for (unsigned int i = 0; i < n_dofs_per_face; ++i)
        {
          std::vector<unsigned int>::const_iterator comp_it =
            std::find(first_vector_components.begin(),
                      first_vector_components.end(),
                      fe.face_system_to_component_index(i).first);
          if (comp_it != first_vector_components.end())
            {
              const unsigned int first_vector_component = *comp_it;

              // find corresponding other components of vector
              DoFTuple vector_dofs;
              vector_dofs[0]       = i;
              unsigned int n_found = 1;

              Assert(
                *comp_it + spacedim <= fe.n_components(),
                ExcMessage(
                  "Error: the finite element does not have enough components "
                  "to define rotated periodic boundaries."));

              for (unsigned int k = 0; k < n_dofs_per_face; ++k)
                if ((k != i) && (quadrature.point(k) == quadrature.point(i)) &&
                    (fe.face_system_to_component_index(k).first >=
                     first_vector_component) &&
                    (fe.face_system_to_component_index(k).first <
                     first_vector_component + spacedim))
                  {
                    vector_dofs[fe.face_system_to_component_index(k).first -
                                first_vector_component] = k;
                    n_found++;
                    if (n_found == dim)
                      break;
                  }

              // ... and rotate all dofs belonging to vector valued components
              // that are selected by first_vector_components:
              for (int i = 0; i < spacedim; ++i)
                {
                  transformation[vector_dofs[i]][vector_dofs[i]] = 0.;
                  for (int j = 0; j < spacedim; ++j)
                    transformation[vector_dofs[i]][vector_dofs[j]] =
                      matrix[i][j];
                }
            }
        }
      return transformation;
    }
  } /*namespace*/


  // Low level interface:


  template <typename FaceIterator, typename number>
  void
  make_periodicity_constraints(
    const FaceIterator &                         face_1,
    const typename identity<FaceIterator>::type &face_2,
    AffineConstraints<number> &                  affine_constraints,
    const ComponentMask &                        component_mask,
    const bool                                   face_orientation,
    const bool                                   face_flip,
    const bool                                   face_rotation,
    const FullMatrix<double> &                   matrix,
    const std::vector<unsigned int> &            first_vector_components)
  {
    static const int dim      = FaceIterator::AccessorType::dimension;
    static const int spacedim = FaceIterator::AccessorType::space_dimension;

    Assert((dim != 1) || (face_orientation == true && face_flip == false &&
                          face_rotation == false),
           ExcMessage("The supplied orientation "
                      "(face_orientation, face_flip, face_rotation) "
                      "is invalid for 1D"));

    Assert((dim != 2) || (face_orientation == true && face_rotation == false),
           ExcMessage("The supplied orientation "
                      "(face_orientation, face_flip, face_rotation) "
                      "is invalid for 2D"));

    Assert(face_1 != face_2,
           ExcMessage("face_1 and face_2 are equal! Cannot constrain DoFs "
                      "on the very same face"));

    Assert(face_1->at_boundary() && face_2->at_boundary(),
           ExcMessage("Faces for periodicity constraints must be on the "
                      "boundary"));

    Assert(matrix.m() == matrix.n(),
           ExcMessage("The supplied (rotation or interpolation) matrix must "
                      "be a square matrix"));

    Assert(first_vector_components.empty() || matrix.m() == spacedim,
           ExcMessage("first_vector_components is nonempty, so matrix must "
                      "be a rotation matrix exactly of size spacedim"));

#ifdef DEBUG
    if (!face_1->has_children())
      {
        Assert(face_1->n_active_fe_indices() == 1, ExcInternalError());
        const unsigned int n_dofs_per_face =
          face_1->get_fe(face_1->nth_active_fe_index(0)).dofs_per_face;

        Assert(matrix.m() == 0 ||
                 (first_vector_components.empty() &&
                  matrix.m() == n_dofs_per_face) ||
                 (!first_vector_components.empty() && matrix.m() == spacedim),
               ExcMessage(
                 "The matrix must have either size 0 or spacedim "
                 "(if first_vector_components is nonempty) "
                 "or the size must be equal to the # of DoFs on the face "
                 "(if first_vector_components is empty)."));
      }

    if (!face_2->has_children())
      {
        Assert(face_2->n_active_fe_indices() == 1, ExcInternalError());
        const unsigned int n_dofs_per_face =
          face_2->get_fe(face_2->nth_active_fe_index(0)).dofs_per_face;

        Assert(matrix.m() == 0 ||
                 (first_vector_components.empty() &&
                  matrix.m() == n_dofs_per_face) ||
                 (!first_vector_components.empty() && matrix.m() == spacedim),
               ExcMessage(
                 "The matrix must have either size 0 or spacedim "
                 "(if first_vector_components is nonempty) "
                 "or the size must be equal to the # of DoFs on the face "
                 "(if first_vector_components is empty)."));
      }
#endif

    // A lookup table on how to go through the child faces depending on the
    // orientation:

    static const int lookup_table_2d[2][2] = {
      //          flip:
      {0, 1}, //  false
      {1, 0}, //  true
    };

    static const int lookup_table_3d[2][2][2][4] = {
      //                    orientation flip  rotation
      {
        {
          {0, 2, 1, 3}, //  false       false false
          {2, 3, 0, 1}, //  false       false true
        },
        {
          {3, 1, 2, 0}, //  false       true  false
          {1, 0, 3, 2}, //  false       true  true
        },
      },
      {
        {
          {0, 1, 2, 3}, //  true        false false
          {1, 3, 0, 2}, //  true        false true
        },
        {
          {3, 2, 1, 0}, //  true        true  false
          {2, 0, 3, 1}, //  true        true  true
        },
      },
    };

    if (face_1->has_children() && face_2->has_children())
      {
        // In the case that both faces have children, we loop over all children
        // and apply make_periodicty_constrains recursively:

        Assert(face_1->n_children() ==
                   GeometryInfo<dim>::max_children_per_face &&
                 face_2->n_children() ==
                   GeometryInfo<dim>::max_children_per_face,
               ExcNotImplemented());

        for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_face;
             ++i)
          {
            // Lookup the index for the second face
            unsigned int j;
            switch (dim)
              {
                case 2:
                  j = lookup_table_2d[face_flip][i];
                  break;
                case 3:
                  j = lookup_table_3d[face_orientation][face_flip]
                                     [face_rotation][i];
                  break;
                default:
                  AssertThrow(false, ExcNotImplemented());
              }

            make_periodicity_constraints(face_1->child(i),
                                         face_2->child(j),
                                         affine_constraints,
                                         component_mask,
                                         face_orientation,
                                         face_flip,
                                         face_rotation,
                                         matrix,
                                         first_vector_components);
          }
      }
    else
      {
        // Otherwise at least one of the two faces is active and we need to do
        // some work and enter the constraints!

        // The finite element that matters is the one on the active face:
        const FiniteElement<dim, spacedim> &fe =
          face_1->has_children() ?
            face_2->get_fe(face_2->nth_active_fe_index(0)) :
            face_1->get_fe(face_1->nth_active_fe_index(0));

        const unsigned int n_dofs_per_face = fe.dofs_per_face;

        // Sometimes we just have nothing to do (for all finite elements, or
        // systems which accidentally don't have any dofs on the boundary).
        if (n_dofs_per_face == 0)
          return;

        const FullMatrix<double> transformation =
          compute_transformation(fe, matrix, first_vector_components);

        if (!face_2->has_children())
          {
            // Performance hack: We do not need to compute an inverse if the
            // matrix is the identity matrix.
            if (first_vector_components.empty() && matrix.m() == 0)
              {
                set_periodicity_constraints(face_2,
                                            face_1,
                                            transformation,
                                            affine_constraints,
                                            component_mask,
                                            face_orientation,
                                            face_flip,
                                            face_rotation);
              }
            else
              {
                FullMatrix<double> inverse(transformation.m());
                inverse.invert(transformation);

                set_periodicity_constraints(face_2,
                                            face_1,
                                            inverse,
                                            affine_constraints,
                                            component_mask,
                                            face_orientation,
                                            face_flip,
                                            face_rotation);
              }
          }
        else
          {
            Assert(!face_1->has_children(), ExcInternalError());

            // Important note:
            // In 3D we have to take care of the fact that face_rotation gives
            // the relative rotation of face_1 to face_2, i.e. we have to invert
            // the rotation when constraining face_2 to face_1. Therefore
            // face_flip has to be toggled if face_rotation is true: In case of
            // inverted orientation, nothing has to be done.
            set_periodicity_constraints(face_1,
                                        face_2,
                                        transformation,
                                        affine_constraints,
                                        component_mask,
                                        face_orientation,
                                        face_orientation ?
                                          face_rotation ^ face_flip :
                                          face_flip,
                                        face_rotation);
          }
      }
  }



  template <typename DoFHandlerType, typename number>
  void
  make_periodicity_constraints(
    const std::vector<
      GridTools::PeriodicFacePair<typename DoFHandlerType::cell_iterator>>
      &                              periodic_faces,
    AffineConstraints<number> &      constraints,
    const ComponentMask &            component_mask,
    const std::vector<unsigned int> &first_vector_components)
  {
    // Loop over all periodic faces...
    for (auto &pair : periodic_faces)
      {
        using FaceIterator        = typename DoFHandlerType::face_iterator;
        const FaceIterator face_1 = pair.cell[0]->face(pair.face_idx[0]);
        const FaceIterator face_2 = pair.cell[1]->face(pair.face_idx[1]);

        Assert(face_1->at_boundary() && face_2->at_boundary(),
               ExcInternalError());

        Assert(face_1 != face_2, ExcInternalError());

        // ... and apply the low level make_periodicity_constraints function to
        // every matching pair:
        make_periodicity_constraints(face_1,
                                     face_2,
                                     constraints,
                                     component_mask,
                                     pair.orientation[0],
                                     pair.orientation[1],
                                     pair.orientation[2],
                                     pair.matrix,
                                     first_vector_components);
      }
  }


  // High level interface variants:


  template <typename DoFHandlerType, typename number>
  void
  make_periodicity_constraints(const DoFHandlerType &             dof_handler,
                               const types::boundary_id           b_id1,
                               const types::boundary_id           b_id2,
                               const unsigned int                 direction,
                               dealii::AffineConstraints<number> &constraints,
                               const ComponentMask &component_mask)
  {
    static const int space_dim = DoFHandlerType::space_dimension;
    (void)space_dim;
    Assert(0 <= direction && direction < space_dim,
           ExcIndexRange(direction, 0, space_dim));

    Assert(b_id1 != b_id2,
           ExcMessage("The boundary indicators b_id1 and b_id2 must be "
                      "different to denote different boundaries."));

    std::vector<
      GridTools::PeriodicFacePair<typename DoFHandlerType::cell_iterator>>
      matched_faces;

    // Collect matching periodic cells on the coarsest level:
    GridTools::collect_periodic_faces(
      dof_handler, b_id1, b_id2, direction, matched_faces);

    make_periodicity_constraints<DoFHandlerType>(matched_faces,
                                                 constraints,
                                                 component_mask);
  }



  template <typename DoFHandlerType, typename number>
  void
  make_periodicity_constraints(const DoFHandlerType &     dof_handler,
                               const types::boundary_id   b_id,
                               const unsigned int         direction,
                               AffineConstraints<number> &constraints,
                               const ComponentMask &      component_mask)
  {
    static const int dim       = DoFHandlerType::dimension;
    static const int space_dim = DoFHandlerType::space_dimension;
    (void)dim;
    (void)space_dim;

    Assert(0 <= direction && direction < space_dim,
           ExcIndexRange(direction, 0, space_dim));

    Assert(dim == space_dim, ExcNotImplemented());

    std::vector<
      GridTools::PeriodicFacePair<typename DoFHandlerType::cell_iterator>>
      matched_faces;

    // Collect matching periodic cells on the coarsest level:
    GridTools::collect_periodic_faces(dof_handler,
                                      b_id,
                                      direction,
                                      matched_faces);

    make_periodicity_constraints<DoFHandlerType>(matched_faces,
                                                 constraints,
                                                 component_mask);
  }



  namespace internal
  {
    namespace Assembler
    {
      struct Scratch
      {};


      template <int dim, int spacedim>
      struct CopyData
      {
        unsigned int                         dofs_per_cell;
        std::vector<types::global_dof_index> parameter_dof_indices;
#ifdef DEAL_II_WITH_MPI
        std::vector<dealii::LinearAlgebra::distributed::Vector<double>>
          global_parameter_representation;
#else
        std::vector<dealii::Vector<double>> global_parameter_representation;
#endif
      };
    } // namespace Assembler

    namespace
    {
      /**
       * This is a function that is called by the _2 function and that operates
       * on one cell only. It is worked in parallel if multhithreading is
       * available.
       */
      template <int dim, int spacedim>
      void
      compute_intergrid_weights_3(
        const typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator
          &cell,
        const Assembler::Scratch &,
        Assembler::CopyData<dim, spacedim> &copy_data,
        const unsigned int                  coarse_component,
        const FiniteElement<dim, spacedim> &coarse_fe,
        const InterGridMap<dealii::DoFHandler<dim, spacedim>>
          &                                        coarse_to_fine_grid_map,
        const std::vector<dealii::Vector<double>> &parameter_dofs)
      {
        // for each cell on the parameter grid: find out which degrees of
        // freedom on the fine grid correspond in which way to the degrees of
        // freedom on the parameter grid
        //
        // since for continuous FEs some dofs exist on more than one cell, we
        // have to track which ones were already visited. the problem is that if
        // we visit a dof first on one cell and compute its weight with respect
        // to some global dofs to be non-zero, and later visit the dof again on
        // another cell and (since we are on another cell) recompute the weights
        // with respect to the same dofs as above to be zero now, we have to
        // preserve them. we therefore overwrite all weights if they are nonzero
        // and do not enforce zero weights since that might be only due to the
        // fact that we are on another cell.
        //
        // example:
        // coarse grid
        //  |     |     |
        //  *-----*-----*
        //  | cell|cell |
        //  |  1  |  2  |
        //  |     |     |
        //  0-----1-----*
        //
        // fine grid
        //  |  |  |  |  |
        //  *--*--*--*--*
        //  |  |  |  |  |
        //  *--*--*--*--*
        //  |  |  |  |  |
        //  *--x--y--*--*
        //
        // when on cell 1, we compute the weights of dof 'x' to be 1/2 from
        // parameter dofs 0 and 1, respectively. however, when later we are on
        // cell 2, we again compute the prolongation of shape function 1
        // restricted to cell 2 to the globla grid and find that the weight of
        // global dof 'x' now is zero. however, we should not overwrite the old
        // value.
        //
        // we therefore always only set nonzero values. why adding up is not
        // useful: dof 'y' would get weight 1 from parameter dof 1 on both cells
        // 1 and 2, but the correct weight is nevertheless only 1.

        // vector to hold the representation of a single degree of freedom on
        // the coarse grid (for the selected fe) on the fine grid

        copy_data.dofs_per_cell = coarse_fe.dofs_per_cell;
        copy_data.parameter_dof_indices.resize(copy_data.dofs_per_cell);

        // get the global indices of the parameter dofs on this parameter grid
        // cell
        cell->get_dof_indices(copy_data.parameter_dof_indices);

        // loop over all dofs on this cell and check whether they are
        // interesting for us
        for (unsigned int local_dof = 0; local_dof < copy_data.dofs_per_cell;
             ++local_dof)
          if (coarse_fe.system_to_component_index(local_dof).first ==
              coarse_component)
            {
              // the how-many-th parameter is this on this cell?
              const unsigned int local_parameter_dof =
                coarse_fe.system_to_component_index(local_dof).second;

              copy_data.global_parameter_representation[local_parameter_dof] =
                0.;

              // distribute the representation of @p{local_parameter_dof} on the
              // parameter grid cell
              // @p{cell} to the global data space
              coarse_to_fine_grid_map[cell]->set_dof_values_by_interpolation(
                parameter_dofs[local_parameter_dof],
                copy_data.global_parameter_representation[local_parameter_dof]);
            }
      }



      /**
       * This is a function that is called by the _2 function and that operates
       * on one cell only. It is worked in parallel if multhithreading is
       * available.
       */
      template <int dim, int spacedim>
      void
      copy_intergrid_weights_3(
        const Assembler::CopyData<dim, spacedim> &  copy_data,
        const unsigned int                          coarse_component,
        const FiniteElement<dim, spacedim> &        coarse_fe,
        const std::vector<types::global_dof_index> &weight_mapping,
        const bool                                  is_called_in_parallel,
        std::vector<std::map<types::global_dof_index, float>> &weights)
      {
        unsigned int pos = 0;
        for (unsigned int local_dof = 0; local_dof < copy_data.dofs_per_cell;
             ++local_dof)
          if (coarse_fe.system_to_component_index(local_dof).first ==
              coarse_component)
            {
              // now that we've got the global representation of each parameter
              // dof, we've only got to clobber the non-zero entries in that
              // vector and store the result
              //
              // what we have learned: if entry @p{i} of the global vector holds
              // the value @p{v[i]}, then this is the weight with which the
              // present dof contributes to @p{i}. there may be several such
              // @p{i}s and their weights' sum should be one. Then, @p{v[i]}
              // should be equal to @p{\sum_j w_{ij} p[j]} with @p{p[j]} be the
              // values of the degrees of freedom on the coarse grid. we can
              // thus compute constraints which link the degrees of freedom
              // @p{v[i]} on the fine grid to those on the coarse grid,
              // @p{p[j]}. Now to use these as real constraints, rather than as
              // additional equations, we have to identify representants among
              // the @p{i} for each @p{j}. this will be done by simply taking
              // the first @p{i} for which @p{w_{ij}==1}.
              //
              // guard modification of the weights array by a Mutex. since it
              // should happen rather rarely that there are several threads
              // operating on different intergrid weights, have only one mutex
              // for all of them
              for (types::global_dof_index i = 0;
                   i < copy_data.global_parameter_representation[pos].size();
                   ++i)
                // set this weight if it belongs to a parameter dof.
                if (weight_mapping[i] != numbers::invalid_dof_index)
                  {
                    // only overwrite old value if not by zero
                    if (copy_data.global_parameter_representation[pos](i) != 0)
                      {
                        const types::global_dof_index
                          wi = copy_data.parameter_dof_indices[local_dof],
                          wj = weight_mapping[i];
                        weights[wi][wj] =
                          copy_data.global_parameter_representation[pos](i);
                      }
                  }
                else if (!is_called_in_parallel)
                  {
                    // Note that when this function operates with distributed
                    // fine grid, this assertion is switched off since the
                    // condition does not necessarily hold
                    Assert(copy_data.global_parameter_representation[pos](i) ==
                             0,
                           ExcInternalError());
                  }

              ++pos;
            }
      }



      /**
       * This is a helper function that is used in the computation of intergrid
       * constraints. See the function for a thorough description of how it
       * works.
       */
      template <int dim, int spacedim>
      void
      compute_intergrid_weights_2(
        const dealii::DoFHandler<dim, spacedim> &coarse_grid,
        const unsigned int                       coarse_component,
        const InterGridMap<dealii::DoFHandler<dim, spacedim>>
          &                                         coarse_to_fine_grid_map,
        const std::vector<dealii::Vector<double>> & parameter_dofs,
        const std::vector<types::global_dof_index> &weight_mapping,
        std::vector<std::map<types::global_dof_index, float>> &weights)
      {
        Assembler::Scratch                 scratch;
        Assembler::CopyData<dim, spacedim> copy_data;

        unsigned int n_interesting_dofs = 0;
        for (unsigned int local_dof = 0;
             local_dof < coarse_grid.get_fe().dofs_per_cell;
             ++local_dof)
          if (coarse_grid.get_fe().system_to_component_index(local_dof).first ==
              coarse_component)
            ++n_interesting_dofs;

        copy_data.global_parameter_representation.resize(n_interesting_dofs);

        bool is_called_in_parallel = false;
        for (std::size_t i = 0;
             i < copy_data.global_parameter_representation.size();
             ++i)
          {
#ifdef DEAL_II_WITH_MPI
            MPI_Comm communicator = MPI_COMM_SELF;
            try
              {
                const typename dealii::parallel::TriangulationBase<dim,
                                                                   spacedim>
                  &tria = dynamic_cast<const typename dealii::parallel::
                                         TriangulationBase<dim, spacedim> &>(
                    coarse_to_fine_grid_map.get_destination_grid()
                      .get_triangulation());
                communicator          = tria.get_communicator();
                is_called_in_parallel = true;
              }
            catch (std::bad_cast &)
              {
                // Nothing bad happened: the user used serial Triangulation
              }


            IndexSet locally_relevant_dofs;
            DoFTools::extract_locally_relevant_dofs(
              coarse_to_fine_grid_map.get_destination_grid(),
              locally_relevant_dofs);

            copy_data.global_parameter_representation[i].reinit(
              coarse_to_fine_grid_map.get_destination_grid()
                .locally_owned_dofs(),
              locally_relevant_dofs,
              communicator);
#else
            const types::global_dof_index n_fine_dofs = weight_mapping.size();
            copy_data.global_parameter_representation[i].reinit(n_fine_dofs);
#endif
          }

        auto worker =
          [coarse_component,
           &coarse_grid,
           &coarse_to_fine_grid_map,
           &parameter_dofs](const typename dealii::DoFHandler<dim, spacedim>::
                              active_cell_iterator &            cell,
                            const Assembler::Scratch &          scratch_data,
                            Assembler::CopyData<dim, spacedim> &copy_data) {
            compute_intergrid_weights_3<dim, spacedim>(cell,
                                                       scratch_data,
                                                       copy_data,
                                                       coarse_component,
                                                       coarse_grid.get_fe(),
                                                       coarse_to_fine_grid_map,
                                                       parameter_dofs);
          };

        auto copier =
          [coarse_component,
           &coarse_grid,
           &weight_mapping,
           is_called_in_parallel,
           &weights](const Assembler::CopyData<dim, spacedim> &copy_data) {
            copy_intergrid_weights_3<dim, spacedim>(copy_data,
                                                    coarse_component,
                                                    coarse_grid.get_fe(),
                                                    weight_mapping,
                                                    is_called_in_parallel,
                                                    weights);
          };

        WorkStream::run(coarse_grid.begin_active(),
                        coarse_grid.end(),
                        worker,
                        copier,
                        scratch,
                        copy_data);

#ifdef DEAL_II_WITH_MPI
        for (std::size_t i = 0;
             i < copy_data.global_parameter_representation.size();
             ++i)
          copy_data.global_parameter_representation[i].update_ghost_values();
#endif
      }



      /**
       * This is a helper function that is used in the computation of integrid
       * constraints. See the function for a thorough description of how it
       * works.
       */
      template <int dim, int spacedim>
      unsigned int
      compute_intergrid_weights_1(
        const dealii::DoFHandler<dim, spacedim> &coarse_grid,
        const unsigned int                       coarse_component,
        const dealii::DoFHandler<dim, spacedim> &fine_grid,
        const unsigned int                       fine_component,
        const InterGridMap<dealii::DoFHandler<dim, spacedim>>
          &coarse_to_fine_grid_map,
        std::vector<std::map<types::global_dof_index, float>> &weights,
        std::vector<types::global_dof_index> &                 weight_mapping)
      {
        // aliases to the finite elements used by the dof handlers:
        const FiniteElement<dim, spacedim> &coarse_fe = coarse_grid.get_fe(),
                                           &fine_fe   = fine_grid.get_fe();

        // global numbers of dofs
        const types::global_dof_index n_coarse_dofs = coarse_grid.n_dofs(),
                                      n_fine_dofs   = fine_grid.n_dofs();

        // local numbers of dofs
        const unsigned int fine_dofs_per_cell = fine_fe.dofs_per_cell;

        // alias the number of dofs per cell belonging to the coarse_component
        // which is to be the restriction of the fine grid:
        const unsigned int coarse_dofs_per_cell_component =
          coarse_fe
            .base_element(
              coarse_fe.component_to_base_index(coarse_component).first)
            .dofs_per_cell;


        // Try to find out whether the grids stem from the same coarse grid.
        // This is a rather crude test, but better than nothing
        Assert(coarse_grid.get_triangulation().n_cells(0) ==
                 fine_grid.get_triangulation().n_cells(0),
               ExcGridsDontMatch());

        // check whether the map correlates the right objects
        Assert(&coarse_to_fine_grid_map.get_source_grid() == &coarse_grid,
               ExcGridsDontMatch());
        Assert(&coarse_to_fine_grid_map.get_destination_grid() == &fine_grid,
               ExcGridsDontMatch());


        // check whether component numbers are valid
        AssertIndexRange(coarse_component, coarse_fe.n_components());
        AssertIndexRange(fine_component, fine_fe.n_components());

        // check whether respective finite elements are equal
        Assert(coarse_fe.base_element(
                 coarse_fe.component_to_base_index(coarse_component).first) ==
                 fine_fe.base_element(
                   fine_fe.component_to_base_index(fine_component).first),
               ExcFiniteElementsDontMatch());

#ifdef DEBUG
        // if in debug mode, check whether the coarse grid is indeed coarser
        // everywhere than the fine grid
        for (typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator
               cell = coarse_grid.begin_active();
             cell != coarse_grid.end();
             ++cell)
          Assert(cell->level() <= coarse_to_fine_grid_map[cell]->level(),
                 ExcGridNotCoarser());
#endif

        /*
         * From here on: the term `parameter' refers to the selected component
         * on the coarse grid and its analogon on the fine grid. The naming of
         * variables containing this term is due to the fact that
         * `selected_component' is longer, but also due to the fact that the
         * code of this function was initially written for a program where the
         * component which we wanted to match between grids was actually the
         * `parameter' variable.
         *
         * Likewise, the terms `parameter grid' and `state grid' refer to the
         * coarse and fine grids, respectively.
         *
         * Changing the names of variables would in principle be a good idea,
         * but would not make things simpler and would be another source of
         * errors. If anyone feels like doing so: patches would be welcome!
         */



        // set up vectors of cell-local data; each vector represents one degree
        // of freedom of the coarse-grid variable in the fine-grid element
        std::vector<dealii::Vector<double>> parameter_dofs(
          coarse_dofs_per_cell_component,
          dealii::Vector<double>(fine_dofs_per_cell));
        // for each coarse dof: find its position within the fine element and
        // set this value to one in the respective vector (all other values are
        // zero by construction)
        for (unsigned int local_coarse_dof = 0;
             local_coarse_dof < coarse_dofs_per_cell_component;
             ++local_coarse_dof)
          for (unsigned int fine_dof = 0; fine_dof < fine_fe.dofs_per_cell;
               ++fine_dof)
            if (fine_fe.system_to_component_index(fine_dof) ==
                std::make_pair(fine_component, local_coarse_dof))
              {
                parameter_dofs[local_coarse_dof](fine_dof) = 1.;
                break;
              }


        // find out how many DoFs there are on the grids belonging to the
        // components we want to match
        unsigned int n_parameters_on_fine_grid = 0;
        {
          // have a flag for each dof on the fine grid and set it to true if
          // this is an interesting dof. finally count how many true's there
          std::vector<bool> dof_is_interesting(fine_grid.n_dofs(), false);
          std::vector<types::global_dof_index> local_dof_indices(
            fine_fe.dofs_per_cell);

          for (typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator
                 cell = fine_grid.begin_active();
               cell != fine_grid.end();
               ++cell)
            if (cell->is_locally_owned())
              {
                cell->get_dof_indices(local_dof_indices);
                for (unsigned int i = 0; i < fine_fe.dofs_per_cell; ++i)
                  if (fine_fe.system_to_component_index(i).first ==
                      fine_component)
                    dof_is_interesting[local_dof_indices[i]] = true;
              }

          n_parameters_on_fine_grid = std::count(dof_is_interesting.begin(),
                                                 dof_is_interesting.end(),
                                                 true);
        }


        // set up the weights mapping
        weights.clear();
        weights.resize(n_coarse_dofs);

        weight_mapping.clear();
        weight_mapping.resize(n_fine_dofs, numbers::invalid_dof_index);

        {
          std::vector<types::global_dof_index> local_dof_indices(
            fine_fe.dofs_per_cell);
          unsigned int next_free_index = 0;
          for (typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator
                 cell = fine_grid.begin_active();
               cell != fine_grid.end();
               ++cell)
            if (cell->is_locally_owned())
              {
                cell->get_dof_indices(local_dof_indices);
                for (unsigned int i = 0; i < fine_fe.dofs_per_cell; ++i)
                  // if this DoF is a parameter dof and has not yet been
                  // numbered, then do so
                  if ((fine_fe.system_to_component_index(i).first ==
                       fine_component) &&
                      (weight_mapping[local_dof_indices[i]] ==
                       numbers::invalid_dof_index))
                    {
                      weight_mapping[local_dof_indices[i]] = next_free_index;
                      ++next_free_index;
                    }
              }

          Assert(next_free_index == n_parameters_on_fine_grid,
                 ExcInternalError());
        }


        // for each cell on the parameter grid: find out which degrees of
        // freedom on the fine grid correspond in which way to the degrees of
        // freedom on the parameter grid
        //
        // do this in a separate function to allow for multithreading there. see
        // this function also if you want to read more information on the
        // algorithm used.
        compute_intergrid_weights_2(coarse_grid,
                                    coarse_component,
                                    coarse_to_fine_grid_map,
                                    parameter_dofs,
                                    weight_mapping,
                                    weights);


        // ok, now we have all weights for each dof on the fine grid. if in
        // debug mode lets see if everything went smooth, i.e. each dof has sum
        // of weights one
        //
        // in other words this means that if the sum of all shape functions on
        // the parameter grid is one (which is always the case), then the
        // representation on the state grid should be as well (division of
        // unity)
        //
        // if the parameter grid has more than one component, then the
        // respective dofs of the other components have sum of weights zero, of
        // course. we do not explicitly ask which component a dof belongs to,
        // but this at least tests some errors
#ifdef DEBUG
        for (unsigned int col = 0; col < n_parameters_on_fine_grid; ++col)
          {
            double sum = 0;
            for (types::global_dof_index row = 0; row < n_coarse_dofs; ++row)
              if (weights[row].find(col) != weights[row].end())
                sum += weights[row][col];
            Assert((std::fabs(sum - 1) < 1.e-12) ||
                     ((coarse_fe.n_components() > 1) && (sum == 0)),
                   ExcInternalError());
          }
#endif


        return n_parameters_on_fine_grid;
      }


    } // namespace
  }   // namespace internal



  template <int dim, int spacedim>
  void
  compute_intergrid_constraints(
    const DoFHandler<dim, spacedim> &              coarse_grid,
    const unsigned int                             coarse_component,
    const DoFHandler<dim, spacedim> &              fine_grid,
    const unsigned int                             fine_component,
    const InterGridMap<DoFHandler<dim, spacedim>> &coarse_to_fine_grid_map,
    AffineConstraints<double> &                    constraints)
  {
    // store the weights with which a dof on the parameter grid contributes to a
    // dof on the fine grid. see the long doc below for more info
    //
    // allocate as many rows as there are parameter dofs on the coarse grid and
    // as many columns as there are parameter dofs on the fine grid.
    //
    // weight_mapping is used to map the global (fine grid) parameter dof
    // indices to the columns
    //
    // in the original implementation, the weights array was actually of
    // FullMatrix<double> type. this wasted huge amounts of memory, but was
    // fast. nonetheless, since the memory consumption was quadratic in the
    // number of degrees of freedom, this was not very practical, so we now use
    // a vector of rows of the matrix, and in each row a vector of pairs
    // (colnum,value). this seems like the best tradeoff between memory and
    // speed, as it is now linear in memory and still fast enough.
    //
    // to save some memory and since the weights are usually (negative) powers
    // of 2, we choose the value type of the matrix to be @p{float} rather than
    // @p{double}.
    std::vector<std::map<types::global_dof_index, float>> weights;

    // this is this mapping. there is one entry for each dof on the fine grid;
    // if it is a parameter dof, then its value is the column in weights for
    // that parameter dof, if it is any other dof, then its value is -1,
    // indicating an error
    std::vector<types::global_dof_index> weight_mapping;

    const unsigned int n_parameters_on_fine_grid =
      internal::compute_intergrid_weights_1(coarse_grid,
                                            coarse_component,
                                            fine_grid,
                                            fine_component,
                                            coarse_to_fine_grid_map,
                                            weights,
                                            weight_mapping);
    (void)n_parameters_on_fine_grid;

    // global numbers of dofs
    const types::global_dof_index n_coarse_dofs = coarse_grid.n_dofs(),
                                  n_fine_dofs   = fine_grid.n_dofs();


    // get an array in which we store which dof on the coarse grid is a
    // parameter and which is not
    std::vector<bool> coarse_dof_is_parameter(coarse_grid.n_dofs());
    if (true)
      {
        std::vector<bool> mask(coarse_grid.get_fe(0).n_components(), false);
        mask[coarse_component] = true;
        extract_dofs<DoFHandler<dim, spacedim>>(coarse_grid,
                                                ComponentMask(mask),
                                                coarse_dof_is_parameter);
      }

    // now we know that the weights in each row constitute a constraint. enter
    // this into the constraints object
    //
    // first task: for each parameter dof on the parameter grid, find a
    // representant on the fine, global grid. this is possible since we use
    // conforming finite element. we take this representant to be the first
    // element in this row with weight identical to one. the representant will
    // become an unconstrained degree of freedom, while all others will be
    // constrained to this dof (and possibly others)
    std::vector<types::global_dof_index> representants(
      n_coarse_dofs, numbers::invalid_dof_index);
    for (types::global_dof_index parameter_dof = 0;
         parameter_dof < n_coarse_dofs;
         ++parameter_dof)
      if (coarse_dof_is_parameter[parameter_dof] == true)
        {
          // if this is the line of a parameter dof on the coarse grid, then it
          // should have at least one dependent node on the fine grid
          Assert(weights[parameter_dof].size() > 0, ExcInternalError());

          // find the column where the representant is mentioned
          std::map<types::global_dof_index, float>::const_iterator i =
            weights[parameter_dof].begin();
          for (; i != weights[parameter_dof].end(); ++i)
            if (i->second == 1)
              break;
          Assert(i != weights[parameter_dof].end(), ExcInternalError());
          const types::global_dof_index column = i->first;

          // now we know in which column of weights the representant is, but we
          // don't know its global index. get it using the inverse operation of
          // the weight_mapping
          types::global_dof_index global_dof = 0;
          for (; global_dof < weight_mapping.size(); ++global_dof)
            if (weight_mapping[global_dof] ==
                static_cast<types::global_dof_index>(column))
              break;
          Assert(global_dof < weight_mapping.size(), ExcInternalError());

          // now enter the representants global index into our list
          representants[parameter_dof] = global_dof;
        }
      else
        {
          // consistency check: if this is no parameter dof on the coarse grid,
          // then the respective row must be empty!
          Assert(weights[parameter_dof].size() == 0, ExcInternalError());
        }



    // note for people that want to optimize this function: the largest part of
    // the computing time is spent in the following, rather innocent block of
    // code. basically, it must be the AffineConstraints::add_entry call which
    // takes the bulk of the time, but it is not known to the author how to make
    // it faster...
    std::vector<std::pair<types::global_dof_index, double>> constraint_line;
    for (types::global_dof_index global_dof = 0; global_dof < n_fine_dofs;
         ++global_dof)
      if (weight_mapping[global_dof] != numbers::invalid_dof_index)
        // this global dof is a parameter dof, so it may carry a constraint note
        // that for each global dof, the sum of weights shall be one, so we can
        // find out whether this dof is constrained in the following way: if the
        // only weight in this row is a one, and the representant for the
        // parameter dof of the line in which this one is is the present dof,
        // then we consider this dof to be unconstrained. otherwise, all other
        // dofs are constrained
        {
          const types::global_dof_index col = weight_mapping[global_dof];
          Assert(col < n_parameters_on_fine_grid, ExcInternalError());

          types::global_dof_index first_used_row = 0;

          {
            Assert(weights.size() > 0, ExcInternalError());
            std::map<types::global_dof_index, float>::const_iterator col_entry =
              weights[0].end();
            for (; first_used_row < n_coarse_dofs; ++first_used_row)
              {
                col_entry = weights[first_used_row].find(col);
                if (col_entry != weights[first_used_row].end())
                  break;
              }

            Assert(col_entry != weights[first_used_row].end(),
                   ExcInternalError());

            if ((col_entry->second == 1) &&
                (representants[first_used_row] == global_dof))
              // dof unconstrained or constrained to itself (in case this cell
              // is mapped to itself, rather than to children of itself)
              continue;
          }


          // otherwise enter all constraints
          constraints.add_line(global_dof);

          constraint_line.clear();
          for (types::global_dof_index row = first_used_row;
               row < n_coarse_dofs;
               ++row)
            {
              const std::map<types::global_dof_index, float>::const_iterator j =
                weights[row].find(col);
              if ((j != weights[row].end()) && (j->second != 0))
                constraint_line.emplace_back(representants[row], j->second);
            }

          constraints.add_entries(global_dof, constraint_line);
        }
  }



  template <int dim, int spacedim>
  void
  compute_intergrid_transfer_representation(
    const DoFHandler<dim, spacedim> &              coarse_grid,
    const unsigned int                             coarse_component,
    const DoFHandler<dim, spacedim> &              fine_grid,
    const unsigned int                             fine_component,
    const InterGridMap<DoFHandler<dim, spacedim>> &coarse_to_fine_grid_map,
    std::vector<std::map<types::global_dof_index, float>>
      &transfer_representation)
  {
    // store the weights with which a dof on the parameter grid contributes to a
    // dof on the fine grid. see the long doc below for more info
    //
    // allocate as many rows as there are parameter dofs on the coarse grid and
    // as many columns as there are parameter dofs on the fine grid.
    //
    // weight_mapping is used to map the global (fine grid) parameter dof
    // indices to the columns
    //
    // in the original implementation, the weights array was actually of
    // FullMatrix<double> type. this wasted huge amounts of memory, but was
    // fast. nonetheless, since the memory consumption was quadratic in the
    // number of degrees of freedom, this was not very practical, so we now use
    // a vector of rows of the matrix, and in each row a vector of pairs
    // (colnum,value). this seems like the best tradeoff between memory and
    // speed, as it is now linear in memory and still fast enough.
    //
    // to save some memory and since the weights are usually (negative) powers
    // of 2, we choose the value type of the matrix to be @p{float} rather than
    // @p{double}.
    std::vector<std::map<types::global_dof_index, float>> weights;

    // this is this mapping. there is one entry for each dof on the fine grid;
    // if it is a parameter dof, then its value is the column in weights for
    // that parameter dof, if it is any other dof, then its value is -1,
    // indicating an error
    std::vector<types::global_dof_index> weight_mapping;

    internal::compute_intergrid_weights_1(coarse_grid,
                                          coarse_component,
                                          fine_grid,
                                          fine_component,
                                          coarse_to_fine_grid_map,
                                          weights,
                                          weight_mapping);

    // now compute the requested representation
    const types::global_dof_index n_global_parm_dofs =
      std::count_if(weight_mapping.begin(),
                    weight_mapping.end(),
                    [](const types::global_dof_index dof) {
                      return dof != numbers::invalid_dof_index;
                    });

    // first construct the inverse mapping of weight_mapping
    std::vector<types::global_dof_index> inverse_weight_mapping(
      n_global_parm_dofs, numbers::invalid_dof_index);
    for (types::global_dof_index i = 0; i < weight_mapping.size(); ++i)
      {
        const types::global_dof_index parameter_dof = weight_mapping[i];
        // if this global dof is a parameter
        if (parameter_dof != numbers::invalid_dof_index)
          {
            Assert(parameter_dof < n_global_parm_dofs, ExcInternalError());
            Assert((inverse_weight_mapping[parameter_dof] ==
                    numbers::invalid_dof_index),
                   ExcInternalError());

            inverse_weight_mapping[parameter_dof] = i;
          }
      }

    // next copy over weights array and replace respective numbers
    const types::global_dof_index n_rows = weight_mapping.size();

    transfer_representation.clear();
    transfer_representation.resize(n_rows);

    const types::global_dof_index n_coarse_dofs = coarse_grid.n_dofs();
    for (types::global_dof_index i = 0; i < n_coarse_dofs; ++i)
      {
        std::map<types::global_dof_index, float>::const_iterator j =
          weights[i].begin();
        for (; j != weights[i].end(); ++j)
          {
            const types::global_dof_index p = inverse_weight_mapping[j->first];
            Assert(p < n_rows, ExcInternalError());

            transfer_representation[p][i] = j->second;
          }
      }
  }



  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  make_zero_boundary_constraints(
    const DoFHandlerType<dim, spacedim> &dof,
    const types::boundary_id             boundary_id,
    AffineConstraints<number> &          zero_boundary_constraints,
    const ComponentMask &                component_mask)
  {
    Assert(component_mask.represents_n_components(dof.get_fe(0).n_components()),
           ExcMessage("The number of components in the mask has to be either "
                      "zero or equal to the number of components in the finite "
                      "element."));

    const unsigned int n_components = DoFTools::n_components(dof);

    Assert(component_mask.n_selected_components(n_components) > 0,
           ComponentMask::ExcNoComponentSelected());

    // a field to store the indices on the face
    std::vector<types::global_dof_index> face_dofs;
    face_dofs.reserve(max_dofs_per_face(dof));
    // a field to store the indices on the cell
    std::vector<types::global_dof_index> cell_dofs;
    cell_dofs.reserve(max_dofs_per_cell(dof));

    typename DoFHandlerType<dim, spacedim>::active_cell_iterator
      cell = dof.begin_active(),
      endc = dof.end();
    for (; cell != endc; ++cell)
      if (!cell->is_artificial() && cell->at_boundary())
        {
          const FiniteElement<dim, spacedim> &fe = cell->get_fe();

          // get global indices of dofs on the cell
          cell_dofs.resize(fe.dofs_per_cell);
          cell->get_dof_indices(cell_dofs);

          for (unsigned int face_no = 0;
               face_no < GeometryInfo<dim>::faces_per_cell;
               ++face_no)
            {
              const typename DoFHandlerType<dim, spacedim>::face_iterator face =
                cell->face(face_no);

              // if face is on the boundary and satisfies the correct boundary
              // id property
              if (face->at_boundary() &&
                  ((boundary_id == numbers::invalid_boundary_id) ||
                   (face->boundary_id() == boundary_id)))
                {
                  // get indices and physical location on this face
                  face_dofs.resize(fe.dofs_per_face);
                  face->get_dof_indices(face_dofs, cell->active_fe_index());

                  // enter those dofs into the list that match the component
                  // signature.
                  for (const types::global_dof_index face_dof : face_dofs)
                    {
                      // Find out if a dof has a contribution in this component,
                      // and if so, add it to the list
                      const std::vector<types::global_dof_index>::iterator
                        it_index_on_cell = std::find(cell_dofs.begin(),
                                                     cell_dofs.end(),
                                                     face_dof);
                      Assert(it_index_on_cell != cell_dofs.end(),
                             ExcInvalidIterator());
                      const unsigned int index_on_cell =
                        std::distance(cell_dofs.begin(), it_index_on_cell);
                      const ComponentMask &nonzero_component_array =
                        cell->get_fe().get_nonzero_components(index_on_cell);
                      bool nonzero = false;
                      for (unsigned int c = 0; c < n_components; ++c)
                        if (nonzero_component_array[c] && component_mask[c])
                          {
                            nonzero = true;
                            break;
                          }

                      if (nonzero)
                        zero_boundary_constraints.add_line(face_dof);
                    }
                }
            }
        }
  }



  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  make_zero_boundary_constraints(
    const DoFHandlerType<dim, spacedim> &dof,
    AffineConstraints<number> &          zero_boundary_constraints,
    const ComponentMask &                component_mask)
  {
    make_zero_boundary_constraints(dof,
                                   numbers::invalid_boundary_id,
                                   zero_boundary_constraints,
                                   component_mask);
  }


} // end of namespace DoFTools



// explicit instantiations

#include "dof_tools_constraints.inst"



DEAL_II_NAMESPACE_CLOSE
