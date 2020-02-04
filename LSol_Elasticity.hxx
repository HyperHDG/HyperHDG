/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 * 
 * Definition local solver class: In this case for diffusion.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef ELASTICITYSOLVER_HXX
#define ELASTICITYSOLVER_HXX

#include <FuncAndQuad.hxx>
#include <Hypercube.hxx>

#include <array>


template<unsigned int poly_deg, unsigned int quad_deg, class GeomT>
class ElasticRods
{
  public:
    static constexpr unsigned int hyEdge_dimension() { return GeomT::hyEdge_dim(); }
    static constexpr bool use_geometry() { return true; }
    static constexpr unsigned int n_glob_dofs_per_node()
    {
      unsigned int amount = 1;
      for (unsigned int iteration = 0; iteration < hyEdge_dim - 1; ++ iteration)  amount *= poly_deg + 1;
      return GeomT::space_dim() * amount;
    }
  private:
    static constexpr unsigned int hyEdge_dim   = GeomT::space_dim(),
                                  space_dim    = GeomT::space_dim(),
                                  n_quads_     = FuncQuad::compute_n_quad_points(quad_deg),
                                  n_shape_fct_ = n_glob_dofs_per_node() * (poly_deg + 1),
                                  n_shape_bdr_ = n_glob_dofs_per_node();
    const double tau_;
    const std::array<double, n_quads_> q_weights_;
    const std::array< std::array<double, n_quads_ > , poly_deg + 1 > trial_;
    const std::array< std::array<double, 2 > , poly_deg + 1 > trial_bdr_;
    const std::array<double, (hyEdge_dim+1) * n_shape_fct_ * (hyEdge_dim+1) * n_shape_fct_> loc_mat_;
    
    inline std::array<double, (hyEdge_dim+1) * n_shape_fct_> assemble_rhs(const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    inline std::array<double, (hyEdge_dim+1) * n_shape_fct_> solve_local_problem(const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    
    inline std::array< std::array<double, n_shape_bdr_> , 2 * hyEdge_dim > dual_at_boundary(const std::array<double, (hyEdge_dim+1) * n_shape_fct_>& coeffs) const;
    inline std::array< std::array<double, n_shape_bdr_> , 2 * hyEdge_dim > primal_at_boundary(const std::array<double, (hyEdge_dim+1) * n_shape_fct_>& coeffs) const;
    
    inline std::array< std::array<double, n_shape_bdr_>, 2*hyEdge_dim > node_dof_to_edge_dof(const std::array< std::array<double, GeomT::space_dim() * n_shape_bdr_>, 2*hyEdge_dim > lambda, const GeomT& geom) const;
    inline std::array< std::array<double, GeomT::space_dim() * n_shape_bdr_>, 2*hyEdge_dim > edge_dof_to_node_dof(const std::array< std::array<double, n_shape_bdr_>, 2*hyEdge_dim > lambda, const GeomT& geom) const;
  public:
    typedef double constructor_value_type;
    ElasticRods(const constructor_value_type& tau);
    
    // Function for matrix--vector multiply    
    std::array< std::array<double, GeomT::space_dim() * n_shape_bdr_> , 2 * hyEdge_dim >
      numerical_flux_from_lambda(const std::array< std::array<double, GeomT::space_dim() * n_shape_bdr_> , 2*hyEdge_dim >& lambda_values, const GeomT& geom) const;
    
    // Plotting functions
    template<unsigned int sizeT> std::array<double, Hypercube<hyEdge_dim>::pow(sizeT)> primal_at_dyadic
      (const std::array<double, sizeT>& abscissas, const std::array< std::array<double, GeomT::space_dim() * n_shape_bdr_> , 2*hyEdge_dim >& lambda_values, const GeomT& geom) const;
    template<unsigned int sizeT> std::array< std::array<double,hyEdge_dim> , Hypercube<hyEdge_dim>::pow(sizeT) > dual_at_dyadic
      (const std::array<double, sizeT>& abscissas, const std::array< std::array<double, GeomT::space_dim() * n_shape_bdr_> , 2*hyEdge_dim >& lambda_values, const GeomT& geom) const;
};

#endif
