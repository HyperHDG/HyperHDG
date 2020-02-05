#ifndef HYPERCUBE_HXX
#define HYPERCUBE_HXX

/**
 * \brief Helper class containing numbers and functions related to hypercubes.
 *
 * \todo All functions related to dim-dimensional tensor can go here.
 */
template <int dim>
struct Hypercube
{
  /// Number of faces of the hypercube
  static constexpr unsigned int n_faces ()
  {
    return 2*dim;
  }

  /// Number of vertices of the hypercube
  static constexpr unsigned int n_vertices ()
  {
    return 1 << dim;
  }
  
  /// Return `n` to the power `dim`, which is the size of the `dim`-dimensional tensor of dimension `n`.
  static constexpr unsigned int pow(unsigned int n)
  {
    int result = 1;
    for (unsigned int i=0;i<dim;++i)
      result *= n;
    return result;
  }

};

#endif // end of ifndef HYPERCUBE_HXX
