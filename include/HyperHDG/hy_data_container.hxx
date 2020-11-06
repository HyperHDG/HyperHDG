#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>

/*!*************************************************************************************************
 * \brief   Class for saving some (abstract) data per hyperedge.
 *
 * For parabolic PDEs, for example, it might be necessary to save some data of the old time step.
 * This data might be related to hypernodes and therefore be saved in a global vector that is
 * administrated by the HypernodeFactory, or it might be related to hypernodes and therefore be
 * administrated by the HyDataContainer.
 *
 * \tparam  data_t          The class name of the data which is saved per hyperedge.
 * \tparam  vectorT         Data structure in which the datas are stored. Defaults to std::vector.
 * \tparam  hyEdge_index_t  Index type which is used to identify entries in this data structure.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <typename data_t,
          typename vectorT = std::vector<data_t>,
          typename hyEdge_index_t = decltype(std::declval<vectorT>().size())>
class HyDataContainer
{
 private:
  /*!***********************************************************************************************
   * \brief   The internal overall data container which holds all data of a hypergraph's hyperedges.
   ************************************************************************************************/
  vectorT data_container;

 public:
  /*!***********************************************************************************************
   * \brief   Defines the return value of the class.
   ************************************************************************************************/
  typedef data_t value_type;
  /*!***********************************************************************************************
   * \brief   Defines the value type of input argument for standard constructor.
   ************************************************************************************************/
  typedef decltype(std::declval<vectorT>().size()) constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Construct a topology from a given filename.
   *
   * \param   n_hyEdges     Number of hyperedges.
   ************************************************************************************************/
  HyDataContainer(const constructor_value_type n_hyEdges) : data_container(n_hyEdges) {}

  /*!***********************************************************************************************
   * \brief   Get data of hyperedge of given index.
   *
   * \param   index         The index of the hyperedge whose data is to be returned.
   * \retval  data          The data of the selected hyperedge.
   ************************************************************************************************/
  value_type& operator[](const hyEdge_index_t index) { return get_hyEdge(index); }
  /*!***********************************************************************************************
   * \brief   Get data of hyperedge of given index.
   *
   * \param   index         The index of the hyperedge whose data is to be returned.
   * \retval  data          The data of the selected hyperedge.
   ************************************************************************************************/
  value_type& get_hyEdge(const hyEdge_index_t index)
  {
    hy_assert(index < data_container.size() && index >= 0,
              "Index must be non-negative and smaller than "
                << data_container.size() << " (which is the amount of hyperedges). It was " << index
                << "!");
    return data_container[index];
  }
};  // end of class File
