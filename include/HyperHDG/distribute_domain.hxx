#pragma once  // Ensure that file is included only once in a single compilation.

#ifdef HYPERHDG_PETSC

#include <HyperHDG/read_domain.hxx>

#include <petsc.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

/*!*************************************************************************************************
 * \brief   Assign each (global) hypernode to an owning rank.
 *
 * Placeholder partitioner: round-robin over ranks. It maximally scatters the numbering and so
 * exercises the ghost/scatter/local-to-global machinery as hard as possible, while requiring no
 * external library. This is intended to be replaced by a graph partitioner (sequential KaHIP
 * \c kaffpa, vertex-weighted by degree) that minimises the edge cut.
 *
 * \param   v       Global hypernode index.
 * \param   size    Number of ranks (= number of parts).
 * \retval  rank    Owning rank in [0, size).
 *
 * \authors   Joseph Holten, Karlsruhe Institute of Technology, 2026.
 **************************************************************************************************/
template <typename hyNode_index_t>
int compute_partition_placeholder(const hyNode_index_t v, const int size)
{
  return static_cast<int>(v % static_cast<hyNode_index_t>(size));
}

/*!*************************************************************************************************
 * \brief   Read, partition, renumber and distribute a domain across ranks (hyEdge_dim == 1).
 *
 * First draft: rank 0 reads the whole domain (sequential I/O), computes a partition and an
 * owned-contiguous global renumbering, then scatters to each rank its locally owned hyperedges
 * together with the points and per-edge data they reference. Each rank receives a \b local
 * \c DomainInfo whose hypernodes are numbered owned-first (ghosts last) and which carries the
 * local-to-global map (\c lgmap) plus \c n_owned_hyNodes / \c n_global_hyNodes.
 *
 * For \c hyEdge_dim == 1 a hypernode is a graph vertex and coincides with a point, so partitioning
 * the points partitions the hypernodes. Each hyperedge is owned by exactly one rank
 * (\c owner = min part of its two endpoints); a referenced endpoint owned by another rank becomes a
 * ghost on the owning rank.
 *
 * \authors   Joseph Holten, Karlsruhe Institute of Technology, 2026.
 **************************************************************************************************/
template <unsigned int hyEdge_dim,
          unsigned int space_dim,
          template <typename...> typename vectorT = std::vector,
          typename pointT = Point<space_dim, double>,
          typename hyEdge_index_t = unsigned int,
          typename hyNode_index_t = hyEdge_index_t,
          typename pt_index_t = hyNode_index_t>
DomainInfo<hyEdge_dim, space_dim, vectorT, pointT, hyEdge_index_t, hyNode_index_t, pt_index_t>
distribute_domain(const std::string& filename, MPI_Comm comm)
{
  static_assert(hyEdge_dim == 1,
                "distribute_domain only supports hyEdge_dim == 1 (networks) so far.");
  static_assert(sizeof(hyEdge_index_t) == sizeof(PetscInt),
                "HyperHDG index type must have the same size as PetscInt, as matrix COO indices "
                "are passed to PETSc by reinterpretation.");

  using DI = DomainInfo<hyEdge_dim, space_dim, vectorT, pointT, hyEdge_index_t, hyNode_index_t,
                        pt_index_t>;
  using node_pair_t = std::array<hyNode_index_t, 2 * hyEdge_dim>;  // == array<.,2>
  using point_pair_t = std::array<pt_index_t, 1 << hyEdge_dim>;    // == array<.,2>
  using value_t = typename pointT::value_type;

  enum Tag
  {
    TAG_HDR = 100,
    TAG_LGMAP,
    TAG_PTS,
    TAG_EDGES,
    TAG_TYPES,
    TAG_PROPS
  };

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // The local arrays defining one rank's chunk.
  struct Chunk
  {
    hyNode_index_t n_owned = 0, n_local = 0, n_global = 0;
    hyEdge_index_t n_oe = 0;
    unsigned int nprop = 0;
    std::vector<hyNode_index_t> lgmap;  // local -> global hypernode
    std::vector<pointT> pts;            // owned + ghost point coordinates (local order)
    std::vector<node_pair_t> edges;     // owned edges in local node numbering
    std::vector<node_pair_t> types;     // hyFaces_hyEdge for owned edges
    std::vector<value_t> props;         // n_oe * nprop, flattened
  };

  // Build the local DomainInfo from a (received or locally produced) chunk.
  auto build_domain_info = [](const Chunk& c) -> DI
  {
    DI di(c.n_local, c.n_oe, c.n_local, c.n_local);
    di.n_owned_hyNodes = c.n_owned;
    di.n_global_hyNodes = c.n_global;
    di.n_properties = c.nprop;
    di.lgmap = vectorT<hyNode_index_t>(c.lgmap.begin(), c.lgmap.end());
    for (hyNode_index_t l = 0; l < c.n_local; ++l)
      di.points[l] = c.pts[l];
    for (hyEdge_index_t i = 0; i < c.n_oe; ++i)
    {
      di.hyNodes_hyEdge[i] = c.edges[i];
      di.points_hyEdge[i] = point_pair_t{c.edges[i][0], c.edges[i][1]};
      di.hyFaces_hyEdge[i] = c.types[i];
    }
    if (c.nprop)
    {
      di.hyEdge_properties.resize(c.n_oe);
      for (hyEdge_index_t i = 0; i < c.n_oe; ++i)
        di.hyEdge_properties[i].assign(c.props.begin() + i * c.nprop,
                                       c.props.begin() + (i + 1) * c.nprop);
    }
    return di;
  };

  if (rank != 0)
  {
    // Receive this rank's chunk from rank 0.
    Chunk c;
    long long hdr[5];
    MPI_Recv(hdr, 5, MPI_LONG_LONG, 0, TAG_HDR, comm, MPI_STATUS_IGNORE);
    c.n_local = (hyNode_index_t)hdr[0];
    c.n_oe = (hyEdge_index_t)hdr[1];
    c.n_owned = (hyNode_index_t)hdr[2];
    c.n_global = (hyNode_index_t)hdr[3];
    c.nprop = (unsigned int)hdr[4];

    c.lgmap.resize(c.n_local);
    c.pts.resize(c.n_local);
    c.edges.resize(c.n_oe);
    c.types.resize(c.n_oe);
    c.props.resize((size_t)c.n_oe * c.nprop);

    MPI_Recv(c.lgmap.data(), (int)(c.n_local * sizeof(hyNode_index_t)), MPI_BYTE, 0, TAG_LGMAP,
             comm, MPI_STATUS_IGNORE);
    MPI_Recv(c.pts.data(), (int)(c.n_local * sizeof(pointT)), MPI_BYTE, 0, TAG_PTS, comm,
             MPI_STATUS_IGNORE);
    MPI_Recv(c.edges.data(), (int)(c.n_oe * sizeof(node_pair_t)), MPI_BYTE, 0, TAG_EDGES, comm,
             MPI_STATUS_IGNORE);
    MPI_Recv(c.types.data(), (int)(c.n_oe * sizeof(node_pair_t)), MPI_BYTE, 0, TAG_TYPES, comm,
             MPI_STATUS_IGNORE);
    if (c.nprop)
      MPI_Recv(c.props.data(), (int)(c.props.size() * sizeof(value_t)), MPI_BYTE, 0, TAG_PROPS,
               comm, MPI_STATUS_IGNORE);

    return build_domain_info(c);
  }

  // ---- rank 0: read, partition, renumber ----
  DI full = read_domain_hdf5<hyEdge_dim, space_dim, vectorT, pointT, hyEdge_index_t, hyNode_index_t,
                             pt_index_t>(filename, /*serialize=*/false);
  const hyNode_index_t n = full.n_hyNodes;  // == n_points for hyEdge_dim == 1
  const hyEdge_index_t ne = full.n_hyEdges;
  const unsigned int nprop = full.n_properties;

  std::vector<int> part(n);
  for (hyNode_index_t v = 0; v < n; ++v)
    part[v] = compute_partition_placeholder(v, size);

  // Owned-contiguous global renumbering: stable-sort node ids by part (keeps id order within part).
  std::vector<hyNode_index_t> perm(n);  // perm[k] = old node id at new global position k
  std::iota(perm.begin(), perm.end(), hyNode_index_t{0});
  std::stable_sort(perm.begin(), perm.end(),
                   [&](hyNode_index_t a, hyNode_index_t b) { return part[a] < part[b]; });
  std::vector<hyNode_index_t> new_global(n);
  for (hyNode_index_t k = 0; k < n; ++k)
    new_global[perm[k]] = k;
  std::vector<hyNode_index_t> offset(size + 1, 0);  // node range start per rank
  for (hyNode_index_t v = 0; v < n; ++v)
    ++offset[part[v] + 1];
  for (int r = 0; r < size; ++r)
    offset[r + 1] += offset[r];

  // Bucket each edge to its owning rank = min part of its endpoints.
  std::vector<std::vector<hyEdge_index_t> > edges_of(size);
  for (hyEdge_index_t e = 0; e < ne; ++e)
  {
    const auto& nodes = full.hyNodes_hyEdge[e];
    const int owner = std::min(part[nodes[0]], part[nodes[1]]);
    edges_of[owner].push_back(e);
  }

  auto build_chunk = [&](const int r) -> Chunk
  {
    Chunk c;
    c.n_owned = offset[r + 1] - offset[r];
    c.n_global = n;
    c.nprop = nprop;
    const auto& myedges = edges_of[r];
    c.n_oe = (hyEdge_index_t)myedges.size();

    // Ghosts: endpoints of owned edges that belong to another rank, ordered by global index.
    std::vector<hyNode_index_t> ghosts;
    for (hyEdge_index_t e : myedges)
      for (int s = 0; s < 2; ++s)
      {
        const hyNode_index_t v = full.hyNodes_hyEdge[e][s];
        if (part[v] != r)
          ghosts.push_back(v);
      }
    std::sort(ghosts.begin(), ghosts.end(),
              [&](hyNode_index_t a, hyNode_index_t b) { return new_global[a] < new_global[b]; });
    ghosts.erase(std::unique(ghosts.begin(), ghosts.end()), ghosts.end());
    const hyNode_index_t n_ghost = (hyNode_index_t)ghosts.size();
    c.n_local = c.n_owned + n_ghost;

    std::unordered_map<hyNode_index_t, hyNode_index_t> ghost_local;
    for (hyNode_index_t g = 0; g < n_ghost; ++g)
      ghost_local[ghosts[g]] = c.n_owned + g;
    auto local_of = [&](const hyNode_index_t v) -> hyNode_index_t
    { return part[v] == r ? new_global[v] - offset[r] : ghost_local[v]; };

    c.lgmap.resize(c.n_local);
    c.pts.resize(c.n_local);
    for (hyNode_index_t l = 0; l < c.n_owned; ++l)
    {
      c.lgmap[l] = offset[r] + l;
      c.pts[l] = full.points[perm[offset[r] + l]];
    }
    for (hyNode_index_t g = 0; g < n_ghost; ++g)
    {
      c.lgmap[c.n_owned + g] = new_global[ghosts[g]];
      c.pts[c.n_owned + g] = full.points[ghosts[g]];
    }

    c.edges.resize(c.n_oe);
    c.types.resize(c.n_oe);
    c.props.resize((size_t)c.n_oe * nprop);
    for (hyEdge_index_t i = 0; i < c.n_oe; ++i)
    {
      const hyEdge_index_t e = myedges[i];
      c.edges[i] = node_pair_t{local_of(full.hyNodes_hyEdge[e][0]),
                               local_of(full.hyNodes_hyEdge[e][1])};
      c.types[i] = full.hyFaces_hyEdge[e];
      for (unsigned int p = 0; p < nprop; ++p)
        c.props[(size_t)i * nprop + p] = full.hyEdge_properties[e][p];
    }
    return c;
  };

  // Send chunks to ranks 1..size-1, keep rank 0's own.
  for (int r = 1; r < size; ++r)
  {
    Chunk c = build_chunk(r);
    long long hdr[5] = {(long long)c.n_local, (long long)c.n_oe, (long long)c.n_owned,
                        (long long)c.n_global, (long long)c.nprop};
    MPI_Send(hdr, 5, MPI_LONG_LONG, r, TAG_HDR, comm);
    MPI_Send(c.lgmap.data(), (int)(c.n_local * sizeof(hyNode_index_t)), MPI_BYTE, r, TAG_LGMAP,
             comm);
    MPI_Send(c.pts.data(), (int)(c.n_local * sizeof(pointT)), MPI_BYTE, r, TAG_PTS, comm);
    MPI_Send(c.edges.data(), (int)(c.n_oe * sizeof(node_pair_t)), MPI_BYTE, r, TAG_EDGES, comm);
    MPI_Send(c.types.data(), (int)(c.n_oe * sizeof(node_pair_t)), MPI_BYTE, r, TAG_TYPES, comm);
    if (c.nprop)
      MPI_Send(c.props.data(), (int)(c.props.size() * sizeof(value_t)), MPI_BYTE, r, TAG_PROPS,
               comm);
  }

  return build_domain_info(build_chunk(0));
}

#endif  // HYPERHDG_PETSC
