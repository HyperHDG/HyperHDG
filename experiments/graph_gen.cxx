#include <stdio.h>
#include <petsc.h>
#include <petscviewerhdf5.h>
#include <assert.h>

PetscInt compute_vertex_type(const PetscReal* vertex, const PetscReal* max_p, const PetscReal* min_p) {
    if (vertex[0] - min_p[0] < 1e-6 * (max_p[0] - min_p[0]) || max_p[0] - vertex[0] < 1e-6 * (max_p[0] - min_p[0]) ||
        vertex[1] - min_p[1] < 1e-6 * (max_p[1] - min_p[1]) || max_p[1] - vertex[1] < 1e-6 * (max_p[1] - min_p[1]))
      return 1;
    else
      return 0;
}

void compute_types(const PetscReal* vertices, size_t n, const PetscInt* edges, size_t m, PetscInt* node_types, PetscInt* types) {
  // Calculate the bounding box (min/max x, y, z) for all vertices
  PetscReal min_p[3] = {1e10, 1e10, 1e10};
  PetscReal max_p[3] = {1e-10, 1e-10, 1e-10};
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < 3; j++) {
      min_p[j] = std::min(min_p[j], vertices[3*i+j]);
      max_p[j] = std::max(max_p[j], vertices[3*i+j]);
    }
  }
  for (size_t i = 0; i < n; i++)
    node_types[i] = compute_vertex_type(vertices+3*i, max_p, min_p);

  for (size_t i = 0; i < 2*m; i++)
    types[i] = node_types[edges[i]];
}

PetscErrorCode generate_grid_graph(PetscInt _n, Vec *vpoints, IS *is_edges, IS *is_types_nodes, IS *is_types_faces) {
  size_t k = 0, n = _n, m = 2*(n-1)*n;
  double h = 1. / (n-1);
  PetscInt *edges, *types_nodes, *types_faces;
  PetscReal *points;
  PetscCall(VecCreateSeq(PETSC_COMM_SELF, 3*n*n, vpoints));
  PetscCall(VecSetBlockSize(*vpoints, 3));
  PetscCall(PetscObjectSetName((PetscObject)*vpoints, "points"));
  PetscCall(PetscMalloc3(2*m, &edges, n*n, &types_nodes, 2*m, &types_faces));

  PetscCall(VecGetArray(*vpoints, &points));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      points[3*(i*n+j)+0] = j*h; // x
      points[3*(i*n+j)+1] = i*h; // y
      points[3*(i*n+j)+2] = 0;   // z

      if (j+1 < n) {
        edges[k]   = i*n+j;
        edges[k+1] = i*n+j+1;
        k += 2;
      }
      if (i+1 < n) {
        edges[k]   = i*n+j;
        edges[k+1] = (i+1)*n+j;
        k += 2;
      }
    }
  }
  assert(k == 2*m);
  PetscCall(ISCreateGeneral(PETSC_COMM_SELF, 2*m, edges, PETSC_COPY_VALUES, is_edges));
  PetscCall(ISSetBlockSize(*is_edges, 2));
  PetscCall(PetscObjectSetName((PetscObject)*is_edges, "edges"));

  compute_types(points, n*n, edges, m, types_nodes, types_faces);
  PetscCall(ISCreateGeneral(PETSC_COMM_SELF, n*n, types_nodes, PETSC_COPY_VALUES, is_types_nodes));
  PetscCall(PetscObjectSetName((PetscObject)*is_types_nodes, "types_nodes"));
  PetscCall(ISCreateGeneral(PETSC_COMM_SELF, 2*m, types_faces, PETSC_COPY_VALUES, is_types_faces));
  PetscCall(PetscObjectSetName((PetscObject)*is_types_faces, "types_faces"));
  PetscCall(ISSetBlockSize(*is_types_faces, 2));
  
  PetscCall(VecRestoreArray(*vpoints, &points));
  PetscCall(PetscFree(edges));

  return 0;
}

int main(int argc, char** argv) {
  char out[PATH_MAX] = "graph";
  PetscBool is_set;
  PetscInt n = 2;
  Vec points;
  IS edges, ntypes, ftypes;
  PetscViewer viewer;
  char buf[PATH_MAX];

  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "graph_gen options", NULL);
  PetscCall(PetscOptionsInt("-n", "n", NULL, n, &n, &is_set));
  PetscCall(PetscOptionsString("-o", "output path", NULL, out, out, PATH_MAX, &is_set));
  PetscOptionsEnd();

  snprintf(buf, sizeof(buf), "%s.geo.h5", out);
  PetscCall(generate_grid_graph(n, &points, &edges, &ntypes, &ftypes));
  PetscCall(PetscViewerHDF5Open(PETSC_COMM_SELF, buf, FILE_MODE_WRITE, &viewer));
  PetscCall(PetscViewerHDF5PushGroup(viewer, "/domain"));
  PetscCall(VecView(points, viewer));
  PetscCall(ISView(edges, viewer));
  PetscCall(ISView(ntypes, viewer));
  PetscCall(ISView(ftypes, viewer));

  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(VecDestroy(&points));
  PetscCall(ISDestroy(&edges));

  PetscCall(PetscFinalize());
}
