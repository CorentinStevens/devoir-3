#include "devoir_2.h"
#include "devoir_3.h"
#include "utils.h"
#include "model.h"
#include "utils_gmsh.h"
#include <math.h>
#include <cblas.h>
#include <gmshc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define VERBOSE 1
#define PRECISION 10

void display_sol(FE_Model *model, double *sol) {
    int ierr, n_views, *views;
    double *bounds;
    add_gmsh_views(&views, &n_views, &bounds);

    double *data_forces = malloc(6 * model->n_bd_edge * sizeof(double));
    visualize_disp(model, sol, views[1], 0, &bounds[2]);
    visualize_stress(model, sol, views, 1, 0, data_forces, bounds);
    visualize_bd_forces(model, data_forces, views[0], 1, &bounds[0]);

    create_tensor_aliases(views);
    set_view_options(n_views, views, bounds);
    gmshFltkRun(&ierr);
    gmshFltkFinalize(&ierr);
    free(data_forces);
}

void display_info(FE_Model *model, int step, struct timespec ts[4]) {

    char *m_str[3] = {"Plane stress", "Plane strain", "Axisymmetric"};
    char *r_str[4] = {"No", "X", "Y", "RCMK"};

    if (step == 1) {
        printf(
            "\n===========  Linear elasticity simulation - FEM  ===========\n\n"
        );
        printf("%30s = %s\n", "Model", model->model_name);
        printf("%30s = %s\n", "Model type", m_str[model->m_type]);
        printf("%30s = %.3e\n", "Young's Modulus E", model->E);
        printf("%30s = %.3e\n", "Poisson ratio nu", model->nu);
        printf("%30s = %.3e\n\n", "Density rho", model->rho);
    } else if (step == 2) {
        char *e_str = (model->e_type == TRI) ? "Triangle" : "Quadrilateral";
        printf("%30s = %s\n", "Element type", e_str);
        printf("%30s = %zu\n", "Number of elements", model->n_elem);
        printf("%30s = %zu\n", "Number of nodes", model->n_node);
        printf("%30s = %s\n", "Renumbering", r_str[model->renum]);
        printf("%30s = %zu\n\n", "Matrix bandwidth", 2 * model->node_band + 1);
    }
}


int main(int argc, char *argv[]) {

    int ierr;
    double mesh_size_ratio;
    if ((argc < 3) || (sscanf(argv[2], "%lf", &mesh_size_ratio)) != 1) {
        printf("Usage: \n./deformation <model> <mesh_size_ratio>\n");
        printf("model: one of the model implemented in models/\n");
        printf("mesh_size_ratio: mesh size factor\n");
        return -1;
    }

    // Simulation parameters
    const ElementType e_type = TRI;
    const Renumbering renum = RENUM_NO;  // let gmsh do the RCMK renumbering

    FE_Model *model = create_FE_Model(argv[1], e_type, renum);
    display_info(model, 1, NULL);

    gmshInitialize(argc, argv, 0, 0, &ierr);
    gmshOptionSetNumber("General.Verbosity", 2, &ierr);
    model->mesh_model(mesh_size_ratio, e_type);

    load_mesh(model);
    renumber_nodes(model);
    display_info(model, 2, NULL);
    assemble_system(model);
    double *rhs = (double *)calloc(2 * model->n_node, sizeof(double));
    double *sol = (double *)calloc(2 * model->n_node, sizeof(double));
    add_bulk_source(model, rhs);
    enforce_bd_conditions(model, rhs);
    

    double *coord = model->coords;

    FILE *coord_file = fopen("modelCoordinates.txt", "w");
    if (!coord_file) {
        perror("Failed to open modelCoordinates.txt");
        return -1;
    }

    for (size_t i = 0; i < model->n_node; i++) {
        fprintf(coord_file, "%.15le, %.15le\n", coord[2 * i], coord[2 * i + 1]);
    }

    fclose(coord_file);
    printf("Model coordinates written to modelCoordinates.txt\n");

    

    // Write edges to a file in the format (node1, node2)
    FILE *edges_file = fopen("edges.txt", "w");
    if (!edges_file) {
        perror("Failed to open edges.txt");
        return -1;
    }

    for (size_t i = 0; i < model->n_bd_edge; i++) {
        size_t node1 = model->bd_edges[4 * i];
        size_t node2 = model->bd_edges[4 * i + 1];
        fprintf(edges_file, "%zu, %zu\n", node1, node2);
    }

    fclose(edges_file);
    printf("Model edges written to edges.txt\n");
    
    SymBandMatrix *Kbd = model->K;
    SymBandMatrix *Mbd = model->M;
    CSRMatrix *Ksp = band_to_sym_csr(Kbd);
    CSRMatrix *Msp = band_to_sym_csr(Mbd);
    double eps = 1e-8;


    const char* filename = "initial_fork_1.0.txt";
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        return 1;
    }

    // Compter le nombre de lignes
    int lines = 0;
    double a, b, c, d;
    while (fscanf(file, "%lf %lf %lf %lf", &a, &b, &c, &d) == 4) {
        lines++;
    }

    // Revenir au début du fichier
    rewind(file);

    // Allouer dynamiquement u et v
    double* u = malloc(2 * lines * sizeof(double));
    double* v = malloc(2 * lines * sizeof(double));
    if (!u || !v) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(file);
        return 1;
    }

    // Lire à nouveau et remplir u et v
    int index = 0;
    while (fscanf(file, "%lf %lf %lf %lf", &a, &b, &c, &d) == 4) {
        u[2 * index]     = a;
        u[2 * index + 1] = b;
        v[2 * index]     = c;
        v[2 * index + 1] = d;
        index++;
    }

    fclose(file);

    State *state_0 = malloc(sizeof(State));
    state_0->u = u;
    state_0->v = v;

    // Print the contents of u and v
    printf("Contents of u:\n");
    for (int i = 0; i < lines; i++) {
        printf("u[%d] = (%lf, %lf)\n", i, u[2 * i], u[2 * i + 1]);
    }

    printf("\nContents of v:\n");
    for (int i = 0; i < lines; i++) {
        printf("v[%d] = (%lf, %lf)\n", i, v[2 * i], v[2 * i + 1]);
    }

    int I = 0;
     
    state_0 = newMark(state_0, 20.0,0.001,0.25,0.5, Ksp, Msp, I);

    const char *output_filename = "final.txt";
    FILE *output_file = fopen(output_filename, "w");
    if (!output_file) {
        perror("Failed to open output file");
        free(state_0->u);
        free(state_0->v);
        free(state_0);
        return 1;
    }

    for (int i = 0; i < lines; i++) {
        fprintf(output_file, "%.15le %.15le %.15le %.15le\n", 
                state_0->u[2 * i], state_0->u[2 * i + 1], 
                state_0->v[2 * i], state_0->v[2 * i + 1]);
    }

    fclose(output_file);
    printf("Final state written to %s\n", output_filename);



    double *newSol = state_0->u;

    CG(Ksp->n, Ksp->row_ptr, Ksp->col_idx, Ksp->data, rhs, sol, eps);    
    display_sol(model, newSol);
    
    // Free stuff
    free_csr(Ksp);
    free_csr(Msp);
    gmshFinalize(&ierr);
    free(sol);
    free(rhs);
    free_FE_Model(model);
    return 0;
}
