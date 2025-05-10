#include "model.h"
#include "devoir_3.h"
#include "devoir_2.h"

void print_csr_dense(const CSRMatrix *A) {
    for (int i = 0; i < A->n; ++i) {
        int row_start = A->row_ptr[i];
        int row_end = A->row_ptr[i + 1];

        int col_idx = 0;
        for (int k = row_start; k < row_end; ++k) {
            int col = A->col_idx[k];
            while (col_idx < col) {
                printf(" 0.00 ");
                ++col_idx;
            }
            printf("%5.2f ", A->data[k]);
            ++col_idx;
        }
        while (col_idx < A->n) {
            printf(" 0.00 ");
            ++col_idx;
        }
        printf("\n");
    }
}
State* newMark(State *state_0, double T, double dt, double beta, double gamma, CSRMatrix *Ksp, CSRMatrix *Msp, int I){
    int numberSteps = (int)(T/dt);

    printf("KSP\n\n\n");
    print_csr_dense(Ksp);
    printf("MSP\n\n\n");
    print_csr_dense(Msp);



    printf("n = %d\n", Ksp->n);

    double *q_old = (double *)malloc(Ksp->n * sizeof(double));
    double *q_new = (double *)malloc(Ksp->n * sizeof(double));
    double *p_old = (double *)malloc(Ksp->n * sizeof(double));
    double *p_new = (double *)malloc(Ksp->n * sizeof(double));
    double *v = (double *)malloc(Ksp->n * sizeof(double));

    for (int i = 0; i < Ksp->n; i++){
        q_old[i] = state_0->u[i];
        v[i] = state_0->v[i];
    }

    Matvec(Msp->n, Msp->row_ptr, Msp->col_idx, Msp->data, state_0->v, p_old);


    // Calcul matrice A equation 1 et matrice du membre de droite
    CSRMatrix *Aq = malloc(sizeof(CSRMatrix));
    Aq->n = Ksp->n;
    Aq->nnz = Ksp->nnz;
    Aq->row_ptr = malloc((Ksp->n + 1) * sizeof(int));
    Aq->row_ptr = memcpy(Aq->row_ptr, Ksp->row_ptr, (Ksp->n + 1) * sizeof(int));
    Aq->col_idx = malloc(Aq->nnz * sizeof(int));
    Aq->col_idx = memcpy(Aq->col_idx, Ksp->col_idx, Aq->nnz * sizeof(int));
    Aq->data = malloc(Aq->nnz * sizeof(double));
    Aq->data = memcpy(Aq->data, Ksp->data, Aq->nnz * sizeof(double));
    for(int i = 0; i < Aq->nnz; i++){
        Aq->data[i] *= dt*dt*beta;
    } 
    CSRMatrix *bq_partial = malloc(sizeof(CSRMatrix));
    bq_partial->n = Ksp->n;
    bq_partial->nnz = Ksp->nnz;
    bq_partial->row_ptr = malloc((Ksp->n + 1) * sizeof(int));
    bq_partial->row_ptr = memcpy(bq_partial->row_ptr, Ksp->row_ptr, (Ksp->n + 1) * sizeof(int));
    bq_partial->col_idx = malloc(bq_partial->nnz * sizeof(int));
    bq_partial->col_idx = memcpy(bq_partial->col_idx, Ksp->col_idx, bq_partial->nnz * sizeof(int));
    bq_partial->data = malloc(bq_partial->nnz * sizeof(double));
    bq_partial->data = memcpy(bq_partial->data, Ksp->data, bq_partial->nnz * sizeof(double));

    for(int i = 0; i < bq_partial->nnz; i++){
        bq_partial->data[i] *= -0.5*dt*dt*(1-2*beta);
    }

    for(int row = 0; row < Aq->n; row++){
        int colsIdxAq = Aq->row_ptr[row];
        for(int j = Msp->row_ptr[row]; j < Msp->row_ptr[row + 1]; j++){
            while(Msp->col_idx[j] != Aq->col_idx[colsIdxAq]){
                colsIdxAq++;
            }
            Aq->data[colsIdxAq] += Msp->data[j];
            bq_partial->data[colsIdxAq] += Msp->data[j];
        }
    }
    printf("Aq\n\n\n");
    print_csr_dense(Aq);
    printf("bq_partial\n\n\n");
    print_csr_dense(bq_partial);

    // Open the file to write the results
    FILE *file = fopen("time.txt", "w");
    if (file == NULL) {
        printf("Error opening file\n");
        return NULL;
    }
    fprintf(file, "%.15le %.15le %.15le %.15le %.15le\n", 0.0, state_0->u[2 * I], state_0->u[2 * I + 1], v[2 * I], v[2 * I + 1]);


    for(int i = 0 ; i<numberSteps; i++){
        double *b = (double *)malloc(Ksp->n * sizeof(double));
        double *b_2 = (double *)malloc(Ksp->n * sizeof(double));
        double *b_3 = (double *)malloc(Ksp->n * sizeof(double));
        Matvec(Aq->n, bq_partial->row_ptr, bq_partial->col_idx, bq_partial->data, q_old, b);
        for(int j = 0; j < Aq->n; j++){
            b[j] += dt*p_old[j];
        }

        CG(Aq->n, Aq->row_ptr, Aq->col_idx, Aq->data, b, q_new, 1e-12);

        // Write the current time step and node I's displacement and velocity to the file

        printf("Iteration %d\n : u[0] : (%f, %f)", i, q_new[0], q_new[1]);


        for(int j = 0; j < Aq->n; j++){
            b_2[j] = (1-gamma)*q_old[j] + gamma*q_new[j];
        }

        Matvec(Ksp->n, Ksp->row_ptr, Ksp->col_idx, Ksp->data, b_2, b_3);

        for(int j = 0; j < Aq->n; j++){
            p_new[j] = p_old[j] - dt*b_3[j];
        }

        for(int j = 0; j < Aq->n; j++){
            q_old[j] = q_new[j];
            p_old[j] = p_new[j];
        }

        CG(Msp->n, Msp->row_ptr, Msp->col_idx, Msp->data, p_new, v, 1e-12);

        fprintf(file, "%.15le %.15le %.15le %.15le %.15le\n", (i+1) * dt, q_new[2 * I], q_new[2 * I + 1], v[2 * I], v[2 * I + 1]);
        
        
        if(1 == 1){
            char output_filename[50];
            snprintf(output_filename, sizeof(output_filename), "./plots/data/time%d.txt", i);
            FILE *output_file = fopen(output_filename, "w");
            for (int l = 0; l < ((int)Ksp->n/2); l++) {
                fprintf(output_file, "%.15le %.15le %.15le %.15le\n", 
                    q_new[2 * l], q_new[2 * l + 1], v[2 * l], v[2 * l + 1]);
            }
            fclose(output_file);
        }
        
        
        free(b);
        free(b_2);
        free(b_3);
    }

    // Close the file
    fclose(file);

    state_0->u = q_new;
    state_0->v = v;
    return state_0;
}