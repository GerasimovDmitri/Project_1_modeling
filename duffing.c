#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 5000
#define DT 0.01

static struct {
    double delta;
    double alpha;
    double beta;
    double F;
    double omega;
} params = {0.0, 1.0, 0.5, 0.0, 1.0}; 

static struct {
    double *t;
    double *x;
    double *v;
    double *E;
    int n;
} trajectory = {NULL, NULL, NULL, NULL, 0};

static void rhs(double t, double x, double v, double *dx, double *dv) {
    double force = params.F * cos(params.omega * t);
    double nonlinear = params.alpha * x + params.beta * x * x * x;
    
    *dx = v;
    *dv = -params.delta * v - nonlinear + force;
}

static void rk4_step(double *t, double *x, double *v) {
    double k1x, k1v, k2x, k2v, k3x, k3v, k4x, k4v;
    double dx, dv, t_half;
    
    rhs(*t, *x, *v, &dx, &dv);
    k1x = DT * dx;
    k1v = DT * dv;
    
    t_half = *t + DT/2.0;
    rhs(t_half, *x + k1x/2.0, *v + k1v/2.0, &dx, &dv);
    k2x = DT * dx;
    k2v = DT * dv;
    
    rhs(t_half, *x + k2x/2.0, *v + k2v/2.0, &dx, &dv);
    k3x = DT * dx;
    k3v = DT * dv;
    
    rhs(*t + DT, *x + k3x, *v + k3v, &dx, &dv);
    k4x = DT * dx;
    k4v = DT * dv;
    
    *x = *x + (k1x + 2.0*k2x + 2.0*k3x + k4x) / 6.0;
    *v = *v + (k1v + 2.0*k2v + 2.0*k3v + k4v) / 6.0;
    *t = *t + DT;
}

static double compute_energy(double x, double v) {
    double kinetic = 0.5 * v * v;
    double potential = 0.5 * params.alpha * x * x 
                     + 0.25 * params.beta * x * x * x * x;
    return kinetic + potential;
}

static int allocate_trajectory(int n) {
    trajectory.t = (double*)malloc(n * sizeof(double));
    trajectory.x = (double*)malloc(n * sizeof(double));
    trajectory.v = (double*)malloc(n * sizeof(double));
    trajectory.E = (double*)malloc(n * sizeof(double));
    
    if (!trajectory.t || !trajectory.x || !trajectory.v || !trajectory.E) {
        fprintf(stderr, "Ошибка памяти\n");
        return 0;
    }
    trajectory.n = n;
    return 1;
}

static void free_trajectory(void) {
    free(trajectory.t);
    free(trajectory.x);
    free(trajectory.v);
    free(trajectory.E);
    trajectory.n = 0;
}

static void integrate(void) {
    int i;
    trajectory.t[0] = 0.0;
    trajectory.x[0] = 1.0;
    trajectory.v[0] = 0.0;
    trajectory.E[0] = compute_energy(trajectory.x[0], trajectory.v[0]);
    for (i = 0; i < trajectory.n - 1; i++) {
        trajectory.t[i+1] = trajectory.t[i];
        trajectory.x[i+1] = trajectory.x[i];
        trajectory.v[i+1] = trajectory.v[i];
        
        rk4_step(&trajectory.t[i+1], &trajectory.x[i+1], &trajectory.v[i+1]);
        trajectory.E[i+1] = compute_energy(trajectory.x[i+1], trajectory.v[i+1]);
    }
}
static void save_results(const char *filename)
{
    FILE *fp = fopen(filename, "w");
    int i;
    double E0 = trajectory.E[0];
    
    if (!fp) {
        fprintf(stderr, "Файл не создан %s\n", filename);
        return;
    }
    
    fprintf(fp, "# t x v energy rel_error\n");
    
    for (i = 0; i < trajectory.n; i++) {
        double rel_err = fabs((trajectory.E[i] - E0) / E0);
        fprintf(fp, "%12.6f %12.6f %12.6f %12.6e %12.6e\n",
                trajectory.t[i], trajectory.x[i], 
                trajectory.v[i], trajectory.E[i], rel_err);
    }
    
    fclose(fp);
}
static void print_stats(void) {
    int i;
    double max_err = 0.0;
    double E0 = trajectory.E[0];
    
    for (i = 0; i < trajectory.n; i++) {
        double err = fabs((trajectory.E[i] - E0) / E0);
        if (err > max_err) max_err = err;
    }
    printf("  Максимальная ошибка энергии: %.2e\n", max_err);
    printf("  (ожидается < 1e-11 для RK4)\n");
    printf("\n");
}

static void print_params(void) {
    printf("Параметры:\n");
    printf("Затухание  %.3f\n", params.delta);
    printf("Линейная жсткость   %.3f\n", params.alpha);
    printf("Нелинейная жесткость %.3f\n", params.beta);
    printf("Аплитуда  %.3f\n", params.F);
    printf("Частота   %.3f\n", params.omega);
    printf("Численные параметры:\n");
    printf("dt = %.4f, шагов = %d\n", DT, N);
}

int main(void) {
    print_params();
    
    if (!allocate_trajectory(N)) {
        return 1;
    }
    
    integrate();
    save_results("duffing_data.txt");
    print_stats();
    free_trajectory();
    
    return 0;
}
