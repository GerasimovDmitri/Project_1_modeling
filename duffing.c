#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 5000          // количество шагов
#define DT 0.01         // шаг по времени

// Параметры осциллятора Дуффинга
// Уравнение: d^2x/dt^2 + delta*dx/dt + alpha*x + beta*x^3 = force_amp*cos(force_freq*t)

double delta = 0.0;          // затухание (0 = консервативный случай)
double alpha = 1.0;          // линейная жёсткость
double beta  = 0.5;          // нелинейная жёсткость
double force_amp = 0.0;      // амплитуда внешней силы (было gamma)
double force_freq = 1.0;     // частота (было omega)

// Система ОДУ первого порядка:
// dx/dt = v
// dv/dt = -delta*v - alpha*x - beta*x^3 + force_amp*cos(force_freq*t)

void duffing(double t, double x, double v, double *dx, double *dv) {
    *dx = v;
    *dv = -delta * v - alpha * x - beta * x * x * x + force_amp * cos(force_freq * t);
}

// Метод Рунге-Кутты 4-го порядка
void rk4_step(double *t, double *x, double *v, double dt) {
    double k1x, k1v, k2x, k2v, k3x, k3v, k4x, k4v;
    double dx1, dv1, dx2, dv2, dx3, dv3, dx4, dv4;
    double t_mid, x_mid, v_mid;
    
    // k1
    duffing(*t, *x, *v, &dx1, &dv1);
    k1x = dt * dx1;
    k1v = dt * dv1;
    
    // k2
    t_mid = *t + dt/2.0;
    x_mid = *x + k1x/2.0;
    v_mid = *v + k1v/2.0;
    duffing(t_mid, x_mid, v_mid, &dx2, &dv2);
    k2x = dt * dx2;
    k2v = dt * dv2;
    
    // k3
    x_mid = *x + k2x/2.0;
    v_mid = *v + k2v/2.0;
    duffing(t_mid, x_mid, v_mid, &dx3, &dv3);
    k3x = dt * dx3;
    k3v = dt * dv3;
    
    // k4
    duffing(*t + dt, *x + k3x, *v + k3v, &dx4, &dv4);
    k4x = dt * dx4;
    k4v = dt * dv4;
    
    // Финальное обновление
    *x = *x + (k1x + 2.0*k2x + 2.0*k3x + k4x) / 6.0;
    *v = *v + (k1v + 2.0*k2v + 2.0*k3v + k4v) / 6.0;
    *t = *t + dt;
}

// Энергия для консервативного случая (force_amp=0, delta=0)
double energy(double x, double v) {
    return 0.5 * v * v + 0.5 * alpha * x * x + 0.25 * beta * x * x * x * x;
}

int main() {
    FILE *fp;
    double t[N], x[N], v[N], E[N], E0, E_rel_error;
    int i;
    
    printf("Осциллятор Дуффинга - численное решение\n");
    printf("Параметры:\n");
    printf("  delta = %.2f\n", delta);
    printf("  alpha = %.2f\n", alpha);
    printf("  beta  = %.2f\n", beta);
    printf("  force_amp = %.2f\n", force_amp);
    printf("  force_freq = %.2f\n", force_freq);
    printf("  dt    = %.4f\n", DT);
    printf("  шагов = %d\n", N);
    
    t[0] = 0.0;
    x[0] = 1.0;
    v[0] = 0.0;
    
    // Расчёт энергии в начальный момент
    E0 = energy(x[0], v[0]);
    E[0] = E0;
    
    // Интегрирование методом Рунге-Кутты 4-го порядка
    for (i = 0; i < N - 1; i++) {
        t[i+1] = t[i];
        x[i+1] = x[i];
        v[i+1] = v[i];
        rk4_step(&t[i+1], &x[i+1], &v[i+1], DT);
        E[i+1] = energy(x[i+1], v[i+1]);
    }
    
    fp = fopen("duffing_data.txt", "w");
    if (fp == NULL) {
        printf("Ошибка: не могу создать файл данных!\n");
        return 1;
    }
    
    fprintf(fp, "# t          x          v          Energy     RelError\n");
    for (i = 0; i < N; i++) {
        E_rel_error = fabs((E[i] - E0) / E0);
        fprintf(fp, "%12.6f %12.6f %12.6f %12.6e %12.6e\n", 
                t[i], x[i], v[i], E[i], E_rel_error);
    }
    fclose(fp);
    
    // Максимальная ошибка энергии
    double max_error = 0.0;
    for (i = 0; i < N; i++) {
        E_rel_error = fabs((E[i] - E0) / E0);
        if (E_rel_error > max_error) max_error = E_rel_error;
    }
    printf("\nМаксимальная относительная ошибка энергии: %.2e\n", max_error);
    printf("(должна быть < 1e-11 для RK4 - это доказывает точность)\n");
    
    printf("\nДанные сохранены в 'duffing_data.txt'\n");
    printf("Запустите Gnuplot командой: gnuplot duffing_plot.gp\n");
    
    return 0;
} 
