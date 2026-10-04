#include <cmath>
#include <vector>
#include <cstdio>

const double MASS   = 1.0;
const double K      = 1.0;
const double OMEGA  = std::sqrt(K / MASS);
const double T      = 2.0 * M_PI / OMEGA;
const double DT     = T / 30.0;
const int    N_STEPS = 300;
const double X0 = 1.0, V0 = 0.0;

struct State {
    double x, v;
    State(double x = 0, double v = 0) : x(x), v(v) {}
};

State f(const State& s) {
    return State(s.v, -(K / MASS) * s.x);
}

State rk4_step(const State& s, double dt) {
    State k1 = f(s);
    State k2 = f(State(s.x + 0.5 * dt * k1.x, s.v + 0.5 * dt * k1.v));
    State k3 = f(State(s.x + 0.5 * dt * k2.x, s.v + 0.5 * dt * k2.v));
    State k4 = f(State(s.x + dt * k3.x,       s.v + dt * k3.v));

    State next;
    next.x = s.x + (dt / 6.0) * (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x);
    next.v = s.v + (dt / 6.0) * (k1.v + 2.0 * k2.v + 2.0 * k3.v + k4.v);
    return next;
}

std::vector<State> solve_ode() {
    std::vector<State> trajectory;
    State current(X0, V0);
    trajectory.push_back(current);

    for (int i = 0; i < N_STEPS; ++i) {
        current = rk4_step(current, DT);
        trajectory.push_back(current);
    }
    return trajectory;
}

void dump_data(const std::vector<State>& traj) {
    FILE* fp = fopen("oscdata", "w");
    if (!fp) {
        printf("Не удалось открыть файл oscdata для записи\n");
    }

    for (size_t i = 0; i < traj.size(); ++i) {
        double t = i * DT;

        printf("%lf   %le %le\n", t, traj[i].x, traj[i].v);

        if (fp) {
            fprintf(fp, "%lf   %le %le\n", t, traj[i].x, traj[i].v);
        }
    }

    if (fp) fclose(fp);
}

int main() {
    printf("Данные гармонического осциллятора\n");
    printf("m = %.3f, k = %.3f, omega = %.3f, T = %.3f, dt = %.3f\n",
           MASS, K, OMEGA, T, DT);
    printf("x0 = %.3f, v0 = %.3f, шагов: %d\n", X0, V0, N_STEPS);

    std::vector<State> traj_rk4 = solve_ode();
    dump_data(traj_rk4);

    return 0;
}
