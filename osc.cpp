#include <cmath>
#include <vector>
#include <cstdlib>
#include <cstdio>

const double SIGMA = 10.0;
const double RHO   = 28.0;
const double BETA  = 8.0 / 3.0;
const double DT      = 0.001;
const int    N_STEPS = 200000;
const double X0 = 1.0, Y0 = 1.0, Z0 = 1.0;
const double D0 = 1e-8;

struct State {
    double x, y, z;
    State(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
};

State f(const State& s) {
    return State(SIGMA * (s.y - s.x), s.x * (RHO - s.z) - s.y, s.x * s.y - BETA * s.z);
}

State rk4_step(const State& s, double dt) {
    State k1 = f(s);
    State k2 = f(State(s.x + 0.5 * dt * k1.x, s.v + 0.5 * dt * k1.v));
    State k3 = f(State(s.x + 0.5 * dt * k2.x, s.v + 0.5 * dt * k2.v));
    State k4 = f(State(s.x + dt * k3.x,       s.v + dt * k3.v));

    State next;
    next.x = s.x + (dt / 6.0) * (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x);
    next.y = s.y + (dt / 6.0) * (k1.y + 2.0 * k2.y + 2.0 * k3.y + k4.y);
    next.z = s.z + (dt / 6.0) * (k1.z + 2.0 * k2.z + 2.0 * k3.z + k4.z);
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

double norm(const State& s) {
    return std::sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
}

struct LyapResult {
    double lambda;
    std::vector<double> t;
    std::vector<double> lambda_t;
    std::vector<State> trajectory;
};

LyapResult benettin(const State& s0, double d0, int n_steps) {
    State s = s0;
    State delta(1.0, 0.0, 0.0);
    double nd = norm(delta);
    delta.x /= nd; delta.y /= nd; delta.z /= nd;

    State sp(s.x + d0 * delta.x,
             s.y + d0 * delta.y,
             s.z + d0 * delta.z);

    double sum = 0.0;
    LyapResult res;
    res.t.reserve(n_steps);
    res.lambda_t.reserve(n_steps);
    res.trajectory.reserve(n_steps + 1);
    res.trajectory.push_back(s);

    for (int i = 1; i <= n_steps; ++i) {
        s  = rk4_step(s,  DT);
        sp = rk4_step(sp, DT);
        double dx = sp.x - s.x;
        double dy = sp.y - s.y;
        double dz = sp.z - s.z;
        double d  = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (d == 0.0) {
            dx = delta.x * d0;
            dy = delta.y * d0;
            dz = delta.z * d0;
            d  = d0;
        }
        
        sum += std::log(d / d0);
        
        sp.x = s.x + (dx / d) * d0;
        sp.y = s.y + (dy / d) * d0;
        sp.z = s.z + (dz / d) * d0;
        double time = i * DT;
        res.t.push_back(time);
        res.lambda_t.push_back(sum / time);
        res.trajectory.push_back(s);
    }

    res.lambda = sum / (n_steps * DT);
    return res;
}

int main(int argc, char** argv) {
    double x0 = X0, y0 = Y0, z0 = Z0;
    if (argc >= 4) {
        x0 = std::atof(argv[1]);
        y0 = std::atof(argv[2]);
        z0 = std::atof(argv[3]);
    }
    double d0 = (argc >= 5) ? std::atof(argv[4]) : D0;

    std::printf("Вычисление показателя Ляпунова методом Бенеттина\n");
    std::printf("sigma = %.3f, rho = %.3f, beta = %.3f\n", SIGMA, RHO, BETA);
    std::printf("x0 = %.3f, y0 = %.3f, z0 = %.3f\n", x0, y0, z0);
    std::printf("dt = %.4lf, шагов = %d, d0 = %.1le\n\n", DT, N_STEPS, d0);

    State s0(x0, y0, z0);
    LyapResult res = benettin(s0, d0, N_STEPS);

    dump_lyap(res);
    dump_traj(res.trajectory);

    return 0;
}
