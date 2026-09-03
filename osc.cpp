#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cmath>
#include <vector>
#include <cstdio>

const double OMEGA = 1.0;
const double T = 2.0 * M_PI / OMEGA;
const double DT = T / 30.0;
const int N_STEPS = 300;
const double X0 = 1.0, V0 = 0.0;

struct State {
    double x, v;
    State(double x = 0, double v = 0) : x(x), v(v) {}
};

struct Point2D {
    double x, y;
    Point2D(double x = 0, double y = 0) : x(x), y(y) {}
};

double analytical_x(double t) {
    return X0 * cos(OMEGA * t) + (V0 / OMEGA) * sin(OMEGA * t);
}

double analytical_v(double t) {
    return -X0 * OMEGA * sin(OMEGA * t) + V0 * cos(OMEGA * t);
}

State euler_explicit(const State& s, double dt) {
    double a = -OMEGA * OMEGA * s.x;
    return State(s.x + s.v * dt, s.v + a * dt);
}

State euler_implicit(const State& s, double dt) {
    double x_new = (s.x + s.v * dt) / (1.0 + (OMEGA * dt) * (OMEGA * dt));
    double v_new = s.v - OMEGA * OMEGA * x_new * dt;
    return State(x_new, v_new);
}

State verlet(const State& s, double dt) {
    double v_half = s.v - 0.5 * OMEGA * OMEGA * s.x * dt;
    double x_new = s.x + v_half * dt;
    double v_new = v_half - 0.5 * OMEGA * OMEGA * x_new * dt;
    return State(x_new, v_new);
}

std::vector<State> solve_ode(State (*method)(const State&, double)) {
    std::vector<State> trajectory;
    State current(X0, V0);
    trajectory.push_back(current);
    
    for (int i = 0; i < N_STEPS; ++i) {
        current = method(current, DT);
        trajectory.push_back(current);
    }
    return trajectory;
}

std::vector<State> traj_euler, traj_implicit, traj_verlet;
std::vector<Point2D> analytic_traj;
int current_frame = 0;
bool is_animating = true;

void draw_text(float x, float y, const char* text) {
    glColor3f(1.0, 1.0, 1.0);
    glRasterPos2f(x, y);
    for (const char* c = text; *c != '\0'; ++c) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }
}

void draw_phase_portrait() {
    glViewport(0, 0, 400, 400);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-2.0, 2.0, -2.0, 2.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glColor3f(0.2, 0.2, 0.2);
    glBegin(GL_LINES);
    for (double x = -1.5; x <= 1.5; x += 0.5) {
        glVertex2f(x, -1.5); glVertex2f(x, 1.5);
        glVertex2f(-1.5, x); glVertex2f(1.5, x);
    }
    glEnd();
    
    glColor3f(0.5, 0.5, 0.5);
    glBegin(GL_LINES);
    glVertex2f(-2.0, 0); glVertex2f(2.0, 0);
    glVertex2f(0, -2.0); glVertex2f(0, 2.0);
    glEnd();
    glColor3f(0.7, 0.7, 0.7);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < 100; ++i) {
        double angle = 2.0 * M_PI * i / 100.0;
        double x = X0 * cos(angle);
        double v = -X0 * OMEGA * sin(angle);
        glVertex2f(x, v);
    }
    glEnd();
    
    auto draw_trajectory = [&](const std::vector<State>& traj, float r, float g, float b, int max_idx) {
        glColor3f(r, g, b);
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= max_idx && i < (int)traj.size(); ++i) {
            glVertex2f(traj[i].x, traj[i].v);
        }
        glEnd();
        
        if (max_idx < (int)traj.size()) {
            glPointSize(8.0);
            glBegin(GL_POINTS);
            glVertex2f(traj[max_idx].x, traj[max_idx].v);
            glEnd();
        }
    };
    
    int idx = std::min(current_frame, N_STEPS);
    draw_trajectory(traj_euler, 1.0, 0.2, 0.2, idx);
    draw_trajectory(traj_implicit, 0.2, 0.2, 1.0, idx);
    draw_trajectory(traj_verlet, 0.2, 1.0, 0.2, idx);
    draw_text(-1.9, 1.8, "Phase Portrait");
    glColor3f(1.0, 0.2, 0.2); draw_text(-1.9, 1.6, "Euler (explicit)");
    glColor3f(0.2, 0.2, 1.0); draw_text(-1.9, 1.45, "Euler (implicit)");
    glColor3f(0.2, 1.0, 0.2); draw_text(-1.9, 1.3, "Verlet");
}

void draw_position_vs_time() {
    glViewport(400, 0, 400, 400);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, N_STEPS * DT, -2.0, 2.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glColor3f(0.2, 0.2, 0.2);
    glBegin(GL_LINES);
    for (double t = 0; t <= N_STEPS * DT; t += T/4) {
        glVertex2f(t, -1.5); glVertex2f(t, 1.5);
    }
    for (double x = -1.5; x <= 1.5; x += 0.5) {
        glVertex2f(0, x); glVertex2f(N_STEPS * DT, x);
    }
    glEnd();
    
    glColor3f(0.7, 0.7, 0.7);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i <= current_frame && i <= N_STEPS; ++i) {
        double t = i * DT;
        glVertex2f(t, analytical_x(t));
    }
    glEnd();
    
    auto draw_time_series = [&](const std::vector<State>& traj, float r, float g, float b) {
        glColor3f(r, g, b);
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= current_frame && i < (int)traj.size(); ++i) {
            glVertex2f(i * DT, traj[i].x);
        }
        glEnd();
    };
    
    draw_time_series(traj_euler, 1.0, 0.2, 0.2);
    draw_time_series(traj_implicit, 0.2, 0.2, 1.0);
    draw_time_series(traj_verlet, 0.2, 1.0, 0.2);
    draw_text(0.5, 1.8, "Position x(t)");
    glColor3f(1.0, 1.0, 1.0);
    char time_str[50];
    sprintf(time_str, "t = %.2f", current_frame * DT);
    draw_text(0.5, 1.6, time_str);
}

void draw_energy() {
    glViewport(0, 400, 400, 400);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, N_STEPS * DT, 0.0, 1.2);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glColor3f(0.2, 0.2, 0.2);
    glBegin(GL_LINES);
    for (double t = 0; t <= N_STEPS * DT; t += T/4) {
        glVertex2f(t, 0.0); glVertex2f(t, 1.2);
    }
    for (double E = 0.2; E <= 1.0; E += 0.2) {
        glVertex2f(0, E); glVertex2f(N_STEPS * DT, E);
    }
    glEnd();
    double E_exact = 0.5 * (V0*V0 + OMEGA*OMEGA * X0*X0);
    glColor3f(0.7, 0.7, 0.7);
    glBegin(GL_LINES);
    glVertex2f(0, E_exact); glVertex2f(N_STEPS * DT, E_exact);
    glEnd();
    auto draw_energy_series = [&](const std::vector<State>& traj, float r, float g, float b) {
        glColor3f(r, g, b);
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= current_frame && i < (int)traj.size(); ++i) {
            double E = 0.5 * (traj[i].v * traj[i].v + OMEGA*OMEGA * traj[i].x * traj[i].x);
            glVertex2f(i * DT, E);
        }
        glEnd();
    };
    
    draw_energy_series(traj_euler, 1.0, 0.2, 0.2);
    draw_energy_series(traj_implicit, 0.2, 0.2, 1.0);
    draw_energy_series(traj_verlet, 0.2, 1.0, 0.2);
    
    draw_text(0.5, 1.05, "Energy E(t)");
}

void draw_error() {
    glViewport(400, 400, 400, 400);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, N_STEPS * DT, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glColor3f(0.2, 0.2, 0.2);
    glBegin(GL_LINES);
    for (double t = 0; t <= N_STEPS * DT; t += T/4) {
        glVertex2f(t, -0.8); glVertex2f(t, 0.8);
    }
    for (double e = -0.8; e <= 0.8; e += 0.4) {
        glVertex2f(0, e); glVertex2f(N_STEPS * DT, e);
    }
    glEnd();
    auto draw_error_series = [&](const std::vector<State>& traj, float r, float g, float b) {
        glColor3f(r, g, b);
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i <= current_frame && i < (int)traj.size(); ++i) {
            double err = traj[i].x - analytical_x(i * DT);
            glVertex2f(i * DT, err);
        }
        glEnd();
    };
    draw_error_series(traj_euler, 1.0, 0.2, 0.2);
    draw_error_series(traj_implicit, 0.2, 0.2, 1.0);
    draw_error_series(traj_verlet, 0.2, 1.0, 0.2);
    draw_text(0.5, 0.9, "Error x(t) - analytical");
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    draw_phase_portrait();
    draw_position_vs_time();
    draw_energy();
    draw_error();
    glutSwapBuffers();
}

void timer(int value) {
    if (is_animating) {
        current_frame++;
        if (current_frame > N_STEPS) {
            current_frame = 0;
        }
        glutPostRedisplay();
    }
    glutTimerFunc(50, timer, 0);
}

void keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case ' ':
            is_animating = !is_animating;
            break;
        case 'r':
            current_frame = 0;
            glutPostRedisplay();
            break;
        case 27:
            exit(0);
            break;
    }
}

void init() {
    glClearColor(0.1, 0.1, 0.1, 1.0);
    traj_euler = solve_ode(euler_explicit);
    traj_implicit = solve_ode(euler_implicit);
    traj_verlet = solve_ode(verlet);
    printf("Данные гармонического осцилятора\n");
    printf("Угловая чистота = %.2f, Период = %.2f, Шаг = %.2f\n", OMEGA, T, DT);
    printf("Всего шагов: %d\n", N_STEPS);
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, w, 0, h);
    glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 800);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Sus");
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutTimerFunc(0, timer, 0);
    glutMainLoop();
    return 0;
}
