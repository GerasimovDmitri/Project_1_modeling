#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <GL/glut.h>

#ifdef _WIN32
#include <windows.h>
#endif

// Параметры системы
typedef struct {
    double delta, alpha, beta, gamma, omega;
    double t, x, v;
} DuffingSystem;

DuffingSystem sys = {
    .delta = 0.2,
    .alpha = -1.0,
    .beta  = 1.0,
    .gamma = 0.3,
    .omega = 1.0,
    .t = 0.0,
    .x = 1.0,
    .v = 0.0
};

// Параметры анимации
double dt = 0.01;
int running = 1;
double phase_scale = 2.5;
double time_scale = 2.5;

// Хранение истории для следа (хвоста)
#define TRAIL_LEN 200
double trail_x[TRAIL_LEN];
double trail_v[TRAIL_LEN];
int trail_idx = 0;

// Функция правой части ОДУ
void duffing_rhs(double t, double x, double v, double *dx, double *dv) {
    *dx = v;
    *dv = -sys.delta * v - sys.alpha * x - sys.beta * x*x*x 
          + sys.gamma * cos(sys.omega * t);
}

// Шаг Рунге-Кутты 4-го порядка
void rk4_step(void) {
    double k1x, k1v, k2x, k2v, k3x, k3v, k4x, k4v;
    double dx1, dv1, dx2, dv2, dx3, dv3, dx4, dv4;
    
    duffing_rhs(sys.t, sys.x, sys.v, &dx1, &dv1);
    k1x = dt * dx1;
    k1v = dt * dv1;
    
    duffing_rhs(sys.t + dt/2.0, sys.x + k1x/2.0, sys.v + k1v/2.0, &dx2, &dv2);
    k2x = dt * dx2;
    k2v = dt * dv2;
    
    duffing_rhs(sys.t + dt/2.0, sys.x + k2x/2.0, sys.v + k2v/2.0, &dx3, &dv3);
    k3x = dt * dx3;
    k3v = dt * dv3;
    
    duffing_rhs(sys.t + dt, sys.x + k3x, sys.v + k3v, &dx4, &dv4);
    k4x = dt * dx4;
    k4v = dt * dv4;
    
    sys.x = sys.x + (k1x + 2.0*k2x + 2.0*k3x + k4x) / 6.0;
    sys.v = sys.v + (k1v + 2.0*k2v + 2.0*k3v + k4v) / 6.0;
    sys.t = sys.t + dt;
}

// Функция отрисовки
void display(void) {
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    
    // Отрисовка осей
    glColor3f(0.5, 0.5, 0.5);
    glBegin(GL_LINES);
        glVertex2f(-phase_scale, 0.0);
        glVertex2f( phase_scale, 0.0);
        glVertex2f(0.0, -phase_scale);
        glVertex2f(0.0,  phase_scale);
    glEnd();
    
    // Отрисовка всей траектории
    glColor4f(0.5, 0.5, 0.5, 0.3);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < TRAIL_LEN; i++) {
        int idx = (trail_idx + i) % TRAIL_LEN;
        if (trail_x[idx] != 0 || trail_v[idx] != 0 || i == 0) {
            glVertex2f(trail_x[idx], trail_v[idx]);
        }
    }
    glEnd();
    
    // Отрисовка следа
    glColor3f(0.8, 0.2, 0.2);
    glLineWidth(2.0);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < TRAIL_LEN; i++) {
        int idx = (trail_idx + i) % TRAIL_LEN;
        if (trail_x[idx] != 0 || trail_v[idx] != 0 || i == 0) {
            glVertex2f(trail_x[idx], trail_v[idx]);
        }
    }
    glEnd();
    
    // Отрисовка текущей точки
    glPointSize(8.0);
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_POINTS);
        glVertex2f(sys.x, sys.v);
    glEnd();
    
    // Отрисовка текста с параметрами
    glColor3f(1.0, 1.0, 1.0);
    glRasterPos2f(-phase_scale + 0.2, phase_scale - 0.3);
    char info[256];
    sprintf(info, "t = %.2f   x = %.3f   v = %.3f", sys.t, sys.x, sys.v);
    for (char *c = info; *c; c++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *c);
    
    glRasterPos2f(-phase_scale + 0.2, phase_scale - 0.6);
    sprintf(info, "delta=%.1f  alpha=%.1f  beta=%.1f  gamma=%.1f  omega=%.1f",
            sys.delta, sys.alpha, sys.beta, sys.gamma, sys.omega);
    for (char *c = info; *c; c++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    
    glutSwapBuffers();
}

// Обновление состояния
void update(int value) {
    if (running) {
        trail_x[trail_idx] = sys.x;
        trail_v[trail_idx] = sys.v;
        trail_idx = (trail_idx + 1) % TRAIL_LEN;
        
        for (int i = 0; i < 5; i++) {
            rk4_step();
        }
    }
    
    glutPostRedisplay();
    glutTimerFunc(16, update, 0);  // ~60 FPS
}

// Обработка клавиш
void keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case ' ':
            running = !running;
            break;
        case 'r':
            sys.t = 0.0;
            sys.x = 1.0;
            sys.v = 0.0;
            for (int i = 0; i < TRAIL_LEN; i++) {
                trail_x[i] = 0.0;
                trail_v[i] = 0.0;
            }
            trail_idx = 0;
            break;
        case '+':
        case '=':
            dt *= 0.8;
            if (dt < 0.0001) dt = 0.0001;
            break;
        case '-':
            dt *= 1.2;
            if (dt > 0.1) dt = 0.1;
            break;
        case 'q':
        case 27:
            exit(0);
            break;
    }
    glutPostRedisplay();
}

// Инициализация OpenGL
void init_gl(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 800);
    glutCreateWindow("Осциллятор Дуффинга - Анимация фазового портрета");
    
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-phase_scale, phase_scale, -phase_scale, phase_scale);
    glMatrixMode(GL_MODELVIEW);
    
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutTimerFunc(16, update, 0);
}

int main(int argc, char **argv) {
    printf("Осциллятор Дуффинга - OpenGL анимация\n");
    printf("\n");
    printf("Управление:\n");
    printf("  Пробел - пауза/запуск\n");
    printf("  R      - сброс\n");
    printf("  +/-    - изменение шага dt\n");
    printf("  Q/ESC  - выход\n");
    printf("\n");
    
    init_gl(argc, argv);
    glutMainLoop();
    
    return 0;
}
