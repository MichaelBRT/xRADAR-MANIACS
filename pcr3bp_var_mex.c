#include "mex.h"
#include "math.h"

// Helper: fill in the Jacobian matrix A(t)
void compute_jacobian(double mu, const double *x, double A[4][4]) {
    double x1 = x[0], x2 = x[1];
    double r1 = sqrt((x1 + mu)*(x1 + mu) + x2*x2);
    double r2 = sqrt((x1 - 1 + mu)*(x1 - 1 + mu) + x2*x2);

    double r1_3 = r1*r1*r1, r2_3 = r2*r2*r2;
    double r1_5 = r1_3 * r1 * r1;
    double r2_5 = r2_3 * r2 * r2;

    double dUxx = 1
        - (1 - mu)*(1.0 / r1_3 - 3*(x1 + mu)*(x1 + mu) / r1_5)
        - mu*(1.0 / r2_3 - 3*(x1 - 1 + mu)*(x1 - 1 + mu) / r2_5);

    double dUxy = 3*x2 * (
        (1 - mu)*(x1 + mu) / r1_5 +
        mu*(x1 - 1 + mu) / r2_5
    );

    double dUyy = 1
        - (1 - mu)*(1.0 / r1_3 - 3*x2*x2 / r1_5)
        - mu*(1.0 / r2_3 - 3*x2*x2 / r2_5);

    // Fill Jacobian A
    A[0][0] = 0;   A[0][1] = 0;   A[0][2] = 1;  A[0][3] = 0;
    A[1][0] = 0;   A[1][1] = 0;   A[1][2] = 0;  A[1][3] = 1;
    A[2][0] = dUxx; A[2][1] = dUxy; A[2][2] = 0;  A[2][3] = 2;
    A[3][0] = dUxy; A[3][1] = dUyy; A[3][2] = -2; A[3][3] = 0;
}

void pcr3bp_var(const double mu, const double *X, double *dX) {
    const double *Phi = X;        // First 16 entries: STM
    const double *x = X + 16;     // Last 4 entries: state [x, y, vx, vy]

    double *dPhi = dX;            // First 16 entries of output
    double *dx = dX + 16;         // Last 4 entries of output

    // === Step 1: Compute dynamics
    double x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3];
    double r1 = sqrt((x1 + mu)*(x1 + mu) + x2*x2);
    double r2 = sqrt((x1 - 1 + mu)*(x1 - 1 + mu) + x2*x2);
    double r1_3 = r1*r1*r1;
    double r2_3 = r2*r2*r2;

    dx[0] = x3;
    dx[1] = x4;
    dx[2] = x1 + 2*x4 - (1 - mu)*(x1 + mu)/r1_3 - mu*(x1 - 1 + mu)/r2_3;
    dx[3] = x2 - 2*x3 - (1 - mu)*x2/r1_3 - mu*x2/r2_3;

    // === Step 2: Jacobian and dPhi
    double A[4][4];
    compute_jacobian(mu, x, A);

    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 4; ++row) {
            double sum = 0.0;
            for (int k = 0; k < 4; ++k) {
                sum += A[row][k] * Phi[4*k + col];  // Phi is column-major
            }
            dPhi[4*col + row] = sum;
        }
    }
}

// MEX gateway
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double mu = mxGetScalar(prhs[0]);
    const double *X = mxGetPr(prhs[1]); // [Phi(:); x]

    plhs[0] = mxCreateDoubleMatrix(20, 1, mxREAL);
    double *dX = mxGetPr(plhs[0]);

    pcr3bp_var(mu, X, dX);
}