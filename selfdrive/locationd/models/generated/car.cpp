#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8013343019937853643) {
   out_8013343019937853643[0] = delta_x[0] + nom_x[0];
   out_8013343019937853643[1] = delta_x[1] + nom_x[1];
   out_8013343019937853643[2] = delta_x[2] + nom_x[2];
   out_8013343019937853643[3] = delta_x[3] + nom_x[3];
   out_8013343019937853643[4] = delta_x[4] + nom_x[4];
   out_8013343019937853643[5] = delta_x[5] + nom_x[5];
   out_8013343019937853643[6] = delta_x[6] + nom_x[6];
   out_8013343019937853643[7] = delta_x[7] + nom_x[7];
   out_8013343019937853643[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_9096302114247577957) {
   out_9096302114247577957[0] = -nom_x[0] + true_x[0];
   out_9096302114247577957[1] = -nom_x[1] + true_x[1];
   out_9096302114247577957[2] = -nom_x[2] + true_x[2];
   out_9096302114247577957[3] = -nom_x[3] + true_x[3];
   out_9096302114247577957[4] = -nom_x[4] + true_x[4];
   out_9096302114247577957[5] = -nom_x[5] + true_x[5];
   out_9096302114247577957[6] = -nom_x[6] + true_x[6];
   out_9096302114247577957[7] = -nom_x[7] + true_x[7];
   out_9096302114247577957[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_86536487160361699) {
   out_86536487160361699[0] = 1.0;
   out_86536487160361699[1] = 0;
   out_86536487160361699[2] = 0;
   out_86536487160361699[3] = 0;
   out_86536487160361699[4] = 0;
   out_86536487160361699[5] = 0;
   out_86536487160361699[6] = 0;
   out_86536487160361699[7] = 0;
   out_86536487160361699[8] = 0;
   out_86536487160361699[9] = 0;
   out_86536487160361699[10] = 1.0;
   out_86536487160361699[11] = 0;
   out_86536487160361699[12] = 0;
   out_86536487160361699[13] = 0;
   out_86536487160361699[14] = 0;
   out_86536487160361699[15] = 0;
   out_86536487160361699[16] = 0;
   out_86536487160361699[17] = 0;
   out_86536487160361699[18] = 0;
   out_86536487160361699[19] = 0;
   out_86536487160361699[20] = 1.0;
   out_86536487160361699[21] = 0;
   out_86536487160361699[22] = 0;
   out_86536487160361699[23] = 0;
   out_86536487160361699[24] = 0;
   out_86536487160361699[25] = 0;
   out_86536487160361699[26] = 0;
   out_86536487160361699[27] = 0;
   out_86536487160361699[28] = 0;
   out_86536487160361699[29] = 0;
   out_86536487160361699[30] = 1.0;
   out_86536487160361699[31] = 0;
   out_86536487160361699[32] = 0;
   out_86536487160361699[33] = 0;
   out_86536487160361699[34] = 0;
   out_86536487160361699[35] = 0;
   out_86536487160361699[36] = 0;
   out_86536487160361699[37] = 0;
   out_86536487160361699[38] = 0;
   out_86536487160361699[39] = 0;
   out_86536487160361699[40] = 1.0;
   out_86536487160361699[41] = 0;
   out_86536487160361699[42] = 0;
   out_86536487160361699[43] = 0;
   out_86536487160361699[44] = 0;
   out_86536487160361699[45] = 0;
   out_86536487160361699[46] = 0;
   out_86536487160361699[47] = 0;
   out_86536487160361699[48] = 0;
   out_86536487160361699[49] = 0;
   out_86536487160361699[50] = 1.0;
   out_86536487160361699[51] = 0;
   out_86536487160361699[52] = 0;
   out_86536487160361699[53] = 0;
   out_86536487160361699[54] = 0;
   out_86536487160361699[55] = 0;
   out_86536487160361699[56] = 0;
   out_86536487160361699[57] = 0;
   out_86536487160361699[58] = 0;
   out_86536487160361699[59] = 0;
   out_86536487160361699[60] = 1.0;
   out_86536487160361699[61] = 0;
   out_86536487160361699[62] = 0;
   out_86536487160361699[63] = 0;
   out_86536487160361699[64] = 0;
   out_86536487160361699[65] = 0;
   out_86536487160361699[66] = 0;
   out_86536487160361699[67] = 0;
   out_86536487160361699[68] = 0;
   out_86536487160361699[69] = 0;
   out_86536487160361699[70] = 1.0;
   out_86536487160361699[71] = 0;
   out_86536487160361699[72] = 0;
   out_86536487160361699[73] = 0;
   out_86536487160361699[74] = 0;
   out_86536487160361699[75] = 0;
   out_86536487160361699[76] = 0;
   out_86536487160361699[77] = 0;
   out_86536487160361699[78] = 0;
   out_86536487160361699[79] = 0;
   out_86536487160361699[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2055691188639617943) {
   out_2055691188639617943[0] = state[0];
   out_2055691188639617943[1] = state[1];
   out_2055691188639617943[2] = state[2];
   out_2055691188639617943[3] = state[3];
   out_2055691188639617943[4] = state[4];
   out_2055691188639617943[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2055691188639617943[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2055691188639617943[7] = state[7];
   out_2055691188639617943[8] = state[8];
}
void F_fun(double *state, double dt, double *out_3210748245983153222) {
   out_3210748245983153222[0] = 1;
   out_3210748245983153222[1] = 0;
   out_3210748245983153222[2] = 0;
   out_3210748245983153222[3] = 0;
   out_3210748245983153222[4] = 0;
   out_3210748245983153222[5] = 0;
   out_3210748245983153222[6] = 0;
   out_3210748245983153222[7] = 0;
   out_3210748245983153222[8] = 0;
   out_3210748245983153222[9] = 0;
   out_3210748245983153222[10] = 1;
   out_3210748245983153222[11] = 0;
   out_3210748245983153222[12] = 0;
   out_3210748245983153222[13] = 0;
   out_3210748245983153222[14] = 0;
   out_3210748245983153222[15] = 0;
   out_3210748245983153222[16] = 0;
   out_3210748245983153222[17] = 0;
   out_3210748245983153222[18] = 0;
   out_3210748245983153222[19] = 0;
   out_3210748245983153222[20] = 1;
   out_3210748245983153222[21] = 0;
   out_3210748245983153222[22] = 0;
   out_3210748245983153222[23] = 0;
   out_3210748245983153222[24] = 0;
   out_3210748245983153222[25] = 0;
   out_3210748245983153222[26] = 0;
   out_3210748245983153222[27] = 0;
   out_3210748245983153222[28] = 0;
   out_3210748245983153222[29] = 0;
   out_3210748245983153222[30] = 1;
   out_3210748245983153222[31] = 0;
   out_3210748245983153222[32] = 0;
   out_3210748245983153222[33] = 0;
   out_3210748245983153222[34] = 0;
   out_3210748245983153222[35] = 0;
   out_3210748245983153222[36] = 0;
   out_3210748245983153222[37] = 0;
   out_3210748245983153222[38] = 0;
   out_3210748245983153222[39] = 0;
   out_3210748245983153222[40] = 1;
   out_3210748245983153222[41] = 0;
   out_3210748245983153222[42] = 0;
   out_3210748245983153222[43] = 0;
   out_3210748245983153222[44] = 0;
   out_3210748245983153222[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3210748245983153222[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3210748245983153222[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3210748245983153222[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3210748245983153222[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3210748245983153222[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3210748245983153222[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3210748245983153222[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3210748245983153222[53] = -9.8000000000000007*dt;
   out_3210748245983153222[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3210748245983153222[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3210748245983153222[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3210748245983153222[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3210748245983153222[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3210748245983153222[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3210748245983153222[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3210748245983153222[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3210748245983153222[62] = 0;
   out_3210748245983153222[63] = 0;
   out_3210748245983153222[64] = 0;
   out_3210748245983153222[65] = 0;
   out_3210748245983153222[66] = 0;
   out_3210748245983153222[67] = 0;
   out_3210748245983153222[68] = 0;
   out_3210748245983153222[69] = 0;
   out_3210748245983153222[70] = 1;
   out_3210748245983153222[71] = 0;
   out_3210748245983153222[72] = 0;
   out_3210748245983153222[73] = 0;
   out_3210748245983153222[74] = 0;
   out_3210748245983153222[75] = 0;
   out_3210748245983153222[76] = 0;
   out_3210748245983153222[77] = 0;
   out_3210748245983153222[78] = 0;
   out_3210748245983153222[79] = 0;
   out_3210748245983153222[80] = 1;
}
void h_25(double *state, double *unused, double *out_6302416998799214083) {
   out_6302416998799214083[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2935156384893243290) {
   out_2935156384893243290[0] = 0;
   out_2935156384893243290[1] = 0;
   out_2935156384893243290[2] = 0;
   out_2935156384893243290[3] = 0;
   out_2935156384893243290[4] = 0;
   out_2935156384893243290[5] = 0;
   out_2935156384893243290[6] = 1;
   out_2935156384893243290[7] = 0;
   out_2935156384893243290[8] = 0;
}
void h_24(double *state, double *unused, double *out_8374539201565141696) {
   out_8374539201565141696[0] = state[4];
   out_8374539201565141696[1] = state[5];
}
void H_24(double *state, double *unused, double *out_622769213909079846) {
   out_622769213909079846[0] = 0;
   out_622769213909079846[1] = 0;
   out_622769213909079846[2] = 0;
   out_622769213909079846[3] = 0;
   out_622769213909079846[4] = 1;
   out_622769213909079846[5] = 0;
   out_622769213909079846[6] = 0;
   out_622769213909079846[7] = 0;
   out_622769213909079846[8] = 0;
   out_622769213909079846[9] = 0;
   out_622769213909079846[10] = 0;
   out_622769213909079846[11] = 0;
   out_622769213909079846[12] = 0;
   out_622769213909079846[13] = 0;
   out_622769213909079846[14] = 1;
   out_622769213909079846[15] = 0;
   out_622769213909079846[16] = 0;
   out_622769213909079846[17] = 0;
}
void h_30(double *state, double *unused, double *out_4263054424156288348) {
   out_4263054424156288348[0] = state[4];
}
void H_30(double *state, double *unused, double *out_416823426385994663) {
   out_416823426385994663[0] = 0;
   out_416823426385994663[1] = 0;
   out_416823426385994663[2] = 0;
   out_416823426385994663[3] = 0;
   out_416823426385994663[4] = 1;
   out_416823426385994663[5] = 0;
   out_416823426385994663[6] = 0;
   out_416823426385994663[7] = 0;
   out_416823426385994663[8] = 0;
}
void h_26(double *state, double *unused, double *out_8268058895819400932) {
   out_8268058895819400932[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6676659703767299514) {
   out_6676659703767299514[0] = 0;
   out_6676659703767299514[1] = 0;
   out_6676659703767299514[2] = 0;
   out_6676659703767299514[3] = 0;
   out_6676659703767299514[4] = 0;
   out_6676659703767299514[5] = 0;
   out_6676659703767299514[6] = 0;
   out_6676659703767299514[7] = 1;
   out_6676659703767299514[8] = 0;
}
void h_27(double *state, double *unused, double *out_7142653807249491888) {
   out_7142653807249491888[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2591586738186419574) {
   out_2591586738186419574[0] = 0;
   out_2591586738186419574[1] = 0;
   out_2591586738186419574[2] = 0;
   out_2591586738186419574[3] = 1;
   out_2591586738186419574[4] = 0;
   out_2591586738186419574[5] = 0;
   out_2591586738186419574[6] = 0;
   out_2591586738186419574[7] = 0;
   out_2591586738186419574[8] = 0;
}
void h_29(double *state, double *unused, double *out_3540825009639387287) {
   out_3540825009639387287[0] = state[1];
}
void H_29(double *state, double *unused, double *out_93407917928397521) {
   out_93407917928397521[0] = 0;
   out_93407917928397521[1] = 1;
   out_93407917928397521[2] = 0;
   out_93407917928397521[3] = 0;
   out_93407917928397521[4] = 0;
   out_93407917928397521[5] = 0;
   out_93407917928397521[6] = 0;
   out_93407917928397521[7] = 0;
   out_93407917928397521[8] = 0;
}
void h_28(double *state, double *unused, double *out_5438491754500665626) {
   out_5438491754500665626[0] = state[0];
}
void H_28(double *state, double *unused, double *out_4988991099141133053) {
   out_4988991099141133053[0] = 1;
   out_4988991099141133053[1] = 0;
   out_4988991099141133053[2] = 0;
   out_4988991099141133053[3] = 0;
   out_4988991099141133053[4] = 0;
   out_4988991099141133053[5] = 0;
   out_4988991099141133053[6] = 0;
   out_4988991099141133053[7] = 0;
   out_4988991099141133053[8] = 0;
}
void h_31(double *state, double *unused, double *out_2664967470545191733) {
   out_2664967470545191733[0] = state[8];
}
void H_31(double *state, double *unused, double *out_7302867806000650990) {
   out_7302867806000650990[0] = 0;
   out_7302867806000650990[1] = 0;
   out_7302867806000650990[2] = 0;
   out_7302867806000650990[3] = 0;
   out_7302867806000650990[4] = 0;
   out_7302867806000650990[5] = 0;
   out_7302867806000650990[6] = 0;
   out_7302867806000650990[7] = 0;
   out_7302867806000650990[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_8013343019937853643) {
  err_fun(nom_x, delta_x, out_8013343019937853643);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_9096302114247577957) {
  inv_err_fun(nom_x, true_x, out_9096302114247577957);
}
void car_H_mod_fun(double *state, double *out_86536487160361699) {
  H_mod_fun(state, out_86536487160361699);
}
void car_f_fun(double *state, double dt, double *out_2055691188639617943) {
  f_fun(state,  dt, out_2055691188639617943);
}
void car_F_fun(double *state, double dt, double *out_3210748245983153222) {
  F_fun(state,  dt, out_3210748245983153222);
}
void car_h_25(double *state, double *unused, double *out_6302416998799214083) {
  h_25(state, unused, out_6302416998799214083);
}
void car_H_25(double *state, double *unused, double *out_2935156384893243290) {
  H_25(state, unused, out_2935156384893243290);
}
void car_h_24(double *state, double *unused, double *out_8374539201565141696) {
  h_24(state, unused, out_8374539201565141696);
}
void car_H_24(double *state, double *unused, double *out_622769213909079846) {
  H_24(state, unused, out_622769213909079846);
}
void car_h_30(double *state, double *unused, double *out_4263054424156288348) {
  h_30(state, unused, out_4263054424156288348);
}
void car_H_30(double *state, double *unused, double *out_416823426385994663) {
  H_30(state, unused, out_416823426385994663);
}
void car_h_26(double *state, double *unused, double *out_8268058895819400932) {
  h_26(state, unused, out_8268058895819400932);
}
void car_H_26(double *state, double *unused, double *out_6676659703767299514) {
  H_26(state, unused, out_6676659703767299514);
}
void car_h_27(double *state, double *unused, double *out_7142653807249491888) {
  h_27(state, unused, out_7142653807249491888);
}
void car_H_27(double *state, double *unused, double *out_2591586738186419574) {
  H_27(state, unused, out_2591586738186419574);
}
void car_h_29(double *state, double *unused, double *out_3540825009639387287) {
  h_29(state, unused, out_3540825009639387287);
}
void car_H_29(double *state, double *unused, double *out_93407917928397521) {
  H_29(state, unused, out_93407917928397521);
}
void car_h_28(double *state, double *unused, double *out_5438491754500665626) {
  h_28(state, unused, out_5438491754500665626);
}
void car_H_28(double *state, double *unused, double *out_4988991099141133053) {
  H_28(state, unused, out_4988991099141133053);
}
void car_h_31(double *state, double *unused, double *out_2664967470545191733) {
  h_31(state, unused, out_2664967470545191733);
}
void car_H_31(double *state, double *unused, double *out_7302867806000650990) {
  H_31(state, unused, out_7302867806000650990);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
