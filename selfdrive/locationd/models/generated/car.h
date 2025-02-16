#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_8013343019937853643);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_9096302114247577957);
void car_H_mod_fun(double *state, double *out_86536487160361699);
void car_f_fun(double *state, double dt, double *out_2055691188639617943);
void car_F_fun(double *state, double dt, double *out_3210748245983153222);
void car_h_25(double *state, double *unused, double *out_6302416998799214083);
void car_H_25(double *state, double *unused, double *out_2935156384893243290);
void car_h_24(double *state, double *unused, double *out_8374539201565141696);
void car_H_24(double *state, double *unused, double *out_622769213909079846);
void car_h_30(double *state, double *unused, double *out_4263054424156288348);
void car_H_30(double *state, double *unused, double *out_416823426385994663);
void car_h_26(double *state, double *unused, double *out_8268058895819400932);
void car_H_26(double *state, double *unused, double *out_6676659703767299514);
void car_h_27(double *state, double *unused, double *out_7142653807249491888);
void car_H_27(double *state, double *unused, double *out_2591586738186419574);
void car_h_29(double *state, double *unused, double *out_3540825009639387287);
void car_H_29(double *state, double *unused, double *out_93407917928397521);
void car_h_28(double *state, double *unused, double *out_5438491754500665626);
void car_H_28(double *state, double *unused, double *out_4988991099141133053);
void car_h_31(double *state, double *unused, double *out_2664967470545191733);
void car_H_31(double *state, double *unused, double *out_7302867806000650990);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}