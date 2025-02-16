#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_8390988146680213198);
void live_err_fun(double *nom_x, double *delta_x, double *out_6244149898328306099);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_1055248739154000437);
void live_H_mod_fun(double *state, double *out_2477047662753292834);
void live_f_fun(double *state, double dt, double *out_4364589270647015481);
void live_F_fun(double *state, double dt, double *out_6681937296940613540);
void live_h_4(double *state, double *unused, double *out_4264504845329409309);
void live_H_4(double *state, double *unused, double *out_6518477235901727631);
void live_h_9(double *state, double *unused, double *out_56796046794078120);
void live_H_9(double *state, double *unused, double *out_4641047902543376515);
void live_h_10(double *state, double *unused, double *out_1663102567046158310);
void live_H_10(double *state, double *unused, double *out_1384400565122620472);
void live_h_12(double *state, double *unused, double *out_1618518243035899651);
void live_H_12(double *state, double *unused, double *out_6908810429775862190);
void live_h_35(double *state, double *unused, double *out_7884199092030081406);
void live_H_35(double *state, double *unused, double *out_8561604780435216609);
void live_h_32(double *state, double *unused, double *out_7656551497563665165);
void live_H_32(double *state, double *unused, double *out_4475869092338293509);
void live_h_13(double *state, double *unused, double *out_2857072050657840222);
void live_H_13(double *state, double *unused, double *out_2320253428634298347);
void live_h_14(double *state, double *unused, double *out_56796046794078120);
void live_H_14(double *state, double *unused, double *out_4641047902543376515);
void live_h_33(double *state, double *unused, double *out_202772564568973414);
void live_H_33(double *state, double *unused, double *out_5411047775796359005);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}