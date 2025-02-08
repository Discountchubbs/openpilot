from openpilot.common.numpy_fast import clip, interp
from cereal import car
from typing import Any
from openpilot.common.realtime import DT_CTRL
from openpilot.common.params import Params
from openpilot.selfdrive.controls.lib.longcontrol import LongControl
from openpilot.selfdrive.controls.lib.longitudinal_mpc_lib.long_mpc import LongitudinalMpc
from openpilot.selfdrive.car.hyundai.values import HyundaiFlags, CarControllerParams


LongCtrlState = car.CarControl.Actuators.LongControlState

class JerkOutput:
  def __init__(self, jerk_upper_limit, jerk_lower_limit, cb_upper, cb_lower):
    self.jerk_upper_limit = jerk_upper_limit
    self.jerk_lower_limit = jerk_lower_limit
    self.cb_upper = cb_upper
    self.cb_lower = cb_lower

class HKGLongitudinalTuning:
  def __init__(self, CP) -> None:
    self.CP = CP
    self._setup_controllers()
    self._init_state()

  def _setup_controllers(self) -> None:
    self.mpc = LongitudinalMpc(mode='acc')

    self.long_control = LongControl(self.CP)
    self.DT_CTRL = DT_CTRL
    self.params = Params()
    self.hkg_tuning = self.params.get_bool('HKGtuning')
    self.has_radar = self.params.get_bool("HyundaiRadarTracks")

  def _init_state(self) -> None:
    self.last_accel = 0.0
    self.brake_ramp = 0.0
    self.accel_last = 0.0
    self.using_e2e = False
    self.accel_raw = 0.0
    self.accel_last_jerk = 0.0
    self.jerk = 0.0
    self.jerk_count = 0.0
    self.jerk_upper_limit = 0.0
    self.jerk_lower_limit = 0.0
    self.cb_upper = self.cb_lower = 0.0


  def make_jerk(self, CS, accel, actuators):
    state = getattr(actuators, "longControlState", LongCtrlState.pid)
    if state == LongCtrlState.stopping:
      self.jerk = 0.5 / 2 - CS.out.aEgo
      return self.jerk
    else:
      current_accel = clip(actuators.accel, CarControllerParams.ACCEL_MIN, CarControllerParams.ACCEL_MAX)
      self.jerk = (current_accel - self.accel_last_jerk) / self.DT_CTRL
      self.accel_last_jerk = current_accel

    jerk_max = 5.0
    v_error = abs(CS.out.vEgo - CS.out.cruiseState.speed)

    if self.CP.flags & HyundaiFlags.CANFD.value:
      if v_error < 3.0:
        jerk_reduction = interp(v_error,
                              [0.0, 1.0, 3.0],
                              [0.3, 0.5, 1.0])
        jerk_max *= jerk_reduction

      self.jerk_upper_limit = min(max(0.5, self.jerk * 2.0), jerk_max)
      self.jerk_lower_limit = min(max(1.0, -self.jerk * 4.0), jerk_max)
      self.cb_upper = self.cb_lower = 0.0
    else:
      self.jerk_upper_limit = min(max(0.5, self.jerk * 2.0), jerk_max)
      self.jerk_lower_limit = min(max(1.0, -self.jerk * 2.0), jerk_max)
      self.cb_upper = clip(0.9 + accel * 0.2, 0, 1.2)
      self.cb_lower = clip(0.8 + accel * 0.2, 0, 1.2)

    return self.jerk

  def calculate_limited_accel(self, accel, actuators, CS, LongCtrlState):
    """Limit acceleration based on detailed jerk limits."""
    self.make_jerk(CS, accel, actuators)
    accel_delta = accel - self.accel_last
    brake_aggressiveness = 0.0
    ramp_rate = 0.7
    if accel < 0:
      brake_ratio = clip(abs(accel / CarControllerParams.ACCEL_MIN), 0.0, 1.0)
      brake_aggressiveness = brake_ratio ** 1.5
      ramp_rate = interp(brake_aggressiveness,
                         [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
                         [0.2, 0.3, 0.4, 0.6, 0.9, 1.2])
      if brake_ratio > 0.8:
            ramp_rate *= 0.8
      if self.accel_last >= 0:
        self.brake_ramp = 0.0
        ramp_rate *= 0.5
    self.brake_ramp = min(1.0, self.brake_ramp + (ramp_rate * self.DT_CTRL))
    abs_accel = clip(abs(accel), 0.0, 2.0)

    smooth_factor = interp(abs_accel,
                          [0.0, 0.3, 0.6, 1.0, 1.5, 2.0],
                          [0.95, 0.90, 0.85, 0.75, 0.65, 0.55])

    v_error = abs(CS.out.vEgo - CS.out.cruiseState.speed)
    if v_error < 3.0:  # Within 3 m/s of target
        smooth_factor *= interp(v_error,
                              [0.0, 1.5, 3.0],
                              [0.3, 0.5, 0.8])

    if accel < 0 and brake_aggressiveness > 0.8:
      smooth_factor *= 0.8
    if abs_accel < 0.5:
      smooth_factor *= interp(abs_accel, [0.0, 0.5], [1.3, 1.0])
    accel_delta *= (smooth_factor * self.brake_ramp)
    if abs(accel_delta) > 0.5:
      damping = interp(abs(accel_delta),
                       [0.5, 1.0, 1.5, 2.0],
                       [0.90, 0.80, 0.70, 0.60])
      accel_delta *= damping

    limited_accel_delta = clip(accel_delta,
                               -self.jerk_lower_limit * self.DT_CTRL,
                                self.jerk_upper_limit * self.DT_CTRL)

    if accel < 0 and self.accel_last > 0:
      transition_factor = interp(abs(accel),
                                 [0, 0.5, 1.0, 1.5, 2.0],
                                 [0.3, 0.4, 0.5, 0.6, 0.7])
      limited_accel_delta *= transition_factor
    accel = self.accel_last + limited_accel_delta
    self.accel_last = accel
    return accel

  def update(self, accel, actuators, CS) -> float:
    accel = self.calculate_limited_accel(accel, actuators, CS)
    return accel

  def apply_tune(self, CP: Any) -> None:
    CP.vEgoStopping = 0.2
    CP.vEgoStarting = 0.1
    CP.stoppingDecelRate = 0.01
    CP.startAccel = 1.0 if bool(CP.flags & (HyundaiFlags.HYBRID | HyundaiFlags.EV)) else 1.6
    CP.startingState = True

  def get_jerk(self) -> JerkOutput:
    return JerkOutput(
      self.jerk_upper_limit,
      self.jerk_lower_limit,
      self.cb_upper,
      self.cb_lower,
    )

class HKGLongitudinalController:
    def __init__(self, CP):
        if Params().get_bool("HKGtuning"):
            self.tuning = HKGLongitudinalTuning(CP)
        else:
            self.tuning = None

    def apply_tune(self, CP):
        if self.tuning:
            self.tuning.apply_tune(CP)

    def update(self, accel, CS):
        if self.tuning:
            return self.tuning.update(accel, CS)
        return accel

    def get_jerk(self) -> JerkOutput:
        if self.tuning:
            return self.tuning.get_jerk()
        return 0.0
