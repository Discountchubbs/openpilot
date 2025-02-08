from types import SimpleNamespace
from openpilot.common.params import Params

class ParamManager:
  def __init__(self):
    self._params_list: SimpleNamespace = self._create_namespace({
      "hyundai_radar_tracks_available": False,
      "hyundai_radar_tracks_available_cache": False,
    })

  @staticmethod
  def _create_namespace(data: dict) -> SimpleNamespace:
    return SimpleNamespace(**data)

  def get_params(self) -> SimpleNamespace:
    return self._params_list

  def update(self, params: Params) -> None:
    self._params_list = self._create_namespace({
      "hyundai_radar_tracks_available": params.get_bool("HyundaiRadarTracksAvailable"),
      "hyundai_radar_tracks_available_cache": params.get_bool("HyundaiRadarTracksAvailableCache"),
    })
