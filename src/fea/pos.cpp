#include "pos.hpp"


POS::POS() {}

void POS::Transform(Eigen::Matrix4d r, Eigen::Vector3d t, double s) {}

std::vector<Eigen::Vector3d> POS::GetPoints(int idx) {
  if (idx == -1) {
    idx = points_.size() - 1;
  }
  return points_[idx];
}

Eigen::Vector4d POS::GetPose(int idx) {
  if (idx == -1) {
    idx = pose_.size() - 1;
  }
  return pose_[idx];
}

int POS::LenHistory() {
  return points_.size();
}

void POS::Rotate() {}

void POS::Translate() {}

void POS::Scale() {}

