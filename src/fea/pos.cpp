#include "pos.hpp"


POS::POS() {}

void POS::Transform(Eigen::Vector4d r_im, 
                    Eigen::Vector4d r_pt, 
                    Eigen::Vector3d t, 
                    double s) {
  // Pose: Translate, Rotate, Scale
  std::pair<Eigen::Vector3d, Eigen::Vector4d> pose = GetPose();
  pose.first += t;
  pose.second = ConcatenateQuaternions(pose.second, r_im);
  pose.first *= s;
  pose_.push_back(pose);

  // Points: Translate, Rotate, Scale
  std::vector<Eigen::Vector3d> points = GetPoints();
  for (int i = 0; i < points.size(); i++) {
    points[i] += t;
    points[i] = QuaternionRotatePoint(r_pt, points[i]);
    points[i] *= s;
  }
  points_.push_back(points);

}

std::vector<Eigen::Vector3d> POS::GetPoints(int idx) {
  if (idx == -1) {
    idx = points_.size() - 1;
  }
  return points_[idx];
}

std::pair<Eigen::Vector3d, Eigen::Vector4d> POS::GetPose(int idx) {
  if (idx == -1) {
    idx = pose_.size() - 1;
  }
  return pose_[idx];
}

int POS::LenHistory() {
  return points_.size();
}

Eigen::Vector4d POS::ConcatenateQuaternions(const Eigen::Vector4d& qvec1,
                                            const Eigen::Vector4d& qvec2) {
  const Eigen::Vector4d normalized_qvec1 = NormalizeQuaternion(qvec1);
  const Eigen::Vector4d normalized_qvec2 = NormalizeQuaternion(qvec2);
  const Eigen::Quaterniond quat1(normalized_qvec1(0), normalized_qvec1(1),
                                 normalized_qvec1(2), normalized_qvec1(3));
  const Eigen::Quaterniond quat2(normalized_qvec2(0), normalized_qvec2(1),
                                 normalized_qvec2(2), normalized_qvec2(3));
  const Eigen::Quaterniond cat_quat = quat2 * quat1;
  return NormalizeQuaternion(
      Eigen::Vector4d(cat_quat.w(), cat_quat.x(), cat_quat.y(), cat_quat.z()));

}

Eigen::Vector4d POS::NormalizeQuaternion(const Eigen::Vector4d& qvec) {
  const double norm = qvec.norm();
  if (norm == 0) {
    return Eigen::Vector4d(1.0, qvec(1), qvec(2), qvec(3));
  } else {
    return qvec / norm;
  }
}

Eigen::Vector3d POS::QuaternionRotatePoint(const Eigen::Vector4d& qvec,
                                           const Eigen::Vector3d& point) {
  const Eigen::Vector4d normalized_qvec = NormalizeQuaternion(qvec);
  const Eigen::Quaterniond quat(normalized_qvec(0), normalized_qvec(1),
                                normalized_qvec(2), normalized_qvec(3));
  return quat * point;
}