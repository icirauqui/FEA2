#include <iostream>
#include <vector>

#include <eigen3/Eigen/Dense>

class POS {

public:

  POS();

  void Transform(Eigen::Matrix4d r, Eigen::Vector3d t, double s);

  // Accessors
  // Return latest by default
  // User can request a specific set: 0, 1, 2, ...
  std::vector<Eigen::Vector3d> GetPoints(int idx = -1);
  Eigen::Vector4d GetPose(int idx = -1);
  int LenHistory();

private:

  void Rotate();

  void Translate();

  void Scale();

  std::vector<std::vector<Eigen::Vector3d>> points_;
  std::vector<Eigen::Matrix4d> pose_;

};