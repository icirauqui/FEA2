#include <iostream>
#include <vector>

#include <eigen3/Eigen/Dense>

class POS {

public:

  POS();

  void Transform(Eigen::Vector4d r_im, 
                 Eigen::Vector4d r_pt, 
                 Eigen::Vector3d t, 
                 double s);

  // Accessors
  // Return latest by default
  // User can request a specific set: 0, 1, 2, ...
  std::vector<Eigen::Vector3d> GetPoints(int idx = -1);
  std::pair<Eigen::Vector3d, Eigen::Vector4d> GetPose(int idx = -1);
  int LenHistory();

private:

  Eigen::Vector4d ConcatenateQuaternions(const Eigen::Vector4d& qvec1,
                                         const Eigen::Vector4d& qvec2);

  Eigen::Vector4d NormalizeQuaternion(const Eigen::Vector4d& qvec);

  Eigen::Vector3d QuaternionRotatePoint(const Eigen::Vector4d& qvec,
                                        const Eigen::Vector3d& point);

  std::vector<std::vector<Eigen::Vector3d>> points_;
  std::vector<std::pair<Eigen::Vector3d, Eigen::Vector4d>> pose_;

};

