/*
 * COPYRIGHT AND PERMISSION NOTICE
 * Penn Software MSCKF_VIO
 * Copyright (C) 2017 The Trustees of the University of Pennsylvania
 * All rights reserved.
 */

#ifndef MSCKF_VIO_MATH_UTILS_HPP
#define MSCKF_VIO_MATH_UTILS_HPP

#include <cmath>
#include <Eigen/Dense>

namespace msckf_vio {

/*
 *  @brief Create a skew-symmetric matrix from a 3-element vector.
 *  @note Performs the operation:
 *  w   ->  [  0 -w3  w2]
 *          [ w3   0 -w1]
 *          [-w2  w1   0]
 */ 
inline Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d& w) {
  Eigen::Matrix3d w_hat;
  w_hat(0, 0) = 0;
  w_hat(0, 1) = -w(2);
  w_hat(0, 2) = w(1);
  w_hat(1, 0) = w(2);
  w_hat(1, 1) = 0;
  w_hat(1, 2) = -w(0);
  w_hat(2, 0) = -w(1); 
  w_hat(2, 1) = w(0);
  w_hat(2, 2) = 0;
  return w_hat;
}

/*
 * @brief Normalize the given quaternion to unit quaternion.
 */
inline void quaternionNormalize(Eigen::Vector4d& q) {
  double norm = q.norm();
  q = q / norm;
  return;
}

/*
 * @brief Perform q1 * q2
 */
inline Eigen::Vector4d quaternionMultiplication(
    const Eigen::Vector4d& q1,
    const Eigen::Vector4d& q2) {
  Eigen::Matrix4d L;
  L(0, 0) =  q1(3); L(0, 1) =  q1(2); L(0, 2) = -q1(1); L(0, 3) =  q1(0);
  L(1, 0) = -q1(2); L(1, 1) =  q1(3); L(1, 2) =  q1(0); L(1, 3) =  q1(1);
  L(2, 0) =  q1(1); L(2, 1) = -q1(0); L(2, 2) =  q1(3); L(2, 3) =  q1(2);
  L(3, 0) = -q1(0); L(3, 1) = -q1(1); L(3, 2) = -q1(2); L(3, 3) =  q1(3);

  Eigen::Vector4d q = L * q2;
  quaternionNormalize(q);
  return q;
}

/*
 * @brief Convert the vector part of a quaternion to a
 *    full quaternion.
 * @note This function is useful to convert delta quaternion
 *    which is usually a 3x1 vector to a full quaternion.
 *    For more details, check Equation (238) and (239) in
 *    "Indirect Kalman Filter for 3D Attitude Estimation:
 *    A Tutorial for quaternion Algebra".
 */
inline Eigen::Vector4d smallAngleQuaternion(
    const Eigen::Vector3d& dtheta) {

  Eigen::Vector3d dq = dtheta / 2.0;
  Eigen::Vector4d q;
  double dq_square_norm = dq.squaredNorm();

  if (dq_square_norm <= 1) {
    q.head<3>() = dq;
    q(3) = std::sqrt(1-dq_square_norm);
  } else {
    q.head<3>() = dq;
    q(3) = 1;
    q = q / std::sqrt(1+dq_square_norm);
  }

  return q;
}

/*
 * @brief Convert a quaternion to the corresponding rotation matrix
 * @note Pay attention to the convention used. The function follows the
 *    conversion in "Indirect Kalman Filter for 3D Attitude Estimation:
 *    A Tutorial for Quaternion Algebra", Equation (78).
 *
 *    The input quaternion should be in the form
 *      [q1, q2, q3, q4(scalar)]^T
 */
inline Eigen::Matrix3d quaternionToRotation(
    const Eigen::Vector4d& q) {
  const Eigen::Vector3d& q_vec = q.block(0, 0, 3, 1);
  const double& q4 = q(3);
  Eigen::Matrix3d R =
    (2*q4*q4-1)*Eigen::Matrix3d::Identity() -
    2*q4*skewSymmetric(q_vec) +
    2*q_vec*q_vec.transpose();
  //TODO: Is it necessary to use the approximation equation
  //    (Equation (87)) when the rotation angle is small?
  return R;
}

/*
 * @brief Convert a rotation matrix to a quaternion.
 * @note Pay attention to the convention used. The function follows the
 *    conversion in "Indirect Kalman Filter for 3D Attitude Estimation:
 *    A Tutorial for Quaternion Algebra", Equation (78).
 *
 *    The input quaternion should be in the form
 *      [q1, q2, q3, q4(scalar)]^T
 */
inline Eigen::Vector4d rotationToQuaternion(
    const Eigen::Matrix3d& R) {
  Eigen::Vector4d score;
  score(0) = R(0, 0);
  score(1) = R(1, 1);
  score(2) = R(2, 2);
  score(3) = R.trace();

  int max_row = 0, max_col = 0;
  score.maxCoeff(&max_row, &max_col);

  Eigen::Vector4d q = Eigen::Vector4d::Zero();
  if (max_row == 0) {
    q(0) = std::sqrt(1+2*R(0, 0)-R.trace()) / 2.0;
    q(1) = (R(0, 1)+R(1, 0)) / (4*q(0));
    q(2) = (R(0, 2)+R(2, 0)) / (4*q(0));
    q(3) = (R(1, 2)-R(2, 1)) / (4*q(0));
  } else if (max_row == 1) {
    q(1) = std::sqrt(1+2*R(1, 1)-R.trace()) / 2.0;
    q(0) = (R(0, 1)+R(1, 0)) / (4*q(1));
    q(2) = (R(1, 2)+R(2, 1)) / (4*q(1));
    q(3) = (R(2, 0)-R(0, 2)) / (4*q(1));
  } else if (max_row == 2) {
    q(2) = std::sqrt(1+2*R(2, 2)-R.trace()) / 2.0;
    q(0) = (R(0, 2)+R(2, 0)) / (4*q(2));
    q(1) = (R(1, 2)+R(2, 1)) / (4*q(2));
    q(3) = (R(0, 1)-R(1, 0)) / (4*q(2));
  } else {
    q(3) = std::sqrt(1+R.trace()) / 2.0;
    q(0) = (R(1, 2)-R(2, 1)) / (4*q(3));
    q(1) = (R(2, 0)-R(0, 2)) / (4*q(3));
    q(2) = (R(0, 1)-R(1, 0)) / (4*q(3));
  }

  if (q(3) < 0) q = -q;
  quaternionNormalize(q);
  return q;
}

/*
/*
 * @brief Compute the SE_(3) exponential
 */
inline Eigen::VectorXd expSE_3(const Eigen::Vector3d xi_omega,
             const Eigen::VectorXd xi_x) {
    const double theta = std::sqrt(xi_omega.squaredNorm());
    Eigen::Matrix3d Omega_1 = skewSymmetric(xi_omega);
    Eigen::Matrix3d Omega_2 = Omega_1*Omega_1;
    if (theta < 0.000001) {
      return xi_x + 1/2*Omega_1*xi_x + 1/(2*3)*Omega_2*xi_x;
    } else {
      const double A = std::sin(theta)/theta;
      const double B = (1-std::cos(theta))/(theta*theta);
      const double C = (1-A)/(theta*theta);
      const Eigen::Matrix3d V = Eigen::Matrix3d::Identity() + B*Omega_1 + C*Omega_2;
    return V*xi_x;
    }
}

inline Eigen::Matrix3d expSO_3(const Eigen::Vector3d xi_omega) {
    const double theta = std::sqrt(xi_omega.squaredNorm());
    Eigen::Matrix3d Omega_1 = skewSymmetric(xi_omega);
    Eigen::Matrix3d Omega_2 = Omega_1*Omega_1;
    if (theta < 0.000001) {
      return Eigen::Matrix3d::Identity() + 1/2*Omega_1 + 1/(2*3)*Omega_2;
    } else {
      const double A = std::sin(theta)/theta;
      const double B = (1-std::cos(theta))/(theta*theta);
      const Eigen::Matrix3d dR = Eigen::Matrix3d::Identity() + A*Omega_1 + B*Omega_2;
    return dR;
    }
}

/*
@Additional formulas for Stable Embedding
*/
// Omega(omega) in quaternion dynamics
inline Eigen::Matrix4d SE_Omega(
    const Eigen::Vector3d& w) {
  Eigen::Matrix4d R = Eigen::Matrix4d::Zero();
  R.block<3, 3>(0, 0) = -skewSymmetric(w);
  R.block<3, 1>(0, 3) = w;
  R.block<1, 3>(3, 0) = -w;

  return R;
}

// Derivative of [(norm(q)^2 - 1)q] w.r.t(with respec to) q
inline Eigen::Matrix4d SE_A(
    const Eigen::Vector4d& q) {
  const double& q1 = q(0);
  const double& q2 = q(1);
  const double& q3 = q(2);
  const double& q4 = q(3);
  double q_square = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

  Eigen::Matrix4d R = Eigen::Matrix4d::Zero();
  R(0, 0) = q_square + 2 * q1 * q1 - 1;
  R(0, 1) = 2 * q1 * q2;
  R(0, 2) = 2 * q1 * q3;
  R(0, 3) = 2 * q1 * q4;
  R(1, 0) = 2 * q1 * q2;
  R(1, 1) = q_square + 2 * q2 * q2 - 1;
  R(1, 2) = 2 * q2 * q3;
  R(1, 3) = 2 * q2 * q4;
  R(2, 0) = 2 * q1 * q3;
  R(2, 1) = 2 * q2 * q3;
  R(2, 2) = q_square + 2 * q3 * q3 - 1;
  R(2, 3) = 2 * q3 * q4;
  R(3, 0) = 2 * q1 * q4;
  R(3, 1) = 2 * q2 * q4;
  R(3, 2) = 2 * q3 * q4;
  R(3, 3) = q_square + 2 * q4 * q4 - 1;

  return R;
}

// Derivative of [C(q)^T * a] w.r.t q
inline Eigen::MatrixXd SE_B(
    const Eigen::Vector4d& q, const Eigen::Vector3d& a) {
  const double& q1 = q(0);
  const double& q2 = q(1);
  const double& q3 = q(2);
  const double& q4 = q(3);
  const double& a1 = a(0);
  const double& a2 = a(1);
  const double& a3 = a(2);

  Eigen::Matrix<double, 3, 4> R = Eigen::Matrix<double, 3, 4>::Zero();
  R(0, 0) = 2 * (a1 * q1 + a2 * q2 + a3 * q3);
  R(0, 1) = 2 * (a2 * q1 - a1 * q2 + a3 * q4);
  R(0, 2) = 2 * (a3 * q1 - a1 * q3 - a2 * q4);
  R(0, 3) = 2 * (a1 * q4 - a2 * q3 + a3 * q2);
  R(1, 0) = 2 * (a1 * q2 - a2 * q1 - a3 * q4);
  R(1, 1) = 2 * (a1 * q1 + a2 * q2 + a3 * q3);
  R(1, 2) = 2 * (a1 * q4 - a2 * q3 + a3 * q2);
  R(1, 3) = 2 * (a1 * q3 - a3 * q1 + a2 * q4);
  R(2, 0) = 2 * (a1 * q3 - a3 * q1 + a2 * q4);
  R(2, 1) = 2 * (a2 * q3 - a1 * q4 - a3 * q2);
  R(2, 2) = 2 * (a1 * q1 + a2 * q2 + a3 * q3);
  R(2, 3) = 2 * (a2 * q1 - a1 * q2 + a3 * q4);

  return R;
}

inline Eigen::MatrixXd SE_Xi(
    const Eigen::Vector4d& q) {
  const double& q1 = q(0);
  const double& q2 = q(1);
  const double& q3 = q(2);
  const double& q4 = q(3);

  Eigen::Matrix<double, 4, 3> R = Eigen::Matrix<double, 4, 3>::Zero();
  R(0, 0) = q4;
  R(0, 1) = -q3;
  R(0, 2) = q2;
  R(1, 0) = q3;
  R(1, 1) = q4;
  R(1, 2) = -q1;
  R(2, 0) = -q2;
  R(2, 1) = q1;
  R(2, 2) = q4;
  R(3, 0) = -q1;
  R(3, 1) = -q2;
  R(3, 2) = -q3;

  return R;
}

inline Eigen::Matrix4d SE_M_L(
    const Eigen::Vector4d& q) {
  const double& q1 = q(0);
  const double& q2 = q(1);
  const double& q3 = q(2);
  const double& q4 = q(3);

  Eigen::Matrix4d R = Eigen::Matrix4d::Zero();
  R(0, 0) = q4;
  R(0, 1) = q3;
  R(0, 2) = -q2;
  R(0, 3) = q1;
  R(1, 0) = -q3;
  R(1, 1) = q4;
  R(1, 2) = q1;
  R(1, 3) = q2;
  R(2, 0) = q2;
  R(2, 1) = -q1;
  R(2, 2) = q4;
  R(2, 3) = q3;
  R(3, 0) = -q1;
  R(3, 1) = -q2;
  R(3, 2) = -q3;
  R(3, 3) = q4;

  return R;
}

inline Eigen::Matrix4d SE_M_R(
    const Eigen::Vector4d& q) {
  const double& q1 = q(0);
  const double& q2 = q(1);
  const double& q3 = q(2);
  const double& q4 = q(3);

  Eigen::Matrix4d R = Eigen::Matrix4d::Zero();
  R(0, 0) = q4;
  R(0, 1) = -q3;
  R(0, 2) = q2;
  R(0, 3) = q1;
  R(1, 0) = q3;
  R(1, 1) = q4;
  R(1, 2) = -q1;
  R(1, 3) = q2;
  R(2, 0) = -q2;
  R(2, 1) = q1;
  R(2, 2) = q4;
  R(2, 3) = q3;
  R(3, 0) = -q1;
  R(3, 1) = -q2;
  R(3, 2) = -q3;
  R(3, 3) = q4;

  return R;
}

inline Eigen::MatrixXd SE_E(
    const Eigen::Vector4d& q, const Eigen::Vector3d& a) {
  const double& q1 = q(0);
  const double& q2 = q(1);
  const double& q3 = q(2);
  const double& q4 = q(3);
  const double& a1 = a(0);
  const double& a2 = a(1);
  const double& a3 = a(2);

  Eigen::Matrix<double, 3, 4> R = Eigen::Matrix<double, 3, 4>::Zero();
  R(0, 0) = 2 * (a1 * q1 + a2 * q2 + a3 * q3);
  R(0, 1) = 2 * (a2 * q1 - a1 * q2 - a3 * q4);
  R(0, 2) = 2 * (a3 * q1 - a1 * q3 + a2 * q4);
  R(0, 3) = 2 * (a1 * q4 + a2 * q3 - a3 * q2);
  R(1, 0) = 2 * (a1 * q2 - a2 * q1 + a3 * q4);
  R(1, 1) = 2 * (a1 * q1 + a2 * q2 + a3 * q3);
  R(1, 2) = 2 * (a3 * q2 - a1 * q4 - a2 * q3);
  R(1, 3) = 2 * (a3 * q1 - a1 * q3 + a2 * q4);
  R(2, 0) = 2 * (a1 * q3 - a3 * q1 - a2 * q4);
  R(2, 1) = 2 * (a1 * q4 + a2 * q3 - a3 * q2);
  R(2, 2) = 2 * (a1 * q1 + a2 * q2 + a3 * q3);
  R(2, 3) = 2 * (a1 * q2 - a2 * q1 + a3 * q4);

  return R;
}

inline Eigen::Matrix3d SE_quaternionToRotation(
    const Eigen::Vector4d& q) {
  const double& q1 = q(0);
  const double& q2 = q(1);
  const double& q3 = q(2);
  const double& q4 = q(3);
  Eigen::Matrix<double, 3, 3> R = Eigen::Matrix<double, 3, 3>::Zero();
  R(0, 0) = q1 * q1 - q2 * q2 - q3 * q3 + q4 * q4;
  R(0, 1) = 2 * (q1 * q2 + q3 * q4);
  R(0, 2) = 2 * (q1 * q3 - q2 * q4);
  R(1, 0) = 2 * (q1 * q2 - q3 * q4);
  R(1, 1) = -q1 * q1 + q2 * q2 - q3 * q3 + q4 * q4;
  R(1, 2) = 2 * (q2 * q3 + q1 * q4);
  R(2, 0) = 2 * (q1 * q3 + q2 * q4);
  R(2, 1) = 2 * (q2 * q3 - q1 * q4);
  R(2, 2) = -q1 * q1 - q2 * q2 + q3 * q3 + q4 * q4;
  return R;
}

inline Eigen::Vector4d SE_quaternionMultiplication(
    const Eigen::Vector4d& q1,
    const Eigen::Vector4d& q2) {
  Eigen::Matrix4d L;
  L(0, 0) =  q1(3); L(0, 1) =  q1(2); L(0, 2) = -q1(1); L(0, 3) =  q1(0);
  L(1, 0) = -q1(2); L(1, 1) =  q1(3); L(1, 2) =  q1(0); L(1, 3) =  q1(1);
  L(2, 0) =  q1(1); L(2, 1) = -q1(0); L(2, 2) =  q1(3); L(2, 3) =  q1(2);
  L(3, 0) = -q1(0); L(3, 1) = -q1(1); L(3, 2) = -q1(2); L(3, 3) =  q1(3);

  Eigen::Vector4d q = L * q2;
  return q;
}

inline Eigen::Vector4d SE_quaternionNormalize(Eigen::Vector4d q) {
  double norm = q.norm();
  return q / norm;
}

inline double SE_quaternionsquared(Eigen::Vector4d q) {
  return q(0)*q(0) + q(1)*q(1) + q(2)*q(2) + q(3)*q(3);
}

} // end namespace msckf_vio

#endif // MSCKF_VIO_MATH_UTILS_HPP
