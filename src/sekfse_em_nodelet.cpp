/*
 * COPYRIGHT AND PERMISSION NOTICE
 * Penn Software MSCKF_VIO
 * Copyright (C) 2017 The Trustees of the University of Pennsylvania
 * All rights reserved.
 */

#include <msckf_vio/sekfse_em_nodelet.h>

namespace msckf_vio {
void SekfseEMNodelet::onInit() {
  sekfse_em_ptr.reset(new SekfseEM(getPrivateNodeHandle()));
  if (!sekfse_em_ptr->initialize()) {
    ROS_ERROR("Cannot initialize MSCKF VIO...");
    return;
  }
  return;
}

PLUGINLIB_EXPORT_CLASS(msckf_vio::SekfseEMNodelet, nodelet::Nodelet);

} // end namespace msckf_vio