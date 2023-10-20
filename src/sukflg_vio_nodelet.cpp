/*
 * COPYRIGHT AND PERMISSION NOTICE
 * Penn Software MSCKF_VIO
 * Copyright (C) 2017 The Trustees of the University of Pennsylvania
 * All rights reserved.
 */

#include <msckf_vio/sukflg_vio_nodelet.h>

namespace msckf_vio {
void SukflgVioNodelet::onInit() {
  sukflg_vio_ptr.reset(new SukflgVio(getPrivateNodeHandle()));
  if (!sukflg_vio_ptr->initialize()) {
    ROS_ERROR("Cannot initialize MSCKF VIO...");
    return;
  }
  return;
}

PLUGINLIB_EXPORT_CLASS(msckf_vio::SukflgVioNodelet, nodelet::Nodelet);

} // end namespace msckf_vio