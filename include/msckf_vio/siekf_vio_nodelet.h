/*
 * COPYRIGHT AND PERMISSION NOTICE
 * Penn Software SIEKF_VIO
 * Copyright (C) 2017 The Trustees of the University of Pennsylvania
 * All rights reserved.
 */

#ifndef SIEKF_MSCKF_VIO_NODELET_H
#define SIEKF_MSCKF_VIO_NODELET_H

#include <nodelet/nodelet.h>
#include <pluginlib/class_list_macros.h>
#include <msckf_vio/siekf_vio.h>

namespace msckf_vio {
class SiekfVioNodelet : public nodelet::Nodelet {
public:
  SiekfVioNodelet() { return; }
  ~SiekfVioNodelet() { return; }

private:
  virtual void onInit();
  SiekfVioPtr siekf_vio_ptr;
};
} // end namespace msckf_vio

#endif