/*
 * COPYRIGHT AND PERMISSION NOTICE
 * Penn Software SUKFLG_VIO
 * Copyright (C) 2017 The Trustees of the University of Pennsylvania
 * All rights reserved.
 */

#ifndef SUKFLG_MSCKF_VIO_NODELET_H
#define SUKFLG_MSCKF_VIO_NODELET_H

#include <nodelet/nodelet.h>
#include <pluginlib/class_list_macros.h>
#include <msckf_vio/sukflg_vio.h>

namespace msckf_vio {
class SukflgVioNodelet : public nodelet::Nodelet {
public:
  SukflgVioNodelet() { return; }
  ~SukflgVioNodelet() { return; }

private:
  virtual void onInit();
  SukflgVioPtr sukflg_vio_ptr;
};
} // end namespace msckf_vio

#endif