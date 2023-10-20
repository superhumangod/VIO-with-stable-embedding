/*
 * COPYRIGHT AND PERMISSION NOTICE
 * Penn Software SEKFSE_EM
 * Copyright (C) 2017 The Trustees of the University of Pennsylvania
 * All rights reserved.
 */

#ifndef SEKFSE_EM_NODELET_H
#define SEKFSE_EM_NODELET_H

#include <nodelet/nodelet.h>
#include <pluginlib/class_list_macros.h>
#include <msckf_vio/sekfse_em.h>

namespace msckf_vio {
class SekfseEMNodelet : public nodelet::Nodelet {
public:
  SekfseEMNodelet() { return; }
  ~SekfseEMNodelet() { return; }

private:
  virtual void onInit();
  SekfseEMPtr sekfse_em_ptr;
};
} // end namespace msckf_vio

#endif