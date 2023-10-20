

vector< vector<FeatureIDType> > vec_feature_ids
// Process the features which lose track.
  for (const auto& feature_id : processed_feature_ids) {
    auto& feature = map_server[feature_id];

    for (const auto& measurement : feature.observations){
	int cam_sequence = std::distance(state_server.cam_states.begin(),
        state_server.cam_states.find(measurement.first));
	vec_feature_ids[cam_sequence] = feature_id;
    }
  }



      cam_state_ids.push_back(measurement.first);


    FeatureIDType

    MatrixXd H_xj;
    VectorXd r_j;
    featureJacobian(feature.id, cam_state_ids, H_xj, r_j);

    if (gatingTest(H_xj, r_j, cam_state_ids.size()-1)) {
      H_x.block(stack_cntr, 0, H_xj.rows(), H_xj.cols()) = H_xj;
      r.segment(stack_cntr, r_j.rows()) = r_j;
      stack_cntr += H_xj.rows();
    }

    // Put an upper bound on the row size of measurement Jacobian,
    // which helps guarantee the executation time.
    if (stack_cntr > 1500) break;
  }

   = std::distance(state_server.cam_states.begin(),
        state_server.cam_states.find(cam_id));
