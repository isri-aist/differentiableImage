
#pragma once

#include <array>
#include <ostream>
#include <per/core/prCameraModel.h>
#include <sstream>
#include <string>
#include <vector>

#include "ceres/ceres.h"
#include <ceres/cubic_interpolation.h>

#include <per/prFeaturesSet.h>
#include <per/prPerspective.h>
#include <per/prPhotometricnnGMS.h>
#include <per/prRegularlySampledCPImage.h>

#include <visp/vpDisplayX.h>
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpImageTools.h>

class DifferentiableImage
{
public:
  DifferentiableImage() = default;
  ~DifferentiableImage() = default;

  double start(const std::string & path_des,
               const std::string & path_cur,
               const std::string & mask = "",
               const std::string & path_des_dual = "",
               const std::string & path_cur_dual = "",
               bool enable_display = true);
  void init(const std::string & path_des,
            const std::string & path_cur,
            bool enable_display = true,
            const std::string & mask = "");
  double compute();
  void computeCosts(std::vector<double> & sigmas, std::vector<double> & costs);
  double computeInterpolatedSigma(std::vector<double> & sigmas, std::vector<double> & costs);

private:
  vpImage<unsigned char> I_des;
  vpImage<unsigned char> I_cur;
  vpImage<unsigned char> I_desDual;
  vpImage<unsigned char> I_curDual;
  vpImage<unsigned char> I_mask;
  bool hasMask;
  bool enableDisplay;

  // Desired PGM
  vpImage<unsigned char> PGM_des_u;

  vpImage<unsigned char> I_to_upsub;
  vpImage<unsigned char> upsampled_subsampled_I_to_upsub;

  vpImage<unsigned char> diff_u;
  vpImage<unsigned char> PGM_cur_u;

  vpDisplayX disp_I_des;
  vpDisplayX disp_I_cur;
  vpDisplayX disp_PGM_des;
  vpDisplayX disp_upsampled_subsampled_I_to_upsub;
  vpDisplayX disp_diff_u;
  vpDisplayX disp_PGM_cur;

  unsigned int nbKernels;
  double sigmaMin; // within a pixel
  double sigmaMax; // desired_Gray.cols/6.0 ; //the whole image
  double sigma;
  double sigmaStep;
  double r;

  double dCostSigma;
  double interpolated_sigma;

  std::vector<double> v_sigma, v_cost, v_sigmaDual, v_costDual;
  std::ostringstream s;
  std::string filename;

  // Camera
  prSensorModel * _sensor;
  prCameraModel * _camera;
};
