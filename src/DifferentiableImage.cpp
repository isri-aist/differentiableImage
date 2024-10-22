#include <DifferentiableImage/DifferentiableImage.h>
#include <per/core/prOmni.h>

void DifferentiableImage::init(const std::string & path_des,
                               const std::string & path_cur,
                               bool enable_display,
                               const std::string & mask)
{
  enableDisplay = enable_display;
  try
  {
    vpImageIo::read(I_des, path_des);

    if(I_des.getHeight() == 0)
    {
      std::cerr << "Error while reading des image" << std::endl;
      return;
    }

    vpImageIo::read(I_cur, path_cur);

    if(I_cur.getHeight() == 0)
    {
      std::cerr << "Error while reading current image" << std::endl;
      return;
    }

    if(enableDisplay)
    {
      disp_I_des.init(I_des, 100, 100, "I_des (SharpCapturedDesired)");
      vpDisplay::display(I_des);
      vpDisplay::flush(I_des);

      disp_I_cur.init(I_cur, 100 + I_des.getWidth() + 5, 100, "I_cur (SharpCapturedCurrent)");
      vpDisplay::display(I_cur);
      vpDisplay::flush(I_cur);
    }

    hasMask = false;
    if(!mask.empty())
    {
      vpImageIo::read(I_mask, mask);

      if(I_mask.getHeight() == 0)
      {
        std::cerr << "Reading \"" << mask << "\" failed.\n";
        I_mask.resize(I_des.getHeight(), I_des.getWidth());
        I_mask = 255;
      }
      else
        hasMask = true;
    }
    else
    {
      I_mask.resize(I_des.getHeight(), I_des.getWidth());
      I_mask = 255;
    }

    PGM_des_u.resize(I_des.getHeight(), I_des.getWidth());
    I_to_upsub.resize(I_des.getHeight(), I_des.getWidth());
    upsampled_subsampled_I_to_upsub.resize(I_des.getHeight(), I_des.getWidth());
    diff_u.resize(I_des.getHeight(), I_des.getWidth());
    PGM_cur_u.resize(I_cur.getHeight(), I_cur.getWidth());

    if(enableDisplay)
    {
      disp_PGM_des.init(PGM_des_u, 100 + 2 * I_des.getWidth() + 5, 100, "PGM_des (GaussDesired)");
      disp_upsampled_subsampled_I_to_upsub.init(upsampled_subsampled_I_to_upsub, 100 + 3 * I_des.getWidth() + 5, 100,
                                                "upSubGauss");
      disp_diff_u.init(diff_u, 100 + 4 * I_des.getWidth() + 5, 100, "diff_u (upSubGauss-Gauss)");
      disp_PGM_cur.init(PGM_cur_u, 100 + 5 * I_cur.getWidth() + 5, 100, "PGM_cur (GaussCurrent)");
    }
  }
  catch(vpException e)
  {
    std::cerr << "Unable to load image in init : " << e << std::endl;
    return;
  }

  nbKernels = 10;
  sigmaMin = 0.5 / 3.0; // within a pixel
  sigmaMax = I_des.getWidth() / 6.0; // desired_Gray.cols/6.0 ; //the whole image
  sigma = sigmaMin;
  sigmaStep = (sigmaMax - sigmaMin) / nbKernels;
  sigmaStep /= pow(2, nbKernels-3);

  v_cost.reserve(nbKernels);
  v_sigma.reserve(nbKernels);

  // Create Camera
  _sensor = new prPerspective();
  _camera = (prPerspective *)_sensor;

  // _sensor = new prOmni();
  // _camera = (prOmni *)_sensor;
}

double DifferentiableImage::compute()
{
  std::vector<double> sigmas, costs;
  computeCosts(sigmas, costs);
  return computeInterpolatedSigma(sigmas, costs);
}

void DifferentiableImage::computeCosts(std::vector<double> & sigmas, std::vector<double> & costs)
{
  for(unsigned int kernel = 0; kernel < nbKernels; kernel++, sigma += sigmaStep)
  {
    // prepare the desired image
    prRegularlySampledCPImage<unsigned char> IP_des(I_des.getHeight(), I_des.getWidth());

    IP_des.setInterpType(prInterpType::IMAGEPLANE_BILINEAR);

    if(hasMask)
      IP_des.buildFrom(I_des, _camera, &I_mask);
    else
      IP_des.buildFrom(I_des, _camera);

    prRegularlySampledCPImage<float> GP(I_des.getHeight(),
                                        I_des.getWidth()); // contient tous les pr2DCartesianPointVec (ou
                                                           // prFeaturePoint) u_g et fera
                                                           // GS_sample.buildFrom(IP_des, u_g);
    prFeaturesSet<prCartesian2DPointVec, prPhotometricnnGMS<prCartesian2DPointVec>, prRegularlySampledCPImage> fSet_des;
    double lambda_g = sigma;
    prPhotometricnnGMS<prCartesian2DPointVec> GP_sample_des(lambda_g);

    double t0;
    t0 = vpTime::measureTimeMs();
    fSet_des.buildFrom(IP_des, GP, GP_sample_des, false, true); // Goulot !

    vpImage<float> PGM_des_f(I_des.getHeight(), I_des.getWidth());
    vpPoseVector pp;
    fSet_des.sampler.toImage(PGM_des_f, pp, _sensor);
    vpImageConvert::convert(PGM_des_f, PGM_des_u);

    if(enableDisplay)
    {
      vpDisplay::display(PGM_des_u);
      vpDisplay::flush(PGM_des_u);
    }

    prRegularlySampledCPImage<unsigned char> IP_cur(I_cur.getHeight(),
                                                    I_cur.getWidth()); // the regularly sample planar image to be set
                                                                       // from the acquired/loaded perspective image
    prRegularlySampledCPImage<float> GP_cur(I_cur.getHeight(),
                                            I_cur.getWidth()); // contient tous les pr2DCartesianPointVec (ou
                                                               // prFeaturePoint) u_g et fera
                                                               // GS_sample.buildFrom(IP_des, u_g);
    prFeaturesSet<prCartesian2DPointVec, prPhotometricnnGMS<prCartesian2DPointVec>, prRegularlySampledCPImage> fSet_cur;

    prPhotometricnnGMS<prCartesian2DPointVec> GP_sample_cur(lambda_g);
    vpImage<float> PGM_cur_f(I_cur.getHeight(), I_cur.getWidth());

    IP_cur.setInterpType(prInterpType::IMAGEPLANE_BILINEAR);
    std::cout << "kernel: " << kernel << std::endl << "sigma: " << sigma << std::endl;
    std::cout << "build cur" << std::endl;
    if(hasMask)
      IP_cur.buildFrom(I_cur, _camera, &I_mask);
    else
      IP_cur.buildFrom(I_cur, _camera);
    std::cout << "cur built" << std::endl;

    //    IP_cur.toAbsZN(); //prepare pixels intensities ?
    t0 = vpTime::measureTimeMs();
    fSet_cur.buildFrom(IP_cur, GP_cur, GP_sample_cur, false, true); // Goulot !

    fSet_cur.sampler.toImage(PGM_cur_f, pp, _sensor);
    vpImageConvert::convert(PGM_cur_f, PGM_cur_u);

    if(enableDisplay)
    {
      vpDisplay::display(PGM_cur_u);
      vpDisplay::flush(PGM_cur_u);
    }

    vpImageTools::imageDifference(PGM_cur_u, PGM_des_u, I_to_upsub);

    // subsample the image, either the desired one or current-desired if current is set
    vpImage<unsigned char> subsampled_I_to_upsub;
    I_to_upsub.subsample(2, 2, subsampled_I_to_upsub);

    ceres::Grid2D<unsigned char, 1> arrayDes(subsampled_I_to_upsub.bitmap, 0, subsampled_I_to_upsub.getHeight(), 0,
                                             subsampled_I_to_upsub.getWidth());
    ceres::BiCubicInterpolator<ceres::Grid2D<unsigned char, 1>> subsampled_desiredInterpolator(arrayDes);

    double u_subd, v_subd;
    double cost = 0;
    double I_interp;

    unsigned int nbPix = 0;
    for(unsigned int v = 0; v < I_to_upsub.getHeight(); v++)
    {
      v_subd = v * 0.5;
      for(unsigned int u = 0; u < I_to_upsub.getWidth(); u++)
      {
        u_subd = u * 0.5;

        if(I_mask[v][u] != 0)
        {
          subsampled_desiredInterpolator.Evaluate(v_subd, u_subd, &I_interp);
          upsampled_subsampled_I_to_upsub[v][u] = I_interp;

          cost += pow(I_interp - I_to_upsub[v][u], 2);

          nbPix++;
        }
      }
    }

    cost *= 0.5;

    cost /= (double)nbPix;

    dCostSigma = cost - v_cost.back();
    std::cout << "photometric cost: " << cost << " and it derivative w.r.t. sigma: " << dCostSigma << std::endl;

    if(enableDisplay)
    {
      vpDisplay::display(upsampled_subsampled_I_to_upsub);
      vpDisplay::flush(upsampled_subsampled_I_to_upsub);
    }

    vpImageTools::imageDifference(upsampled_subsampled_I_to_upsub, I_to_upsub, diff_u);

    if(enableDisplay)
    {
      vpDisplay::display(diff_u);
      vpDisplay::flush(diff_u);
    }

    v_sigma.push_back(sigma);
    v_cost.push_back(cost);

    sigmaStep *= 2;
  }

  sigmas = v_sigma;
  costs = v_cost;
}

/*
double DifferentiableImage::computeInterpolatedSigma(std::vector<double> & sigmas, std::vector<double> & costs)
{
  double e_max, e_min, val_95;
  e_max = *(costs.begin());

  e_min = *(costs.end());
  val_95 = e_min + 0.05 * (e_max - e_min); // calculate the 95% of the cost (epsilon)
  std::cout << "val_95: " << val_95 << std::endl;

  for(size_t i = 0; i < costs.size() - 1; ++i)
  {
    if((costs[i] <= val_95 && costs[i + 1] >= val_95) || (costs[i] >= val_95 && costs[i + 1] <= val_95))
    {
      double t = (val_95 - costs[i]) / (costs[i + 1] - costs[i]);
      interpolated_sigma =
          sigmas[i] + t * (sigmas[i + 1] - sigmas[i]); // interpolate to find lambda_g at 95% of the cost
      break;
    }
  }
  std::cout << "Interpolated sigma for val_95: " << interpolated_sigma << std::endl;

  return interpolated_sigma;
}
*/

double DifferentiableImage::computeInterpolatedSigma(std::vector<double> & sigmas, std::vector<double> & costs)
{
  double e_max, e_min, val_95;
  std::vector<double> invCosts;
  invCosts.reserve(costs.size());

  for(size_t i = 0; i < costs.size() - 1; ++i)
    invCosts.push_back(1./costs[i]);

  e_max = *(--(invCosts.end()));
  val_95 = 0.95 * e_max; // calculate the 95% of the cost (epsilon)
  std::cout << "inv val_95: " << val_95 << " (val_95 " << 1./val_95 << ")" << std::endl;

  for(size_t i = 0; i < invCosts.size() - 1; ++i)
  {
    if((invCosts[i] <= val_95 && invCosts[i + 1] >= val_95) || (invCosts[i] >= val_95 && invCosts[i + 1] <= val_95))
    {
      double t = (val_95 - invCosts[i]) / (invCosts[i + 1] - invCosts[i]);
      interpolated_sigma =
          sigmas[i] + t * (sigmas[i + 1] - sigmas[i]); // interpolate to find lambda_g at 95% of the cost
      break;
    }
  }
  std::cout << "Interpolated sigma for val_95: " << interpolated_sigma << std::endl;

  return interpolated_sigma;
}