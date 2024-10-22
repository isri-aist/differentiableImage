#include "gflags/gflags.h"
#include "glog/logging.h"
#include <DifferentiableImage/DifferentiableImage.h>

DEFINE_string(desired, "", "File of the reference image");
DEFINE_string(current, "", "File of the current image");
DEFINE_string(desiredDual, "", "File of the reference dual image");
DEFINE_string(currentDual, "", "File of the current dual image");
DEFINE_string(mask, "", "File of the mask image");
DEFINE_bool(display, true, "Enable/Disable display");

int main(int argc, char ** argv)
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  if(FLAGS_desired.empty())
  {
    std::cerr << "Please provide a reference image file name using -desired.\n";
    return 1;
  }
  if(FLAGS_current.empty())
  {
    std::cerr << "No current image file name using -current.\n";
    return 1;
  }

  if(!FLAGS_desiredDual.empty())
  {
    std::cout << "Reference dual image file provided.\n";
    if(FLAGS_currentDual.empty())
    {
      std::cerr << "But no current one.\n";
      return 1;
    }
  }

  if(!FLAGS_currentDual.empty())
  {
    std::cout << "Current dual image file provided.\n";
    if(FLAGS_desiredDual.empty())
    {
      std::cerr << "but no desired one.\n";
      return 1;
    }
  }

  double lambda_result;

  if(FLAGS_desiredDual.empty())
  {
    DifferentiableImage diff_image;

    if(FLAGS_mask.empty()) { diff_image.init(FLAGS_desired, FLAGS_current, FLAGS_display); }
    else { diff_image.init(FLAGS_desired, FLAGS_current, FLAGS_display, FLAGS_mask); }

    lambda_result = diff_image.compute();
  }
  else
  {
    std::vector<double> v_sigma, v_cost, v_sigmaDual, v_costDual;

    {
      DifferentiableImage diff_image;

      if(FLAGS_mask.empty()) { diff_image.init(FLAGS_desired, FLAGS_current, FLAGS_display); }
      else { diff_image.init(FLAGS_desired, FLAGS_current, FLAGS_display, FLAGS_mask); }

      diff_image.computeCosts(v_sigma, v_cost);
    }

    DifferentiableImage diff_imageDual;

    if(FLAGS_mask.empty()) { diff_imageDual.init(FLAGS_desiredDual, FLAGS_currentDual, FLAGS_display); }
    else { diff_imageDual.init(FLAGS_desiredDual, FLAGS_currentDual, FLAGS_display, FLAGS_mask); }

    diff_imageDual.computeCosts(v_sigmaDual, v_costDual);

    std::vector<double>::iterator it_s = v_sigma.begin();
    std::vector<double>::iterator it_c = v_cost.begin();
    std::vector<double>::iterator it_cD = v_costDual.begin();
    for(; it_s != v_sigma.end(); it_s++, it_c++, it_cD++)
      *it_c = (*it_c + *it_cD)*0.5;

    lambda_result = diff_imageDual.computeInterpolatedSigma(v_sigma, v_cost);
  }

  std::cout << "Computed lambda is : " << lambda_result << std::endl;

  return 0;
}
