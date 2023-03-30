
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "gflags/gflags.h"
//#include "pgm_image.h"

#include <opencv2/opencv.hpp>

#include <ceres/cubic_interpolation.h>

#include <array>

DEFINE_string(desired, "", "File of the reference image");
DEFINE_bool(isEqui, false, "");

using namespace cv;

int main(int argc, char** argv) 
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  if (FLAGS_desired.empty()) {
    std::cerr << "Please provide a reference image file name using -desired.\n";
    return 1;
  }  

  // Read the image
  Mat desired_Gray = imread(FLAGS_desired, IMREAD_GRAYSCALE);
  
  if (desired_Gray.cols == 0) {
    std::cerr << "Reading \"" << FLAGS_desired << "\" failed.\n";
    return 3;
  }
  
  namedWindow("SharpCapturedDesired", WINDOW_AUTOSIZE );
  imshow("SharpCapturedDesired", desired_Gray);
  
  namedWindow("GaussDesired", WINDOW_AUTOSIZE );
  
	namedWindow("upSubGaussDesired", WINDOW_AUTOSIZE );
	
	namedWindow("upSubGaussDesired-Desired", WINDOW_AUTOSIZE );
 

  unsigned int nbKernels = 40;
  double sigmaMin=0.5/3.0; // within a pixel 
  double sigmaMax=desired_Gray.cols/6.0;//desired_Gray.cols/6.0 ; //the whole image
  double sigma=sigmaMin;
  double sigmaStep=(sigmaMax-sigmaMin)/nbKernels;
  
  std::vector<double> v_cost;
  std::vector<double> v_sigma;
  v_cost.reserve(nbKernels);
  v_sigma.reserve(nbKernels);
  
  double dCostSigma;
  
  std::ostringstream s;
  std::string filename;
  
  for(unsigned int kernel = 0 ; kernel < nbKernels ; kernel++, sigma+=sigmaStep)
  {
  	Mat Gauss_desired_Gray = desired_Gray.clone();
  	
  	if(FLAGS_isEqui)
  	{
  		// As "BORDER_WRAP is not supported." and equirectangular images are special... 
  		GaussianBlur(Gauss_desired_Gray, Gauss_desired_Gray, Size(0,0), sigma, sigma, BORDER_ISOLATED);//BORDER_WRAP);
  	}
  	else
  	{
  		GaussianBlur(Gauss_desired_Gray, Gauss_desired_Gray, Size(0,0), sigma, sigma, BORDER_REPLICATE);
  	}	
  
		Mat desired;
		Gauss_desired_Gray.convertTo(desired, CV_64F);		
		
		s.str("");
		s.setf(std::ios::right, std::ios::adjustfield);
		s << "Gauss_" << std::setfill('0') << std::setw(6) << kernel << ".jpg";
		filename = s.str();
		
		imwrite(filename, Gauss_desired_Gray);
		
		
		//subsample the image
		Mat subsampled_desired(desired.rows*0.5, desired.cols*0.5, CV_64FC1, Scalar(0.0));
		resize(desired, subsampled_desired, subsampled_desired.size(), 0, 0, INTER_NEAREST);

		ceres::Grid2D<double, 1> arrayDes(subsampled_desired.ptr<double>(0), 0, subsampled_desired.rows, 0, subsampled_desired.cols);
		ceres::BiCubicInterpolator<ceres::Grid2D<double, 1> > subsampled_desiredInterpolator(arrayDes);
		
		cv::Mat upsampled_subsampled_desired = desired.clone();
		cv::Mat diff = desired.clone();
		double u_subd, v_subd;
		double cost = 0;
		double I_interp;

		for (unsigned int v = 0 ; v < desired.rows ; v++)
		{
			v_subd = v*0.5;
		  for (unsigned int u = 0 ; u < desired.cols ; u++)
		  {
		  	u_subd = u*0.5;
		  	
		  	if(desired.at<double>(v, u) != 255)
		  	{
		  	
					subsampled_desiredInterpolator.Evaluate(v_subd, u_subd, &I_interp);
					
					upsampled_subsampled_desired.at<double>(v, u) = I_interp;
					
					diff.at<double>(v, u) = I_interp - desired.at<double>(v, u);
					
					cost += pow(diff.at<double>(v, u), 2);
				}
		    
		  }
		}
		
		cost *= 0.5;
		dCostSigma = cost-v_cost.back();
		std::cout << "photometric cost: " << cost << " and it derivative w.r.t. sigma: " << dCostSigma << std::endl;
	

		//displays  
		upsampled_subsampled_desired.convertTo(upsampled_subsampled_desired, CV_8UC1);
		

		imshow("upSubGaussDesired", upsampled_subsampled_desired);
		
		s.str("");
		s.setf(std::ios::right, std::ios::adjustfield);
		s << "Interp_" << std::setfill('0') << std::setw(6) << kernel << ".jpg";
		filename = s.str();
		
		imwrite(filename, upsampled_subsampled_desired);
		
		//diff is in [-255 ; 255]	so...
		double m = 10.0;
		double shift = 127.5;
		diff.convertTo(diff, CV_8UC1, m, shift);		

		imshow("upSubGaussDesired-Desired", diff);
		
		s.str("");
		s.setf(std::ios::right, std::ios::adjustfield);
		s << "diff_" << std::setfill('0') << std::setw(6) << kernel << ".jpg";
		filename = s.str();
		
		imwrite(filename, diff);
		
		waitKey(1);  
		
		/*if( (v_cost.size() > 1) && (dCostSigma > 0) )
			break;
		*/
		
		v_sigma.push_back(sigma);
		v_cost.push_back(cost);
	
	}
		
	
	s.str("");
  s.setf(std::ios::right, std::ios::adjustfield);
  s << "residuals.txt";
  filename = s.str();
  std::ofstream ficCost(filename.c_str());
  std::vector<double>::iterator it_cost = v_cost.begin();
  for(;it_cost != v_cost.end() ; it_cost++)
  {
      ficCost << *it_cost << std::endl;
  }
  ficCost.close();
  
  s.str("");
  s.setf(std::ios::right, std::ios::adjustfield);
  s << "sigma.txt";
  filename = s.str();
  std::ofstream ficSigma(filename.c_str());
  std::vector<double>::iterator it_sigma = v_sigma.begin();
  for(;it_sigma != v_sigma.end() ; it_sigma++)
  {
      ficSigma << *it_sigma << std::endl;
  }
  ficSigma.close();

  return 0;
}

