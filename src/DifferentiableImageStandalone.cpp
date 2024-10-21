#include <DifferentiableImage/DifferentiableImage.h>
#include "glog/logging.h"
#include "gflags/gflags.h"

DEFINE_string(desired, "", "File of the reference image");
DEFINE_string(current, "", "File of the current image");
DEFINE_string(mask, "", "File of the mask image");
DEFINE_bool(display, true, "Enable/Disable display");

int main(int argc, char **argv){
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    if(FLAGS_desired.empty()) {
        std::cerr << "Please provide a reference image file name using -desired.\n";
        return 1;
    }
    if(FLAGS_current.empty()) 
    {
        std::cerr << "No current image file name using -current.\n";
        return 1;
    }

    DifferentiableImage diff_image;

    if(FLAGS_mask.empty()) 
    {
        diff_image.init(FLAGS_desired, FLAGS_current, FLAGS_display);
    }
    else{
        diff_image.init(FLAGS_desired, FLAGS_current, FLAGS_display, FLAGS_mask);
    }

    double lambda_result = diff_image.compute();

    std::cout << "Computed lambda is : " << lambda_result << std::endl;

    return 0;
}