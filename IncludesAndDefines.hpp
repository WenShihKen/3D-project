#pragma once
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgcodecs.hpp"
#include <vector>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <thread>
#include "spline.h"

//Can modify define value
#define horizon_thickness 6.5//thickness of horizon image layer 
#define vertical_thickness 8.4//thickness of vertical image layer 
#define horizon_Interpolation_layer 10//Interpolate 10 layer
#define vertical_Interpolation_layer 10//Interpolate 10 layer
#define horizon_interval (horizon_thickness/horizon_Interpolation_layer)//the horizonInterpolation layer thickness
#define vertical_interval (vertical_thickness/vertical_Interpolation_layer)//the vertical Interpolation layer thickness
#define Vertical_file_path "data2/vertical/"
#define Vertical_file_amount 13
#define Horizon_file_path "data2/horizon/"
#define Horizon_file_amount 25


using namespace cv;
using namespace std;
using namespace bspline;
