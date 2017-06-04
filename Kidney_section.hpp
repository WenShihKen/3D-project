#include"IncludesAndDefines.hpp"

class Kidney_section {
public:
	vector<Point>AllImageCenter;//All Image layer mass center
	vector<Mat>Original_Data;//Original_Data
	vector<Mat>ImageContour;//Use findcontour function to get contour
	vector<Mat>ResultContour;//The contours of each layer, including 360 degree
	vector<double>ContoursAngle;//Store the contours of rotate angle, we use 0~39 degree and 141~179 degree
	Point Average_Center;
	int MaxAreaImage = 0;//Which image have the maximum area
	clock_t t1;
public:
	void Get_contour(int, int, double);
	Kidney_section(string, unsigned int);//Constructor, read image and get Origin image data
	~Kidney_section();//Clear all vector
	void GetCenter(unsigned int);//Get All Image layer mass center and Average_Center
	void ParallelOnGetContour(double);//Get the contours of each layer
	Mat rotateImage(Mat, double, Point);//Rotate function
};
