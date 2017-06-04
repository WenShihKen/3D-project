#include"IncludesAndDefines.hpp"
#include"Kidney_section.hpp"

Kidney_section::Kidney_section(string FilePath, unsigned int File_Amount) {
	//read file
	for (unsigned int i = 0; i < File_Amount; i++) {
		string filename;
		stringstream ss;
		ss << (i + 1);
		filename += ss.str();
		Mat image = imread(FilePath + filename + ".tif", 1);
		Original_Data.push_back(image);
	}
	//end read file
}

Kidney_section::~Kidney_section() {
	AllImageCenter.clear();
	Original_Data.clear();
	ImageContour.clear();
	ContoursAngle.clear();
	ResultContour.clear();
}

void Kidney_section::GetCenter(unsigned int File_Amount) {
	double maxarea = 0.0;
	for (unsigned int i = 0; i < File_Amount; i++) {
		Mat image;
		Original_Data[i].copyTo(image);
		if (image.empty()) {
			break;
		}
		else if (image.channels() > 1) {
			cvtColor(image, image, COLOR_BGR2GRAY);
		}
		blur(image, image, Size(3, 3));
		threshold(image, image, 25, 255, THRESH_BINARY);
		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;
		findContours(image, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);

		//Find max area
		double temp_area = 0.0;
		int max_contour_Area = 0;
		for (int i = 0; i != contours.size(); i++) {
			if (temp_area < contourArea(contours[i])) {
				max_contour_Area = i;
				temp_area = max(temp_area, contourArea(contours[i]));
			}
		}
		if (maxarea < temp_area) {
			maxarea = temp_area;
			MaxAreaImage = i;
		}

		/// Get the moments
		vector<Moments> mu(contours.size());
		for (unsigned int i = 0; i < contours.size(); i++) {
			mu[i] = moments(contours[i], false);
		}

		///  Get the mass centers:
		vector<Point> mc(contours.size());
		for (unsigned int i = 0; i < contours.size(); i++) {
			mc[i] = Point(mu[i].m10 / mu[i].m00, mu[i].m01 / mu[i].m00);
		}
		/// Draw contours
		Scalar color = Scalar(255, 255, 255);
		drawContours(image, contours, max_contour_Area, color, 2, 8, hierarchy, 0, Point());

		AllImageCenter.push_back(mc[max_contour_Area]);
		ImageContour.push_back(image);
	}
	for (unsigned int i = 0; i < AllImageCenter.size(); i++) {
		Average_Center += AllImageCenter[i];
	}
	Average_Center /= int(AllImageCenter.size());
}

// get contour works here(it will be re-assign to each thread below)
void Kidney_section::Get_contour(int begin, int end, double thickness) {
	for (int degree = begin; degree < end; degree++) {
		Mat result(Size(ImageContour[0].cols, ImageContour[0].rows), CV_8U, Scalar(0));

		int start_x = 100, layer = 1;
		vector < vector<Point>>left_point(1), right_point(1), all_point(1);//storing polygon point
		//left_point.resize(1), right_point.resize(1), all_point.resize(1);

		for (unsigned int i = 0; i < ImageContour.size(); i++) {
			Point start_point, end_point;//line position
			Mat image = ImageContour[i];
			//Rotate image
			image = rotateImage(image, double(degree), Average_Center);
			//End Rotate image

			int max_len = 0;

			for (int i = 0; i < image.rows; i++) {
				for (int j = 0; j < image.cols; j++) {
					if (image.ptr <uchar>(i, j)[0] == 255 && image.ptr <uchar>(i, j + 1)[0] != 255) {
						for (int final = j + 2; final < image.cols; final++) {
							if (image.ptr <uchar>(i, final)[0] == 255) {
								if (max_len < final - j) {
									start_point = Point(j, start_x + round(thickness*layer)), end_point = Point(final, start_x + round(thickness*layer));//thickness is 6.5pixel
									max_len = final - j;
								}
							}
						}
						break;
					}
				}
			}

			left_point[0].push_back(start_point), right_point[0].push_back(end_point);
			layer++;
		}
		for (signed i = 0; i != left_point[0].size(); i++) {
			all_point[0].push_back(left_point[0][i]);
		}
		for (signed i = right_point[0].size() - 1; i != 0; i--) {
			all_point[0].push_back(right_point[0][i]);
		}

		ContoursAngle.push_back(degree);
		drawContours(result, all_point, 0, Scalar(255, 255, 255), 1, 8);
		ResultContour.push_back(result);
	}
}


// parallel works on get contour
void Kidney_section::ParallelOnGetContour(double thickness) {
	vector<thread> some_threads;
	some_threads.push_back(std::thread(&Kidney_section::Get_contour, this, 0, 20, thickness));
	some_threads.push_back(std::thread(&Kidney_section::Get_contour, this, 21, 39, thickness));
	some_threads.push_back(std::thread(&Kidney_section::Get_contour, this, 141, 160, thickness));
	some_threads.push_back(std::thread(&Kidney_section::Get_contour, this, 161, 179, thickness));
	for (auto& t : some_threads) t.join();
	some_threads.clear();
}

Mat Kidney_section::rotateImage(Mat src, double angle, Point center) {
	Mat rot_mat = getRotationMatrix2D(center, angle, 1.0);
	Mat dst;
	warpAffine(src, dst, rot_mat, src.size());
	return dst;
}
