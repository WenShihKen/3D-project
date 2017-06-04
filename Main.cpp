#include"IncludesAndDefines.hpp"
#include"Kidney_section.hpp"


void Shape_Matching(Kidney_section, Kidney_section, double&,int& ,int&);//Find match angle
void scanner(Mat origin, Mat contour, vector<vector<uchar>> &line); //scanning
void makeSTL_horizon(vector<vector<Point2f>> contour, string,int);
void makeSTL_vertical(vector<vector<Point2f>> contour, string,int);
void fixedPoint(Kidney_section &h, vector<vector<Point2f>> &AllContour);
void dtw(vector<vector<uchar>> horizon_line, vector<vector<uchar>> vertical_line, int&, int&);//DTW
vector<vector<Point2f>> Horizon_Interpolation(Kidney_section, double);//Interpolation
vector<vector<Point2f>> Vertical_Interpolation(Kidney_section);

int main()
{
	// timing
	clock_t t1;
	time_t now_time, end_time;
	now_time = time(NULL);
	// initial steps
	std::cout << "Pre processing" << endl;
	Kidney_section Horizon(Horizon_file_path, Horizon_file_amount), Vertical(Vertical_file_path, Vertical_file_amount);
	Horizon.GetCenter(Horizon_file_amount);
	Horizon.ParallelOnGetContour(horizon_thickness);

	Vertical.GetCenter(Vertical_file_amount);

	// end initial steps
	t1 = clock();
	printf("Runtime of all threads:%.4lf seconds\n", t1 / (double)(CLOCKS_PER_SEC));
	/////////////////////////////////////////////////////////////////////////////

	// start matching shape
	double match_angle = 0.0;
	int horizon_match = 0, vertical_match = 0;
	Shape_Matching(Horizon, Vertical, match_angle, horizon_match, vertical_match);
	printf("Angle: %.0lf\n", match_angle);
	// end matching shape
	t1 = clock();
	printf("Runtime of all threads:%.4lf seconds\n", t1 / (double)(CLOCKS_PER_SEC));
	/////////////////////////////////////////////////////////////////////////////

	// start Interpolation
	printf("Start Interpolation\n");
	vector<vector<Point2f>> contour;
	contour = Horizon_Interpolation(Horizon, match_angle);

	vector<vector<Point2f>> contour_vertical;
	contour_vertical = Vertical_Interpolation(Vertical);

	// end Interpolation
	t1 = clock();
	printf("Runtime of all threads:%.4lf seconds\n", t1 / (double)(CLOCKS_PER_SEC));
	/////////////////////////////////////////////////////////////////////////////
	makeSTL_horizon(contour, "3D",horizon_match);
	makeSTL_vertical(contour_vertical, "3D",vertical_match);
	end_time = time(NULL);
	std::cout << "real runtime:" << end_time - now_time << " seconds" << endl;
}

void Shape_Matching(Kidney_section h, Kidney_section v, double& match_angle, int & horizon_match, int&vertical_match) {
	double min_area = DBL_MAX;
	for (int j = 0; j != h.ResultContour.size(); j++) {//Get match angle
		Mat horizon_temp, vertical_temp;
		h.ResultContour[j].copyTo(horizon_temp);
		v.ImageContour[v.MaxAreaImage].copyTo(vertical_temp);

		vector<vector<Point>>contours1, contours2;
		vector<Vec4i>hierarchy1, hierarchy2;

		findContours(horizon_temp, contours1, hierarchy1, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
		findContours(vertical_temp, contours2, hierarchy2, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

		double area1 = 0.0, area2 = 0.0;
		for (int i = 0; i != contours1.size(); i++) {
			area1 = max(area1, contourArea(contours1[i]));
		}
		for (int i = 0; i != contours2.size(); i++) {
			area2 = max(area2, contourArea(contours2[i]));
		}
		if (abs(area1 - area2) < min_area) {
			match_angle = h.ContoursAngle[j];
			min_area = abs(area1 - area2);
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////
	Mat horizon_origin, vertical_origin;//Original data, use copyTo is to avoid the same memory address
	h.Original_Data[h.MaxAreaImage].copyTo(horizon_origin);
	v.Original_Data[v.MaxAreaImage].copyTo(vertical_origin);
	horizon_origin = h.rotateImage(horizon_origin, match_angle, h.Average_Center);//Rotate the image before DTW

	Mat horizon_contour, vertical_contour;//Contours data, use copyTo is avoiding the same memory address
	h.ImageContour[h.MaxAreaImage].copyTo(horizon_contour);
	v.ImageContour[v.MaxAreaImage].copyTo(vertical_contour);
	horizon_contour = h.rotateImage(horizon_contour, match_angle, h.Average_Center);//Rotate the image before DTW

	vector<vector<uchar>>horizon_line, vertical_line;//Storing the Horizon and Vertical image's line pixel

	scanner(horizon_origin, horizon_contour, horizon_line);
	scanner(vertical_origin, vertical_contour, vertical_line);
	dtw(horizon_line, vertical_line, horizon_match, vertical_match);
	match_angle = 0;
}

vector<vector<Point2f>> Horizon_Interpolation(Kidney_section h, double angle)
{
	vector<vector<Point2f>>AllContour(h.ImageContour.size() + 2);//plus 2 for head and bottom points
	fixedPoint(h, AllContour);

	//main part of interpolation///////////////////////////////////////////////////
	//cubic spline process
	vector<vector<double>>x1(AllContour[0].size()), z1(AllContour[0].size()), y2(AllContour[0].size()), z2(AllContour[0].size());
	for (int i = 0; i != AllContour[0].size(); i++) {
		for (int j = 0; j != AllContour.size(); j++) {
			x1[i].push_back(AllContour[j][i].x);
			z1[i].push_back(horizon_thickness*j);
			y2[i].push_back(AllContour[j][i].y);
			z2[i].push_back(horizon_thickness*j);
		}
	}

	spline *spline1 = new spline[(AllContour[0].size())]();
	spline *spline2 = new spline[(AllContour[0].size())]();
	for (int i = 0; i != AllContour[0].size(); i++) {
		spline1[i].set_points(z1[i], x1[i]);
	}
	for (int i = 0; i != AllContour[0].size(); i++) {
		spline2[i].set_points(z2[i], y2[i]);
	}
	//Interpolation processing
	vector<vector<Point2f>>After_interpolation((AllContour.size() - 1) * horizon_Interpolation_layer);
	for (int layer = 0; layer != AllContour.size() - 1; layer++) {
		for (int i = 0; i < horizon_Interpolation_layer; i++) {
			for (int j = 0; j != AllContour[0].size(); j++) {
				After_interpolation[layer * horizon_Interpolation_layer + i].push_back(Point2f(spline1[j](layer*horizon_thickness + i*horizon_interval), spline2[j](layer*horizon_thickness + i*horizon_interval)));
			}
		}
	}
	return After_interpolation;
}

void scanner(Mat origin, Mat contour, vector<vector<uchar>> &line)
{
	for (int i = 0; i < origin.rows; i++) {//Scan origin image
		for (int j = 0; j < origin.cols; j++) {
			if (contour.ptr <uchar>(i, j)[0] == 255 && contour.ptr <uchar>(i, j + 1)[0] != 255) {//Contour is recorded as white points
				vector<uchar>temp;//temp variable
				temp.push_back(origin.ptr <uchar>(i, j)[0]);
				for (int final = j + 2; final < origin.cols; final++) {
					temp.push_back(origin.ptr <uchar>(i, final)[0]);
					if (contour.ptr <uchar>(i, final)[0] == 255) {
						line.push_back(temp);//Storing pixel value
						break;
					}
				}
				break;
			}
		}
	}
	return;
}

void dtw(vector<vector<uchar>> horizon_line, vector<vector<uchar>> vertical_line, int &horizon_match,int &vertical_match)
{
	int min_cost = INT_MAX;
	//DTW
	for (auto i = horizon_line.size() / 3; i != horizon_line.size() * 2 / 3; i++) {
		for (auto j = vertical_line.size() / 3; j != vertical_line.size() * 2 / 3; j++) {

			int **dtw, cost = 0;
			dtw = new int*[horizon_line[i].size()]();
			for (int i2 = 0; i2 != horizon_line[i].size(); i2++) {
				dtw[i2] = new int[vertical_line[j].size()]();
			}
			for (int i2 = 1; i2 != horizon_line[i].size(); i2++) {
				dtw[i2][0] = 500000;
			}
			for (int i2 = 1; i2 != vertical_line[j].size(); i2++) {
				dtw[0][i2] = 500000;
			}
			dtw[0][0] = 0;
			for (int h_size = 1; h_size != horizon_line[i].size(); h_size++) {
				for (int v_size = 1; v_size != vertical_line[j].size(); v_size++) {
					cost = abs(horizon_line[i][h_size] - vertical_line[j][v_size]);
					int min_value = min(dtw[h_size - 1][v_size], dtw[h_size][v_size - 1]);
					min_value = min(min_value, dtw[h_size - 1][v_size - 1]);
					dtw[h_size][v_size] = min_value + cost;
				}
			}
			if (dtw[horizon_line[i].size() - 1][vertical_line[j].size() - 1] < min_cost) {
				min_cost = dtw[horizon_line[i].size() - 1][vertical_line[j].size() - 1];
				horizon_match = i;
				vertical_match = j;
			}
			for (int i2 = 0; i2 != horizon_line[i].size(); i2++) {
				delete[] dtw[i2];
			}
			delete[] dtw;
		}
	}
	printf("Horizon match position:%d same with Vertical match position:%d DTW value is:%d\n", horizon_match, vertical_match, min_cost);
	return;
}

void fixedPoint(Kidney_section &h, vector<vector<Point2f>> &AllContour) {
	for (int i = 0; i < 60; i++) {
		AllContour[0].push_back(h.AllImageCenter[0]);//Head point
		AllContour[AllContour.size() - 1].push_back(h.AllImageCenter[h.AllImageCenter.size() - 1]);//Bottom point

	}
	for (int image_size = 0; image_size != h.ImageContour.size(); image_size++) {
		Mat temp;
		h.ImageContour[image_size].copyTo(temp);
		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;
		findContours(temp, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);

		//find the max contour
		double temp_area = 0.0;
		int max_area = 0;
		for (int i = 0; i != contours.size(); i++) {
			if (temp_area < contourArea(contours[i])) {
				max_area = i;
				temp_area = contourArea(contours[i]);
			}
		}

		//find angle of vector
		double *angle_relative_to_central = new double[contours[max_area].size()]();//create to store the angle of two vector

		for (int i = 0; i != contours[max_area].size(); i++) {
			double temp_angle;
			double x_vector, y_vector;
			x_vector = contours[max_area][i].x - h.AllImageCenter[image_size].x;
			y_vector = contours[max_area][i].y - h.AllImageCenter[image_size].y;

			//
			if (x_vector == 0.0)
				temp_angle = 180 / 2.0;
			else
				temp_angle = atan(fabs(double(y_vector / x_vector)))*180.0 / CV_PI;

			// find quadrant
			if ((x_vector < 0.0) && (y_vector >= 0.0))
				temp_angle = 180 - temp_angle;
			else if ((x_vector < 0.0) && (y_vector < 0.0))
				temp_angle = 180 + temp_angle;
			else if ((x_vector >= 0.0) && (y_vector < 0.0))
				temp_angle = 180 * 2.0 - temp_angle;
			angle_relative_to_central[i] = temp_angle;
		}

		//make point per 6 degree
		for (int i = 0; i < 360; i += 6) {
			double min_d = DBL_MAX;
			int match = 0;
			for (int j = 0; j != contours[max_area].size(); j++) {
				if (abs(angle_relative_to_central[j] - i) < min_d) {
					match = j;
					min_d = abs(angle_relative_to_central[j] - i);
				}
			}
			AllContour[image_size + 1].push_back(contours[max_area][match]);
		}

		delete[] angle_relative_to_central;
	}
}

//here will be change to use 3d delaunay triangulation for constructing 3D model
void makeSTL_horizon(vector<vector<Point2f>> contour, string filename,int horizon_match) {
	fstream file;
	file.open(filename + ".stl", ios::out);//You can modify stl file name
	auto mod = contour[0].size();//Avoid overflow
	file << "solid Default\n";

	float dis_x = FLT_MAX, dis_y = FLT_MAX;

	for (int i = 0; i != contour.size(); i++) {
		for (int j = 0; j != contour[i].size(); j++) {
			dis_x = min(dis_x, contour[i][j].x);
			dis_y = min(dis_y, contour[i][j].y);
		}
	}
	for (int i = 0; i != contour.size(); i++) {
		for (int j = 0; j != contour[i].size(); j++) {
			contour[i][j].x -= dis_x;
			contour[i][j].y -= dis_y;
		}
	}
	//head point
	for (int i = 0; i != contour[0].size(); i++) {
		file << "facet normal 0.000000 0.000000 0.000000\nouter loop\n";

		file << "vertex " << fixed << setprecision(6) << contour[1][i % mod].x << " " << fixed << setprecision(6) << contour[1][i % mod].y << " " << fixed << setprecision(6) << horizon_interval << "\n";

		file << "vertex " << fixed << setprecision(6) << contour[0][0].x << " " << fixed << setprecision(6) << contour[0][0].y << " " << fixed << setprecision(6) << 0 << "\n";

		file << "vertex " << fixed << setprecision(6) << contour[1][(i + 1) % mod].x << " " << fixed << setprecision(6) << contour[1][(i + 1) % mod].y << " " << fixed << setprecision(6) << horizon_interval << "\n";

		file << "endloop\nendfacet\n";
	}

	//processing
	for (int i = 1; i != contour.size() - 2; i++) {
		for (int j = 0; j != contour[i].size(); j++) {
			//The stl file method will explain in Readme
			file << "facet normal 0.000000 0.000000 0.000000\nouter loop\n";
			file << "vertex " << fixed << setprecision(6) << contour[i][j%mod].x << " " << fixed << setprecision(6) << contour[i][j%mod].y << " " << fixed << setprecision(6) << i*horizon_interval << "\n";

			file << "vertex " << fixed << setprecision(6) << contour[i][(j + 2) % mod].x << " " << fixed << setprecision(6) << contour[i][(j + 2) % mod].y << " " << fixed << setprecision(6) << i*horizon_interval << "\n";

			file << "vertex " << fixed << setprecision(6) << contour[i + 1][j%mod].x << " " << fixed << setprecision(6) << contour[i + 1][j%mod].y << " " << fixed << setprecision(6) << (i + 1)*horizon_interval << "\n";

			file << "endloop\nendfacet\n";

			file << "facet normal 0.000000 0.000000 0.000000\nouter loop\n";
			file << "vertex " << fixed << setprecision(6) << contour[i][(j + 2) % mod].x << " " << fixed << setprecision(6) << contour[i][(j + 2) % mod].y << " " << fixed << setprecision(6) << i*horizon_interval << "\n";

			file << "vertex " << fixed << setprecision(6) << contour[i + 1][(j + 2) % mod].x << " " << fixed << setprecision(6) << contour[i + 1][(j + 2) % mod].y << " " << fixed << setprecision(6) << (i + 1)*horizon_interval << "\n";

			file << "vertex " << fixed << setprecision(6) << contour[i + 1][j%mod].x << " " << fixed << setprecision(6) << contour[i + 1][j%mod].y << " " << fixed << setprecision(6) << (i + 1)*horizon_interval << "\n";
			file << "endloop\nendfacet\n";
		}
	}

	//bottom point
	for (int i = 0; i != contour[0].size(); i++) {
		file << "facet normal 0.000000 0.000000 0.000000\nouter loop\n";
		file << "vertex " << fixed << setprecision(6) << contour[contour.size() - 1][0].x << " " << fixed << setprecision(6) << contour[contour.size() - 1][0].y << " " << fixed << setprecision(6) << (contour.size() - 1) * horizon_interval << "\n";

		file << "vertex " << fixed << setprecision(6) << contour[contour.size() - 2][i % mod].x << " " << fixed << setprecision(6) << contour[contour.size() - 2][i % mod].y << " " << fixed << setprecision(6) << (contour.size() - 2)*horizon_interval << "\n";

		file << "vertex " << fixed << setprecision(6) << contour[contour.size() - 2][(i + 1) % mod].x << " " << fixed << setprecision(6) << contour[contour.size() - 2][(i + 1) % mod].y << " " << fixed << setprecision(6) << (contour.size() - 2)*horizon_interval << "\n";

		file << "endloop\nendfacet\n";
	}
	//file << "endsolid Default\n";
	//end file
	file.close();
}

vector<vector<Point2f>>Vertical_Interpolation(Kidney_section v)
{
	vector<vector<Point2f>>AllContour(v.ImageContour.size() + 2);//plus 2 for head and bottom points

	fixedPoint(v, AllContour);

	//main part of interpolation///////////////////////////////////////////////////
	//cubic spline process
	vector<vector<double>>x1(AllContour[0].size()), z1(AllContour[0].size()), y2(AllContour[0].size()), z2(AllContour[0].size());
	for (int i = 0; i != AllContour[0].size(); i++) {
		for (int j = 0; j != AllContour.size(); j++) {
			x1[i].push_back(AllContour[j][i].x);
			z1[i].push_back(vertical_thickness*j);
			y2[i].push_back(AllContour[j][i].y);
			z2[i].push_back(vertical_thickness*j);
		}
	}

	spline *spline1 = new spline[(AllContour[0].size())]();
	spline *spline2 = new spline[(AllContour[0].size())]();
	for (int i = 0; i != AllContour[0].size(); i++) {
		spline1[i].set_points(z1[i], x1[i]);
	}
	for (int i = 0; i != AllContour[0].size(); i++) {
		spline2[i].set_points(z2[i], y2[i]);
	}
	//Interpolation processing
	vector<vector<Point2f>>After_interpolation((AllContour.size() - 1) * vertical_Interpolation_layer);
	for (int layer = 0; layer != AllContour.size() - 1; layer++) {
		for (int i = 0; i < vertical_Interpolation_layer; i++) {
			for (int j = 0; j != AllContour[0].size(); j++) {
				After_interpolation[layer * vertical_Interpolation_layer + i].push_back(Point2f(spline1[j](layer*vertical_thickness + i*vertical_interval), spline2[j](layer*vertical_thickness + i*vertical_interval)));
			}
		}
	}
	return After_interpolation;
}

void makeSTL_vertical(vector<vector<Point2f>> contour, string filename, int vertical_match) {
	fstream file;
	file.open(filename + ".stl", ios::app);//You can modify stl file name
	auto mod = contour[0].size();//Avoid overflow

	file << "solid Default\n";
	float dis_z = FLT_MAX, dis_x = FLT_MAX;
	for (int i = 0; i != contour.size(); i++) {
		for (int j = 0; j != contour[i].size(); j++) {
			dis_z = min(dis_z, contour[i][j].y);
			dis_x = min(dis_x, contour[i][j].x);
		}
	}
	for (int i = 0; i != contour.size(); i++) {
		for (int j = 0; j != contour[i].size(); j++) {
			contour[i][j].y -= dis_z;
			contour[i][j].x -= dis_x;
		}
	}
	//head point
	for (int i = 0; i != contour[0].size(); i++) {
		file << "facet normal 0.000000 0.000000 0.000000\nouter loop\n";
		file << "vertex " << fixed << setprecision(6) << contour[0][0].x << " " << fixed << setprecision(6) << 0 << " " << fixed << setprecision(6) << contour[0][0].y << "\n";

		file << "vertex " << fixed << setprecision(6) << contour[1][i % mod].x << " " << fixed << setprecision(6) << vertical_interval << " " << fixed << setprecision(6) << contour[1][i % mod].y << "\n";

		file << "vertex " << fixed << setprecision(6) << contour[1][(i + 1) % mod].x << " " << fixed << setprecision(6) << vertical_interval << " " << fixed << setprecision(6) << contour[1][(i + 1) % mod].y << "\n";

		file << "endloop\nendfacet\n";
	}

	//processing
	for (int i = 1; i != contour.size() - 2; i++) {
		for (int j = 0; j != contour[i].size(); j++) {
			//The stl file method will explain in Readme
			file << "facet normal 0.000000 0.000000 0.000000\nouter loop\n";

			file << "vertex " << fixed << setprecision(6) << contour[i][(j + 2) % mod].x << " " << fixed << setprecision(6) << i*vertical_interval << " " << fixed << setprecision(6) << contour[i][(j + 2) % mod].y << "\n";

			file << "vertex " << fixed << setprecision(6) << contour[i][j%mod].x << " " << fixed << setprecision(6) << i*vertical_interval << " " << fixed << setprecision(6) << contour[i][j%mod].y << "\n";

			file << "vertex " << fixed << setprecision(6) << contour[i + 1][j%mod].x << " " << fixed << setprecision(6) << (i + 1)*vertical_interval << " " << fixed << setprecision(6) << contour[i + 1][j%mod].y << "\n";

			file << "endloop\nendfacet\n";

			file << "facet normal 0.000000 0.000000 0.000000\nouter loop\n";

			file << "vertex " << fixed << setprecision(6) << contour[i + 1][(j + 2) % mod].x << " " << fixed << setprecision(6) << (i + 1)*vertical_interval << " " << fixed << setprecision(6) << contour[i + 1][(j + 2) % mod].y << "\n";

			file << "vertex " << fixed << setprecision(6) << contour[i][(j + 2) % mod].x << " " << fixed << setprecision(6) << i*vertical_interval << " " << fixed << setprecision(6) << contour[i][(j + 2) % mod].y << "\n";

			file << "vertex " << fixed << setprecision(6) << contour[i + 1][j%mod].x << " " << fixed << setprecision(6) << (i + 1)*vertical_interval << " " << fixed << setprecision(6) << contour[i + 1][j%mod].y << "\n";
			file << "endloop\nendfacet\n";
		}
	}

	//bottom point
	for (int i = 0; i != contour[0].size(); i++) {
		file << "facet normal 0.000000 0.000000 0.000000\nouter loop\n";
		file << "vertex " << fixed << setprecision(6) << contour[contour.size() - 2][i % mod].x << " " << fixed << setprecision(6) << (contour.size() - 2) * vertical_interval << " " << fixed << setprecision(6) << contour[contour.size() - 2][i % mod].y << "\n";

		file << "vertex " << fixed << setprecision(6) << contour[contour.size() - 1][0].x << " " << fixed << setprecision(6) << (contour.size() - 1) * vertical_interval << " " << fixed << setprecision(6) << contour[contour.size() - 1][0].y << "\n";

		file << "vertex " << fixed << setprecision(6) << contour[contour.size() - 2][(i + 1) % mod].x << " " << fixed << setprecision(6) << (contour.size() - 2) * vertical_interval << " " << fixed << setprecision(6) << contour[contour.size() - 2][(i + 1) % mod].y << "\n";



		file << "endloop\nendfacet\n";
	}

	//end file
	file << "endsolid Default\n";
	file.close();
	cout << "stl file success\n";
}