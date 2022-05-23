#include "bmp.h"		//	Simple .bmp library
#include<iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

using namespace std;

#define Baseline 30.0
#define Focal_Length 100
#define Image_Width 35.0
#define Image_Height 35.0
#define Resolution_Row 512
#define Resolution_Col 512
#define View_Grid_Row 9
#define View_Grid_Col 9

struct Point3d
{
	double x;
	double y;
	double z;
	Point3d(double x_, double y_, double z_) :x(x_), y(y_), z(z_) {}
};

struct Point2d
{
	double x;
	double y;
	Point2d(double x_, double y_) :x(x_), y(y_) {}
};


int main(int argc, char** argv)
{
	if(argc < 5 || argc > 7)
	{
		cout << "Arguments prompt: viewSynthesis.exe <LF_dir> <X Y Z> OR: viewSynthesis.exe <LF_dir> <X Y Z> <focal_length> <opacity>" << endl;
		return 0;
	}
	string LFDir = argv[1];
	double Vx = stod(argv[2]), Vy = stod(argv[3]), Vz = stod(argv[4]);
	double targetFocalLen = 100;
	if(argc >= 6)
	{
		targetFocalLen = stod(argv[5]);
	}
	

	vector<Bitmap> viewImageList;
	//! loading light field views
	for (int i = 0; i < View_Grid_Col*View_Grid_Row; i++)
	{
		char name[128];
		sprintf(name, "/cam%03d.bmp", i+1);
		string filePath = LFDir + name;
		Bitmap view_i(filePath.c_str());
		viewImageList.push_back(view_i);
	}

	Bitmap targetView(Resolution_Col, Resolution_Row);
	cout << "Synthesizing image from viewpoint: (" << Vx << "," << Vy << "," << Vz << ") with focal length: " << targetFocalLen << endl;
	//! resample pixels of the target view one by one'

			Color  vp1,  vp2,  vp3,  vp4;
			double x,y,a,b;			
			double dx, dy;
			int u1,u2,v1,v2 , c11, c12, c21, c22, r11, r12, r21 ,r22;
	for (int r = 0; r < Resolution_Row; r++)
	{
		for (int c = 0; c < Resolution_Col; c++)
		{
			Point3d rayRGB(0, 0, 0);
			dx = ( (c - ( (double) Resolution_Col - 1) / 2) * (Image_Width / Resolution_Col)) / targetFocalLen; // dx of the ray from new viewpt to pixel(c,r)
			dy = ( (r - ( (double) Resolution_Row - 1) / 2) * (Image_Height/ Resolution_Row)) / targetFocalLen; // dy of the ray from new viewpt to pixel(c,r)

			x = Vx + Vz*dx; // x-coordinates of the projection on the viewplane by the ray from new viewpt to pixel(c,r)
			y = Vy + Vz*dy; // y-coordinates of the projection on the viewplane by the ray from new viewpt to pixel(c,r)
			if ( ( x >= -Baseline*4 ) && ( x <= Baseline*4 ) && ( y >= -Baseline*4 ) && ( y <= Baseline*4 ) ) { 
			u1 = (int) floor((Baseline*4 + x) / Baseline); // view point ui
			u2 = (int) ceil((Baseline*4 + x) / Baseline); // view point ui+1
			v1 = (int) floor((Baseline*4 - y) / Baseline); // view point vi
			v2 = (int) ceil((Baseline*4 - y) / Baseline); // view point vi+1	
			
			a = (x- (-Baseline*4 + u1*Baseline) ) / Baseline; // alpha
			b = ((Baseline*4 - v1*Baseline) - y) / Baseline; // beta

			Point2d u1v1( (-Baseline*4 + u1 * Baseline) + Focal_Length * dx , (Baseline*4 - v1 * Baseline ) + Focal_Length * dy); // x,y coordinates of point u1v1
			Point2d u1v2( (-Baseline*4 + u1 * Baseline) + Focal_Length * dx , (Baseline*4 - v2 * Baseline ) + Focal_Length * dy); // x,y coordinates of point u1v2
			Point2d u2v1( (-Baseline*4 + u2 * Baseline) + Focal_Length * dx , (Baseline*4 - v1 * Baseline ) + Focal_Length * dy); // x,y coordinates of point u2v1
			Point2d u2v2( (-Baseline*4 + u2 * Baseline) + Focal_Length * dx , (Baseline*4 - v2 * Baseline ) + Focal_Length * dy); // x,y coordinates of point u2v2

			c11 = round((u1v1.x - (-120 + u1*Baseline)) * Resolution_Col / Image_Width + ( (double) Resolution_Col - 1) / 2); // pixel column of point u1v1
			c12 = round((u1v2.x - (-120 + u1*Baseline)) * Resolution_Col / Image_Width + ( (double) Resolution_Col - 1) / 2); // pixel column of point u1v2
			c21 = round((u2v1.x - (-120 + u2*Baseline)) * Resolution_Col / Image_Width + ( (double) Resolution_Col - 1) / 2); // pixel column of point u2v1
			c22 = round((u2v2.x - (-120 + u2*Baseline)) * Resolution_Col / Image_Width + ( (double) Resolution_Col - 1) / 2); // pixel column of point u2v2
			r11 = round((u1v1.y - (120 - v1*Baseline)) * Resolution_Row / Image_Height + ( (double) Resolution_Row - 1) / 2); // pixel row of point u1v1
			r12 = round((u1v2.y - (120 - v2*Baseline)) * Resolution_Row / Image_Height + ( (double) Resolution_Row - 1) / 2); // pixel row of point u1v2
			r21 = round((u2v1.y - (120 - v1*Baseline)) * Resolution_Row / Image_Height + ( (double) Resolution_Row - 1) / 2); // pixel row of point u2v1
			r22 = round((u2v2.y - (120 - v2*Baseline)) * Resolution_Row / Image_Height + ( (double) Resolution_Row - 1) / 2); // pixel row of point u2v2

			if (c11 <= Resolution_Col - 1 && c11 >= 0 && r11 <= Resolution_Row - 1 && r11 >= 0 && c12 <= Resolution_Col - 1 && c12 >= 0 && r12 <= Resolution_Row - 1 && r12 >= 0 &&c21 <= Resolution_Col - 1 && c21 >= 0 && r21 <= Resolution_Row - 1 && r21 >= 0 &&c22 <= Resolution_Col - 1 && c22 >= 0 && r22 <= Resolution_Row - 1 && r22 >= 0 ) // check if the pixels exist
			{
			viewImageList.at(u1+v1*9).getColor(c11, r11, vp1.R, vp1.G, vp1.B); // ui vi
			viewImageList.at(u2+v1*9).getColor(c12 ,r12, vp2.R, vp2.G, vp2.B); // ui+1 vi
			viewImageList.at(u1+v2*9).getColor(c21, r21, vp3.R, vp3.G, vp3.B); // ui vi+1
			viewImageList.at(u2+v2*9).getColor(c22, r22, vp4.R, vp4.G, vp4.B); // ui+1 vi+1
			
			rayRGB.x = (1-b) * ( (1-a) * (vp1.R) + a * (vp2.R) ) + b * ( (1-a) * (vp3.R) + a * (vp4.R));
			rayRGB.y = (1-b) * ( (1-a) * (vp1.G) + a * (vp2.G) ) + b * ( (1-a) * (vp3.G) + a * (vp4.G));
			rayRGB.z = (1-b) * ( (1-a) * (vp1.B) + a * (vp2.B) ) + b * ( (1-a) * (vp3.B) + a * (vp4.B));
			}

			//! record the resampled pixel value
			
			}
			targetView.setColor(c, r, (unsigned char) rayRGB.x, (unsigned char) rayRGB.y, (unsigned char) rayRGB.z);
		}
	}

	if (argc == 7) {
		double opacity = stod(argv[6]);
			Color b;
			for (int r = 0; r < Resolution_Row; r++ ) {
				for (int c = 0; c < Resolution_Col; c++) {
					targetView.getColor(c , r , b.R, b.G, b.B);
					b.R = round(b.R * opacity / 100);
					b.G = round(b.G * opacity / 100);
					b.B = round(b.B * opacity / 100);
					targetView.setColor(c , r , b.R, b.G, b.B);
				}
			}
	}
	string savePath = "newView.bmp";
	targetView.save(savePath.c_str());
	cout << "Result saved!" << endl;
	return 0;
}