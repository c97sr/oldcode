// Copyright Steven Riley (sr@stevenriley.net)

// This file is part of the library idsource.

// idsource is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This work is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this idsource.  If not, see <http://www.gnu.org/licenses/>.

using namespace std;

#include"utility.h"

class AsciiPop {
private:
	double* vals;
	double minx,miny,maxx,maxy,stepx,stepy,maxval;
	int nox,noy; // number of grid points
	void CalcMaxVal();
public:
	// nox and noy are numbers of intervals in the constructors, not number of points
	AsciiPop() {srerror("No default constructor for AsciiPop");};
	// DensityField(double d, int nox_in, int noy_in, double minx_in, double maxx_in, double miny_in, double maxy_in);
	// DensityField(double f(double,double), int nox_in, int noy_in, double minx_in, double maxx_in, double miny_in, double maxy_in);
	// DensityField(string filename, double stepx_in, double stepy_in, double minx_in, double maxx_in, double miny_in, double maxy_in);
	AsciiPop(string filename, bool report);
	// ~DensityField(){delete [] vals;};
	~AsciiPop(){};
	double Value(double x, double y);
	double Value(int x, int y);
	string Table();
	inline double GetMinX(){return minx;};
	inline double GetMaxX(){return maxx;};
	inline double GetMinY(){return miny;};
	inline double GetMaxY(){return maxy;};
	inline double GetMaxVal(){return maxval;};
	inline int GetNoX(){return nox;};
	inline int GetNoY(){return noy;};
};
