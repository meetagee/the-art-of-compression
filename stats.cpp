
#include "stats.h"

stats::stats(PNG & im){

/* your code here */
	height = im.height();
	width = im.width();

	sumHueX.resize(height,vector<double>(width));
	sumHueY.resize(height,vector<double>(width));
	sumSat.resize(height,vector<double>(width));
	sumLum.resize(height,vector<double>(width));
	hist.resize(height,vector<vector<int>>(width));

	for(unsigned int i = 0; i < im.height(); i++){
		//cout << "i is: " << i << endl;
        for(unsigned int j = 0; j < im.width(); j++){        //right bottom pixel [i][j]
            sumHueX[i][j] = 0.0;
            sumHueY[i][j] = 0.0;
            sumSat[i][j] = 0.0;
            sumLum[i][j] = 0.0;
            for(int k = 0; k < 36; k++){            //set everything to 0
                hist[i][j].push_back(0);
            }
			if(i > 0 && j > 0){
				HSLAPixel * pixel = im.getPixel(j,i);
				sumHueX[i][j] = cos(pixel->h * PI/180) + sumHueX[i][j-1] + sumHueX[i-1][j] - sumHueX[i-1][j-1];
				sumHueY[i][j] = sin(pixel->h * PI/180) + sumHueY[i][j-1] + sumHueY[i-1][j] - sumHueY[i-1][j-1];
				sumSat[i][j] = pixel->s + sumSat[i][j-1] + sumSat[i-1][j] - sumSat[i-1][j-1];
				sumLum[i][j] = pixel->l + sumLum[i][j-1] + sumLum[i-1][j] - sumLum[i-1][j-1];
				for(int k = 0; k < 36; k++){
					if(k*10 <= pixel->h &&  pixel->h < (k+1)*10){
						hist[i][j][k] = hist[i][j-1][k] + hist[i-1][j][k] - hist[i-1][j-1][k] + 1;
					}
					else
						hist[i][j][k] = hist[i][j-1][k] + hist[i-1][j][k] - hist[i-1][j-1][k];
				}
			}
			else if(i > 0 && j == 0) {
				HSLAPixel * pixel = im.getPixel(j,i);
				sumHueX[i][j] = cos(pixel->h * PI/180) + sumHueX[i-1][j];
				sumHueY[i][j] = sin(pixel->h * PI/180) + sumHueY[i-1][j];
				sumSat[i][j] = pixel->s + sumSat[i-1][j];
				sumLum[i][j] = pixel->l + sumLum[i-1][j];
				for(int k = 0; k < 36; k++){
					if(k*10 <= pixel->h &&  pixel->h < (k+1)*10){
						hist[i][j][k] = hist[i-1][j][k] + 1;
					}
					else
						hist[i][j][k] = hist[i-1][j][k];
				}
			}
			else if(i == 0 && j > 0){
				HSLAPixel * pixel = im.getPixel(j,i);
				sumHueX[i][j] = cos(pixel->h * PI/180) + sumHueX[i][j-1];
				sumHueY[i][j] = sin(pixel->h * PI/180) + sumHueY[i][j-1];
				sumSat[i][j] = pixel->s + sumSat[i][j-1];
				sumLum[i][j] = pixel->l + sumLum[i][j-1];
				for(int k = 0; k < 36; k++){
					if(k*10 <= pixel->h &&  pixel->h < (k+1)*10){
						hist[i][j][k] = hist[i][j-1][k] + 1;
					}
					else
						hist[i][j][k] = hist[i][j-1][k];
				}
			}
			else if(i == 0 && j ==0){
				HSLAPixel * pixel = im.getPixel(j,i);
				sumHueX[i][j] = cos(pixel->h * PI/180);
				sumHueY[i][j] = sin(pixel->h * PI/180);
				sumSat[i][j] = pixel->s;
				sumLum[i][j] = pixel->l;
				for(int k = 0; k < 36; k++){
					if(k*10 <= pixel->h &&  pixel->h < (k+1)*10){
						hist[i][j][k]++;
					}
				}
			}
        }
    }
	//cout << "Finished init" << endl;
}

long stats::rectArea(pair<int,int> ul, pair<int,int> lr){
/* your code here */
		int x1 = ul.first;
		int y1 = ul.second;
		int x2 = lr.first;
		int y2 = lr.second;
		if(x2 >= x1 && y2 >= y1){
			return (x2 - x1 + 1)*(y2 - y1 + 1);
		}
		if(x2 >= x1 && y2 < y1){
			long area1 = (x2 - x1 + 1)*(height - y1);
			long area2 = (x2 - x1 + 1)*(y2 + 1);
			return area1 + area2;
		}
		if(x2 < x1 && y2 >= y1){
			long area1 = (width - x1)*(y2 - y1 + 1);
			long area2 = (x2 + 1)*(y2 - y1 + 1);
			return area1 + area2;
		}
		if(x2 < x1 && y2 < y1){
			long area1 = (x2 + 1)*(y2 + 1);		//up left
			long area2 = (width - x1)*(y2 + 1);				//up right
			long area3 = (x2 + 1)*(height - y1); // low left
			long area4 = (width - x1)*(height - y1); //low right
			return area1 + area2 + area3 + area4;
		}
        return long();
}

HSLAPixel stats::getAvg(pair<int,int> ul, pair<int,int> lr){

/* your code here */
    int x1 = ul.first;
    int y1 = ul.second;
    int x2 = lr.first;
    int y2 = lr.second;
    long pixelNum = rectArea(ul,lr);
    double avgHueX = 0;
    double avgHueY = 0;
    double avgSat = 0;
    double avgLum = 0;
    if(x1 > 0 && y1 > 0){
        avgHueX = (sumHueX[y2][x2] - sumHueX[y2][x1-1] - sumHueX[y1-1][x2] + sumHueX[y1-1][x1-1])/pixelNum;
        avgHueY = (sumHueY[y2][x2] - sumHueY[y2][x1-1] - sumHueY[y1-1][x2] + sumHueY[y1-1][x1-1])/pixelNum;
        avgSat = (sumSat[y2][x2] - sumSat[y2][x1-1] - sumSat[y1-1][x2] + sumSat[y1-1][x1-1])/pixelNum;
        avgLum = (sumLum[y2][x2] - sumLum[y2][x1-1] - sumLum[y1-1][x2] + sumLum[y1-1][x1-1])/pixelNum;
    }
    else if(x1 == 0 && y1 > 0){
        avgHueX = (sumHueX[y2][x2]  - sumHueX[y1-1][x2])/pixelNum;
        avgHueY = (sumHueY[y2][x2]  - sumHueY[y1-1][x2])/pixelNum;
        avgSat = (sumSat[y2][x2] - sumSat[y1-1][x2])/pixelNum;
        avgLum = (sumLum[y2][x2] - sumLum[y1-1][x2])/pixelNum;
    }
    else if(x1 > 0 && y1 == 0){
        avgHueX = (sumHueX[y2][x2]  - sumHueX[y2][x1-1])/pixelNum;
        avgHueY = (sumHueY[y2][x2]  - sumHueY[y2][x1-1])/pixelNum;
        avgSat = (sumSat[y2][x2] - sumSat[y2][x1-1])/pixelNum;
        avgLum = (sumLum[y2][x2] - sumLum[y2][x1-1])/pixelNum;
    }
    else if(x1 == 0 && y1 == 0){
        avgHueX = (sumHueX[y2][x2])/pixelNum;
        avgHueY = (sumHueY[y2][x2])/pixelNum;
        avgSat = (sumSat[y2][x2])/pixelNum;
        avgLum = (sumLum[y2][x2])/pixelNum;
    }
    double avgHue = atan2(avgHueY, avgHueX) * 180 / PI;
    if(avgHue < 0)
        avgHue += 360;
    return HSLAPixel(avgHue,avgSat, avgLum);
}

vector<int> stats::buildHist(pair<int,int> ul, pair<int,int> lr){

/* your code here */
    int x1 = ul.first;
    int y1 = ul.second;
    int x2 = lr.first;
    int y2 = lr.second;
    vector<int> histogram(36);

    if(x2 >= x1 && y2 >= y1){
      return buildHistHelper(ul, lr);
    }
    else if(x2 >= x1 && y2 < y1){
      pair<int, int> lr2(x2, height - 1);
      pair<int, int> ul2(x1, 0);
      vector<int> upperVect = buildHistHelper(ul2, lr);
      vector<int> lowerVect = buildHistHelper(ul, lr2);
      for(int i = 0; i < 36; i++){
        histogram[i] = upperVect[i] + lowerVect[i];
      }
      return histogram;
    }
    else if(x2 < x1 && y2 >= y1){
      pair<int, int> lr2(width - 1, y2);
      pair<int, int> ul2(0, y1);
      vector<int> upperVect = buildHistHelper(ul2, lr);
      vector<int> lowerVect = buildHistHelper(ul, lr2);
      for(int i = 0; i < 36; i++){
        histogram[i] = upperVect[i] + lowerVect[i];
      }
      return histogram;
    }
    else if(x2 < x1 && y2 < y1) {
      pair<int, int> ul1(0,0);
      pair<int, int> ul2(0, y1);
      pair<int, int> lr2(x2, height - 1);
      pair<int, int> ul3(x1, 0);
      pair<int, int> lr3(width - 1, y2);
      pair<int, int> lr4(width - 1, height - 1);
      vector<int> upperLeftVect = buildHistHelper(ul1, lr);
      vector<int> lowerLeftVect = buildHistHelper(ul2, lr2);
      vector<int> upperRightVect = buildHistHelper(ul3, lr3);
      vector<int> lowerRightVect = buildHistHelper(ul, lr4);
      for(int i = 0; i < 36; i++){
        histogram[i] = upperLeftVect[i] + lowerLeftVect[i] + upperRightVect[i] + lowerRightVect[i];
      }
      return histogram;
    }
	return histogram;
}

vector<int> stats::buildHistHelper(pair<int,int> ul, pair<int,int> lr){
  int x1 = ul.first;
  int y1 = ul.second;
  int x2 = lr.first;
  int y2 = lr.second;
  vector<int> histogram(36);
  for(int k = 0; k < 36; k++){
      if(x1 > 0 && y1 > 0)
          histogram[k] = hist[y2][x2][k] - hist[y2][x1-1][k] - hist[y1-1][x2][k] + hist[y1-1][x1-1][k];
      else if(x1 <= 0 && y1 > 0)
          histogram[k] = hist[y2][x2][k] - hist[y1-1][x2][k];
      else if(x1 > 0 && y1 <= 0)
          histogram[k] = hist[y2][x2][k] - hist[y2][x1-1][k];
      else if(x1 <= 0 && y1 <= 0)
          histogram[k] = hist[y2][x2][k];
  }
  return histogram;
}

// takes a distribution and returns entropy
// partially implemented so as to avoid rounding issues.
double stats::entropy(vector<int> & distn,int area){

    double entropy = 0.;

/* your code here */

    for (int i = 0; i < 36; i++) {
        if (distn[i] > 0 )
            entropy += ((double) distn[i]/(double) area)
                                    * log2((double) distn[i]/(double) area);
    }

    return  -1 * entropy;

}

double stats::entropy(pair<int,int> ul, pair<int,int> lr){

/* your code here */
    int area = (int) rectArea(ul,lr);
    vector<int> distribution = buildHist(ul,lr);
    return entropy(distribution, area);
}
