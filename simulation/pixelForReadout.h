
#ifndef TP_PIXELFORREADOUT_H

#define TP_PIXELFORREADOUT_H

// Define the pixel structure for raw data readout.

struct pixel {
  int col;
  int row;
  int ana;
  float anaVcal;
  int roc;
  int colROC;
  int rowROC;
  int raw[6];
  //DP float xy[2];
};

struct cluster {
  vector <pixel> vpix;
  int size;
  int sumA;//DP
  float charge;
  float col,row;
  int layer;
  //DP double xy[2]; // local coordinates
  //DP double xyz[3];
};

#endif // TP_PIXELFORREADOUT_H
