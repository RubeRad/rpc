#include "rpc_ocv.h"

#include <opencv2/core/core_c.h>
#include <opencv2/calib3d/calib3d_c.h>
#include <opencv2/calib3d/calib3d.hpp>

using namespace RPC_NS;

double // return fit RMS>0 (pixels) if OK, else <0
fit_dlt(RPC<double>& rpc,
        double* coeffs11)
{
   std::vector<cv::Point3d> gps;
   std::vector<cv::Point2d> ips;
   // assume same ground/image scaling as rpc, so generate 3x3x3 fit points in
   // +/-1 range
   for (int ix=-1; ix<=1; ix++) {
   for (int iy=-1; iy<=1; iy++) {
   for (int iz=-1; iz<=1; iz++) {
      double nx(ix), ny(iy), nz(iz), ns(0), nl(0);
      xyz2xy(rpc.coeffs, nx,ny,nz, ns,nl);
      gps.emplace_back(nx,ny,nz);
      ips.emplace_back(ns,nl);
   }}}

   cv::Mat fpp(3,3, CV_64F);
   for (int i=0; i<3; ++i) fpp.at<double>(i,i) = 1.0;

   cv::Mat rvec, tvec, no_distortion;
   bool ok = cv::solvePnP(gps, ips, fpp, no_distortion, rvec, tvec);

   if (!ok)
      return -1.0;

   cv::Mat R;
   cv::Rodrigues(rvec, R);
   cv::Mat P;
   cv::hconcat(R, tvec, P);
   P *= 1.0/P.at<double>(2,0);

   cout << P << endl;

   for (size_t i=0; i<gps.size(); ++i) {
      cv::Mat gp(4,1, CV_64F);
      gp.at<double>(0,0) = gps[i].x;
      gp.at<double>(1,0) = gps[i].y;
      gp.at<double>(2,0) = gps[i].z;
      gp.at<double>(3,0) = 1.0;

      cv::Mat iph = P * gp;
      iph *= 1.0/iph.at<double>(2,0);
      cv::Point2d ip(iph.at<double>(0,0), iph.at<double>(1,0));
      cv::Point2d dip = ip-ips[i];

      cout << dip << endl;

   }

   
   return (ok ? 1.0 : -1.0);
}
