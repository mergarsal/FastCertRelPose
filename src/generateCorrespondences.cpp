#include "generateCorrespondences.h"

#define PI 3.14159

namespace Essential
{
 /* from random_generators.cpp */

 Vector3 generateRandomPointTruncated(double FoV, double min_depth, double max_depth)
        {
                generateRandomPointTruncated(FoV, min_depth, max_depth, false, 20);
        }

Vector3 generateRandomPointNaive(double min_depth, double max_depth)
{
  /* Naive points: no Fov */
  Eigen::Vector3d cleanPoint;
  cleanPoint[0] = (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0;
  cleanPoint[1] = (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0;
  cleanPoint[2] = (((double) rand()) / ((double) RAND_MAX));  // always non-negative Z
  Eigen::Vector3d direction = cleanPoint.normalized();
  cleanPoint =
      (max_depth-min_depth) * cleanPoint + min_depth * direction;


  return cleanPoint;

}


Vector3 generateRandomPointTruncated(double FoV, double min_depth, double max_depth, bool allow_coplanar, double max_X)
{

  Vector3 point3D;
  // Note: tan(X) applies mod(X, pi) before computing the actual tan

  // depth
  point3D[2] = min_depth + (max_depth - min_depth) * ( ((double) std::rand() / (double) RAND_MAX));

  double xmax = (tan(FoV * 0.5 * PI / 180)) * point3D[2];
  if (xmax <= 0) xmax *= -1;  

  if (!allow_coplanar)
    // update xmax
    xmax = std::min(xmax, max_X);

  point3D[0] = xmax * (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0;
  point3D[1] = xmax * (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0;

  return point3D;

}

/* From opengv */ 
Vector3 addNoise( double noiseLevel, Vector3 cleanPoint, double focal_length )
{
  //compute a vector in the normal plane (based on good conditioning)
  Eigen::Vector3d normalVector1;
  if(
      (fabs(cleanPoint[0]) > fabs(cleanPoint[1])) &&
      (fabs(cleanPoint[0]) > fabs(cleanPoint[2])) )
  {
    normalVector1[1] = 1.0;
    normalVector1[2] = 0.0;
    normalVector1[0] = -cleanPoint[1]/cleanPoint[0];
  }
  else
  {
    if(
        (fabs(cleanPoint[1]) > fabs(cleanPoint[0])) &&
        (fabs(cleanPoint[1]) > fabs(cleanPoint[2])) )
    {
      normalVector1[2] = 1.0;
      normalVector1[0] = 0.0;
      normalVector1[1] = -cleanPoint[2]/cleanPoint[1];
    }
    else
    {
      normalVector1[0] = 1.0;
      normalVector1[1] = 0.0;
      normalVector1[2] = -cleanPoint[0]/cleanPoint[2];
    }
  }

  normalVector1 = normalVector1 / normalVector1.norm();
  Eigen::Vector3d normalVector2 = cleanPoint.cross(normalVector1);
  double noiseX =
      noiseLevel * (((double) std::rand())/ ((double) RAND_MAX)-0.5)*2.0 / 1.4142;
  double noiseY =
      noiseLevel * (((double) std::rand())/ ((double) RAND_MAX)-0.5)*2.0 / 1.4142;

  Eigen::Vector3d noisyPoint =
      focal_length * cleanPoint + noiseX *normalVector1 + noiseY * normalVector2;
  noisyPoint = noisyPoint / noisyPoint.norm();
  return noisyPoint;

}

Vector3 generateRandomTranslation( double parallax )
{

  Vector3 translation;
  translation[0] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  translation[1] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  translation[2] = -(((double) rand())/ ((double) RAND_MAX));


  return (parallax * translation.normalized());
}


/* From opengv */ 
Matrix3 generateRandomRotation( double maxAngle )
{

  Vector3 rpy;
  rpy[0] = ((double) std::rand())/ ((double) RAND_MAX);
  rpy[1] = ((double) std::rand())/ ((double) RAND_MAX);
  rpy[2] = ((double) std::rand())/ ((double) RAND_MAX);

  rpy[0] = maxAngle*2.0*(rpy[0]-0.5);
  rpy[1] = maxAngle*2.0*(rpy[1]-0.5);
  rpy[2] = maxAngle*2.0*(rpy[2]-0.5);

  Eigen::Matrix3d R1;
  R1(0,0) = 1.0;
  R1(0,1) = 0.0;
  R1(0,2) = 0.0;
  R1(1,0) = 0.0;
  R1(1,1) = cos(rpy[0]);
  R1(1,2) = -sin(rpy[0]);
  R1(2,0) = 0.0;
  R1(2,1) = -R1(1,2);
  R1(2,2) = R1(1,1);

  Eigen::Matrix3d R2;
  R2(0,0) = cos(rpy[1]);
  R2(0,1) = 0.0;
  R2(0,2) = sin(rpy[1]);
  R2(1,0) = 0.0;
  R2(1,1) = 1.0;
  R2(1,2) = 0.0;
  R2(2,0) = -R2(0,2);
  R2(2,1) = 0.0;
  R2(2,2) = R2(0,0);

  Eigen::Matrix3d R3;
  R3(0,0) = cos(rpy[2]);
  R3(0,1) = -sin(rpy[2]);
  R3(0,2) = 0.0;
  R3(1,0) =-R3(0,1);
  R3(1,1) = R3(0,0);
  R3(1,2) = 0.0;
  R3(2,0) = 0.0;
  R3(2,1) = 0.0;
  R3(2,2) = 1.0;

  Matrix3 rotation = R3 * R2 * R1;

  rotation.col(0) = rotation.col(0) / rotation.col(0).norm();
  rotation.col(2) = rotation.col(0).cross(rotation.col(1));
  rotation.col(2) = rotation.col(2) / rotation.col(2).norm();
  rotation.col(1) = rotation.col(2).cross(rotation.col(0));
  rotation.col(1) = rotation.col(1) / rotation.col(1).norm();

  return rotation;
}




void createSyntheticExperiment(
    const size_t n_points,
    double noise,
    double outlier_fraction,
    double FoV,
    double parallax,
    double min_depth,
    double max_depth,
    Vector3 & translation,
    Matrix3 & rotation,
    bearing_vectors_t & points_correspondences,
    Eigen::MatrixXd & gt,
    std::vector<int> & indices_outliers )
    {
        createSyntheticExperiment(n_points, noise, outlier_fraction, FoV,
        parallax, min_depth, max_depth, translation,
        rotation, points_correspondences, gt, indices_outliers, false, 20, 1000000000);

    }

    void createSyntheticExperiment(
    const size_t n_points,
    double noise,
    double outlier_fraction,
    double FoV,
    double parallax,
    double min_depth,
    double max_depth,
    Vector3 & translation,
    Matrix3 & rotation,
    bearing_vectors_t & points_correspondences,
    Eigen::MatrixXd & gt,
    std::vector<int> & indices_outliers,
    double min_epipolar_error_sq   // min epipolar error for the outliers
    )
    {
        createSyntheticExperiment(n_points, noise, outlier_fraction, FoV,
        parallax, min_depth, max_depth, translation,
        rotation, points_correspondences, gt, indices_outliers, false, 20, min_epipolar_error_sq);

    }


void createSyntheticExperiment(
    const size_t n_points,
    double noise,
    double outlier_fraction,
    double FoV,
    double parallax,
    double min_depth,
    double max_depth,
    Vector3 & translation,
    Matrix3 & rotation,
    bearing_vectors_t & points_correspondences,
    Eigen::MatrixXd & gt,
    std::vector<int> & indices_outliers,
    bool allow_coplanar,   // allow to generate synthetic scenes with coplanar points
    double max_X  // this value is the absolute value. Max value for X-axis allowed
    )
    {
    // call the other function
    // with default params
    createSyntheticExperiment(n_points, noise, outlier_fraction, FoV,
        parallax, min_depth, max_depth, translation,
        rotation, points_correspondences, gt, indices_outliers,
        allow_coplanar, max_X, 1000000000);
    }


void createSyntheticExperiment(
    const size_t n_points,
    double noise,
    double outlier_fraction,
    double FoV,
    double parallax,
    double min_depth,
    double max_depth,
    Vector3 & translation,
    Matrix3 & rotation,
    bearing_vectors_t & points_correspondences,
    Eigen::MatrixXd & gt,
    std::vector<int> & indices_outliers,
    bool allow_coplanar,   // allow to generate synthetic scenes with coplanar points
    double max_X,  // this value is the absolute value. Max value for X-axis allowed
    double min_epipolar_error_sq   // min epipolar error for the outliers
    )
{
        createSyntheticExperiment(n_points, noise, outlier_fraction, FoV,
        parallax, min_depth, max_depth, translation,
        rotation, points_correspondences, gt, indices_outliers,
        allow_coplanar, max_X, min_epipolar_error_sq, 800);


}

void createSyntheticExperiment(
    const size_t n_points,
    double noise,
    double outlier_fraction,
    double FoV,
    double parallax,
    double min_depth,
    double max_depth,
    Vector3 & translation,
    Matrix3 & rotation,
    bearing_vectors_t & points_correspondences,
    Eigen::MatrixXd & gt,
    std::vector<int> & indices_outliers,
    bool allow_coplanar,   // allow to generate synthetic scenes with coplanar points
    double max_X,  // this value is the absolute value. Max value for X-axis allowed
    double min_epipolar_error_sq,   // min epipolar error for the outliers
    double focal_length    // focal length in pixelss
    )
{


  // generate a random pointcloud
  gt.setZero();
  for( size_t i = 0; i < n_points; i++ )
        gt.col(i) = generateRandomPointTruncated( FoV, min_depth,
                                max_depth, allow_coplanar, max_X );

  
  double max_rotation = 0.5;

  // 2. Create the relative rotation and translation
  // 2.1 Fix first frame
  Matrix3 R1 = Matrix3::Identity();
  Vector3 T1 = Vector3::Zero();

  // 2.2 Fix second frame
  Matrix3 R2 = Matrix3::Zero();
  Vector3 T2 = Vector3::Zero();

  bool valid_position = false;
  Eigen::MatrixXd points3D_frame2 = gt;
  int max_n_iters = 100, n_iters = 0;

  // create random translation
  // note: we always create a backward movement
  T2 = generateRandomTranslation(parallax);
  
  do
  {
    // Generate a random pose for the second camera
    // such that all the points lie inside its FOV
    R2 = generateRandomRotation(max_rotation);

    // compute relative coordinates
    points3D_frame2.setZero();
    for( size_t i = 0; i < n_points; i++ )
    {
      Vector3 temp = R2.transpose() * (gt.col(i)- T2);
      points3D_frame2.col(i) = Eigen::Map<Vector3>(temp.data(), 3, 1);
    }

    
    // check condition for FoV

    Eigen::VectorXd ratio_X = (points3D_frame2.row(0)).cwiseQuotient((points3D_frame2.row(2)));
    double max_ratio_x = (ratio_X.array().abs()).maxCoeff();
    // check if any point was outside the FoV

    if (max_ratio_x > tan(FoV * 0.5 * PI / 180))
    {
        n_iters++;
        continue;
    }


    // std::cout << "X-coordinates of points is valid \n";
    Eigen::VectorXd ratio_Y = (points3D_frame2.row(1)).cwiseQuotient((points3D_frame2.row(2)));
    double max_ratio_y = (ratio_Y.array().abs()).maxCoeff();
    // check if any point was outside the FoV

    if (max_ratio_y > tan(FoV * 0.5 * PI / 180))
    {
        n_iters++;
        continue;
    }


    // the position is valid
    valid_position = true;

  }while(!valid_position && (n_iters <= max_n_iters));

  // check invalid rotation 
  if ((!valid_position) && (n_iters > max_n_iters))
        std::cout << "[ERROR] Pose is not valid.\nSome points do not lie in the FOV\n";
  
  
  rotation = R2;
  translation = T2.normalized();


  // Generate the correspondences
  for( size_t i = 0; i < n_points; i++ )
  {
    Vector3 obs1, obs2;
    obs1 = gt.col(i);
    obs1.normalize();

    obs2 = points3D_frame2.col(i);
    obs2.normalize();

    //add noise
    if( noise > 0.0 )
    {
      obs1 = addNoise(noise, obs1, focal_length);
      obs2 = addNoise(noise, obs2, focal_length);
    }

    CorrespondingFeatures correspondence;
    correspondence.bearing_vector_0 = obs1.normalized();
    correspondence.bearing_vector_1 = obs2.normalized();
    correspondence.weight_ = 1.0;

    // add it to the std::vector
    points_correspondences.push_back(correspondence);
  }

  /*   OUTLIERS   */

  // add outliers
  size_t number_outliers = (size_t) floor(outlier_fraction*n_points);
  size_t max_number_iters = 50;
  for(size_t i = 0; i < number_outliers; i++)
  {
    size_t i_iter = 0;
    bool valid_outlier = false;
    do
    {
                
            Vector3 outlier_correspondence = generateRandomPointTruncated(170, min_depth,
            max_depth, true, max_X);
           

            //normalize the bearing vector
            outlier_correspondence = rotation.transpose() * (outlier_correspondence - T2);
            outlier_correspondence.normalize();

            // check residual
            Matrix3 t_cross;
                t_cross.setZero();
                t_cross(0, 1) = -translation(2); t_cross(0, 2) = translation(1);
                t_cross(1, 0) = translation(2); t_cross(1, 2) = -translation(0);
                t_cross(2, 0) = -translation(1); t_cross(2, 1) = translation(0);

            double res = (pow(
                        ((points_correspondences[i].bearing_vector_0).transpose() * t_cross * rotation * outlier_correspondence),
                         2));


  
            if (res > min_epipolar_error_sq)
            {
                points_correspondences[i].bearing_vector_1 = outlier_correspondence;
                indices_outliers.push_back(i);

                valid_outlier = true;
            }
            
            if ((!valid_outlier) && (i_iter > max_number_iters))
            {
                // break the loop
                outlier_correspondence[0] = 1.0;
                outlier_correspondence[1] = 1.0;
                outlier_correspondence[2] = 1.0;

                outlier_correspondence.normalize();

                points_correspondences[i].bearing_vector_1 = outlier_correspondence;
                indices_outliers.push_back(i);

                valid_outlier = true;
            }

            // increase the counter
            i_iter += 1;


    }while(!valid_outlier);
  }


}



}  // end of essential namespace
