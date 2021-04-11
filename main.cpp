#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <cmath>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

VectorXd calculateError(VectorXd state,double gt_x,double gt_y,double gt_vx,double gt_vy){
  VectorXd gt(4);
  gt << gt_x,gt_y,gt_vx,gt_vy;
  VectorXd error = state-gt; 
  return error;
}

MatrixXd calculateJacobian(VectorXd state){
   double px = state(0);
   double py = state(1);
   double vx = state(2);
   double vy = state(3);
   MatrixXd jacobian(3,4);
   double sqrt_power = sqrt((pow(px,2)+pow(py,2)));
   double constant = (px*px + py*py)*sqrt_power;
   jacobian << (px/sqrt_power),(py/sqrt_power),0,0,
               (-1*py/(pow(px,2)+pow(py,2))),(px/(pow(px,2)+pow(py,2))),0,0,
               py*(vx*py - vy*px)/constant, px*(px*vy - py*vx)/constant, px/sqrt_power, py/sqrt_power;
   return jacobian;
}

VectorXd radar2State(VectorXd m_radar){
   double r_m_rho = m_radar(0);
   double r_m_pi = m_radar(1);
   double r_m_rdot = m_radar(2);
   double r_m_px = r_m_rho*cos(r_m_pi);
        if ( r_m_px < 0.0001 ) {
        r_m_px = 0.0001;
      }
        double r_m_py = r_m_rho*sin(r_m_pi);
        if(r_m_py< 0.0001){
          r_m_py = 0.0001;
        }  
        double r_m_vx = r_m_rdot*cos(r_m_pi);
        double r_m_vy = r_m_rdot*sin(r_m_pi);
  VectorXd r_x(4) ;
  r_x << r_m_px,r_m_py,r_m_vx,r_m_vy;
  return r_x;
}
 
int main()
{ 
  ifstream reader("C:\\DEV\\website\\data\\obj_pose-laser-radar-synthetic-input.txt");
  int tick = -10;
  //Matrix for laser
   MatrixXd R_laser_ = MatrixXd(2, 2);//measurement noise
    //measurement covariance matrix - laser
   R_laser_ << 0.0225, 0,
        0, 0.0225;

    MatrixXd l_H(2,4);// used to convert predicted state to measurement state, shape wise
     l_H << 1,0,0,0,
          0,1,0,0;

    VectorXd l_m_x(2);//measured state

  //Matrix for Radar
   MatrixXd R_radar_ = MatrixXd(3, 3);// measurement noise
    //measurement covariance matrix - radar
   R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
  
   VectorXd R_m_x(3);
   MatrixXd Hj(3,4);


  string line;
  string part;

  //Global Variables , same for laser and radar
  MatrixXd P(4,4);//process covariance matrix
  MatrixXd p_P(4,4);//predicted covariance matrix
  MatrixXd F(4,4);//prediction matrix
  F << 1,0,0,0,
         0,1,0,0,
         0,0,1,0,
         0,0,0,1; 
  MatrixXd Q(4,4);//process covariance matrix
  VectorXd y;// difference between predicted and measured
  MatrixXd I = MatrixXd::Identity(4,4);
  VectorXd x(4);//state
  VectorXd p_x(4);//predicted state
  VectorXd noise(4);//noise due to model simplicity
  MatrixXd S;
  MatrixXd K;
  VectorXd p_R_x(3);//state measurement in radar space

  char type;
  int count = 0;//for initializing
  double gt_x;//ground truth
  double gt_y;//ground truth
  double gt_vy;//ground truth
  double gt_vx;//ground truth

  double noise_ax = 9.0;//noise due to model simplicity
  double noise_ay = 9.0;//noise due to model simplicity
 
  long long c_timestamp;
  long long p_timestamp;
  double dt;

  float m_px,m_py,m_vx,m_vy;
  float m_rho,m_pi,m_rdot;

  while(getline(reader,line)){
    istringstream my_stream(line);
    my_stream>>type;
    if(type=='L'){
      my_stream>>m_px;
      my_stream>>m_py;
      l_m_x<< m_px,m_py;
      my_stream>>c_timestamp;
      //for intializing state and covariance matrix
      if(count==0){
        x <<m_px,m_py,0,0;
        P << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;
        p_timestamp = c_timestamp;
        count++;
        continue;
      }
      dt = (c_timestamp - p_timestamp)/1000000.0;
       p_timestamp = c_timestamp;
      F(0,2) = dt;
      F(1,3) = dt;
      noise<< (noise_ax*dt*dt/2),(noise_ay*dt*dt/2),(noise_ax*dt),(noise_ay*dt);
      p_x = F*x + noise;
      
      Q << (dt*dt*dt*dt*noise_ax*noise_ax/4),0,(dt*dt*dt*noise_ax*noise_ax/2),0,
            0,(dt*dt*dt*dt*noise_ay*noise_ay/4),0,(dt*dt*dt*noise_ay*noise_ay/2),
            (dt*dt*dt*noise_ax*noise_ax/2),0,(dt*dt*noise_ax*noise_ax),0,
            0,(dt*dt*dt*noise_ay*noise_ay/2),0,(dt*dt*noise_ay*noise_ay);
      p_P = F*P*F.transpose() + Q;
      y = l_m_x-(l_H*p_x);
      S = l_H*p_P*l_H.transpose() + R_laser_;
      K = p_P*l_H.transpose()*S.inverse();
      x = p_x + K*y;
      P = (I-K*l_H)*p_P;
      my_stream>>gt_x;
      my_stream>>gt_y;
      my_stream>>gt_vx;
      my_stream>>gt_vy;
      cout<<"Error"<<calculateError(x,gt_x,gt_y,gt_vx,gt_vy)<<endl;
       VectorXd gt(4);
       gt << gt_x,gt_y,gt_vx,gt_vy;

      cout<<"predicted state = "<<x<<endl;
      cout<<"ground truth"<<gt<<endl;
    }
    if(type=='R'){
     my_stream>>m_rho;
      my_stream>>m_pi;
      my_stream>>m_rdot;
      R_m_x << m_rho,m_pi,m_rdot;
      my_stream>>c_timestamp;
       if(count==0){
        x  = radar2State(R_m_x);
        P << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;
        p_timestamp = c_timestamp;
        count++;
        continue;
      }
      dt = (c_timestamp - p_timestamp)/1000000.0;
      p_timestamp = c_timestamp;
      F(0,2) = dt;
      F(1,3) = dt;
      noise<< (noise_ax*dt*dt/2),(noise_ay*dt*dt/2),(noise_ax*dt),(noise_ay*dt);
      p_x = F*x + noise;
      Q << (dt*dt*dt*dt*noise_ax*noise_ax/4),0,(dt*dt*dt*noise_ax*noise_ax/2),0,
            0,(dt*dt*dt*dt*noise_ay*noise_ay/4),0,(dt*dt*dt*noise_ay*noise_ay/2),
            (dt*dt*dt*noise_ax*noise_ax/2),0,(dt*dt*noise_ax*noise_ax),0,
            0,(dt*dt*dt*noise_ay*noise_ay/2),0,(dt*dt*noise_ay*noise_ay);
      p_P = F*P*F.transpose() + Q;
      double px = p_x(0);
      double py = p_x(1);
      double vx = p_x(2);
      double vy = p_x(3);
      double rho = sqrt(px*px + py*py);
      double theta = atan2(py, px);
      double rho_dot = (px*vx + py*vy) / rho;
      p_R_x << rho,theta,rho_dot;
      y = R_m_x - p_R_x;
      while ( y(1) > M_PI || y(1) < -M_PI ) {
    if ( y(1) > M_PI ) {
      y(1) -= M_PI;
    } else {
      y(1) += M_PI;
    }
  }
      Hj = calculateJacobian(p_x);
      S = Hj*p_P*Hj.transpose() + R_radar_;
      K = p_P*Hj.transpose()*S.inverse();
      x = p_x + K*y;
      P = (I-K*Hj)*p_P;
      my_stream>>gt_x;
      my_stream>>gt_y;
      my_stream>>gt_vx;
      my_stream>>gt_vy;
       VectorXd gt(4);
       gt << gt_x,gt_y,gt_vx,gt_vy;
      cout<<"Error"<<calculateError(x,gt_x,gt_y,gt_vx,gt_vy)<<endl;
      cout<<"predicted state = "<<x<<endl;
      cout<<"ground truth"<<gt<<endl;
    } 
  }
  return 0;
}
