/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::default_random_engine gen;

  num_particles = 100; // TODO: Set the number of particles
  
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  for(int i =0;i<num_particles;i++){
      Particle p;
      p.id = i;
      p.x = dist_x(gen);
      p.y = dist_y(gen);
      p.theta = dist_theta(gen);
      p.weight = 1.0;
      particles.push_back(p);
      //weights.push_back(1);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  
  for(int i =0;i<num_particles;i++){
    //double n_x,n_y,n_theta;
    if(fabs(yaw_rate) < 0.00001){
      //n_x = particles[i].x + (velocity*delta_t)*cos(particles[i].theta);
      particles[i].x += velocity*delta_t*cos(particles[i].theta);
      //n_y = particles[i].y + (velocity*delta_t)*sin(particles[i].theta);
      particles[i].y += velocity*delta_t*sin(particles[i].theta);
	  //n_theta = particles[i].theta;
    }
    else{
   
      particles[i].x += velocity*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta))/yaw_rate;
      particles[i].y += velocity*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t))/yaw_rate;
      particles[i].theta += yaw_rate*delta_t;
    }
    
    std::normal_distribution<double> x_noise(particles[i].x, std_pos[0]);
    std::normal_distribution<double> y_noise(particles[i].y, std_pos[1]);
    std::normal_distribution<double> theta_noise(particles[i].theta, std_pos[2]);
      
    particles[i].x = x_noise(gen);
    particles[i].y = y_noise(gen);
    particles[i].theta = theta_noise(gen);                                                   
                                                   
  }
  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
   for(unsigned int i=0;i<observations.size();i++){
     double min_dist = 1000000;
     double distances = 0;
     int id = 0;
     for(unsigned int j=0;j<predicted.size();j++){
       distances = dist(observations[i].x,observations[i].y, predicted[j].x, predicted[j].y);
       
       if(distances < min_dist){
         id = predicted[j].id;
         min_dist = distances;
       }
     }
     observations[i].id = id;
   }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   //Transforming coordinates to Map coordinates.
  double std_0, std_1;
  
  std_0 = 2.0 * std_landmark[0] * std_landmark[0];
  std_1 = 2.0 * std_landmark[1] * std_landmark[1];
  
  for(int i=0;i<num_particles;i++){
    vector<LandmarkObs> within_sensor_rng, tf_obs;
    for(unsigned int j=0;j<observations.size();j++){
      LandmarkObs obs;
      obs.id = -1;
      //Consulted SDCN Google Hangout Video on youtube to fix this transformation.
      obs.x = particles[i].x + (observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta));
      obs.y =  particles[i].y + (observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta));
      tf_obs.push_back(obs);
  	}
    
    for(unsigned int k=0;k<map_landmarks.landmark_list.size();k++){
      LandmarkObs valid_lm;
      double particleToLmDist;
      particleToLmDist = dist(particles[i].x,particles[i].y,map_landmarks.landmark_list[k].x_f,map_landmarks.landmark_list[k].y_f);
      if(particleToLmDist < sensor_range){
        valid_lm.id = map_landmarks.landmark_list[k].id_i;
        valid_lm.x  = map_landmarks.landmark_list[k].x_f;
        valid_lm.y  = map_landmarks.landmark_list[k].y_f;
        within_sensor_rng.push_back(valid_lm);
      }
    }
    
    dataAssociation(within_sensor_rng,tf_obs);
    
    particles[i].weight = 1.0;
    
    for(unsigned int j=0;j<tf_obs.size();j++){
      double mu_x, mu_y;
      //Udacity knowledge hub suggestion fixed indexing issues.
   	  for(unsigned int k=0;k<within_sensor_rng.size();k++){ 
        if(within_sensor_rng[k].id == tf_obs[j].id){
          mu_x = within_sensor_rng[k].x;
          mu_y = within_sensor_rng[k].y; 
        }
      }
      double x_side = pow((tf_obs[j].x-mu_x),2.0)/std_0;
      double y_side = pow((tf_obs[j].y-mu_y),2.0)/std_1;
      double tmp_wgt = exp(-(x_side + y_side))/(2.0*M_PI*std_landmark[0]*std_landmark[1]);
      if(tmp_wgt>0){
          particles[i].weight *= tmp_wgt;
      }
    }
    weights.push_back(particles[i].weight);

  }  
   
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  vector<Particle> r_particles;
  std::discrete_distribution<int> random_select(weights.begin(),weights.end());
  
  for(int j=0;j<num_particles;++j){
    r_particles.push_back(particles[random_select(gen)]);
  }
  
  particles = r_particles;
  weights.clear();
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
 
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}