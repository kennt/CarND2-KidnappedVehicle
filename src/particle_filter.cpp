/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <limits>

#include "map.h"
#include "particle_filter.h"


using namespace std;

// distance  Returns the squared distance between two points
static double distance(double x, double y, double x2, double y2)
{
  return (x-x2)*(x-x2) + (y-y2)*(y-y2);
}


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position
  //   (based on estimates of x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  // CHECK: this doesn't compile on my Mac
  //std::random_device rd;
  //std::mt19937 gen(rd);
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (unsigned int i=0; i<num_particles; i++)
  {
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
  std::default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  for (auto& particle : particles)
  {
    double theta0 = particle.theta;
    if (fabs(yaw_rate) < 0.000001)
    {
      particle.x += velocity*delta_t*cos(theta0);
      particle.y += velocity*delta_t*sin(theta0);
      // particle.theta doesn't change
    }
    else
    {
      double thetaf = theta0 + yaw_rate*delta_t;
      particle.x += velocity*(sin(thetaf) - sin(theta0))/yaw_rate;
      particle.y += velocity*(cos(theta0) - cos(thetaf))/yaw_rate;
      particle.theta = thetaf;
    }

    // CHECK: should we add the noise before or after we calculate x,y
    // Add in random noise
    particle.x += dist_x(gen);
    particle.y += dist_y(gen);
    particle.theta += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

    // For now, do this the naive way (i.e. brute force)
    // Maybe add a KD-tree for searching

}

static double multi_gaussian(double x, double y, double x2, double y2, double std_landmark[])
{
  // CHECK : divide by zero?
  double expx = (x - x2)*(x - x2) / (2*std_landmark[0]*std_landmark[0]);
  double expy = (y - y2)*(y - y2) / (2*std_landmark[1]*std_landmark[1]);
  return exp(-(expx + expy)) / (2*M_PI*std_landmark[0]*std_landmark[1]);
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
        const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html

  for (auto& particle : particles)
    {
      // The observations are in the CAR coordinate system
      // Transform them into MAP coordinates with respect to the particle
      // In other words, assume that they have been observed from the particle
      // So, where does the particle think the observations are in MAP coordinates

      // Initialize the particle for the new observations
      double total_weight = 1.0;
      particle.associations.clear();
      particle.associations.reserve(observations.size());
      particle.sense_x.clear();
      particle.sense_x.reserve(observations.size());
      particle.sense_y.clear();
      particle.sense_y.reserve(observations.size());

      for(const auto& observation: observations)
      {
        double rotation = particle.theta;
        double cosr = cos(rotation);
        double sinr = sin(rotation);
        // observation coords in the CAR (local) coordinate system
        double xc = observation.x;
        double yc = observation.y;
        // final calulation, transfrom into MAP (global) coordinates
        double xm = particle.x + cosr*xc - sinr*yc;
        double ym = particle.y + sinr*xc + cosr*yc;

        // find the nearest landmark to xm, ym
        // TODO: cam we narrow down the sensor range?
        // Maybe reuse the previous distance + distance traveled ?
        int id = map_landmarks.find_nearest_neighbor(xm, ym, sensor_range);

        // TODO: what to do if we can't find a nearest neighbor
        // (i.e. if all landmarks are outside of sensor range)

        // Update the associations
        particle.associations.push_back(id);
        particle.sense_x.push_back(xm);
        particle.sense_y.push_back(ym);

        // Adjust the IDs since the landmark id starts at 1

        // TODO: this assumes that the landmark ids start at 1, thus the offset
        // Maybe this should use a map of id's to array positions to accomdate
        // arbitrary ids.
        total_weight *= multi_gaussian(xm, ym,
                                       map_landmarks.landmark_list[id-1].x_f,
                                       map_landmarks.landmark_list[id-1].y_f,
                                       std_landmark);
      }

      particle.weight = total_weight;
    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // Gather all the weights together
  temp_weights.clear();
  temp_weights.reserve(num_particles);
  for (const auto& particle : particles)
    temp_weights.push_back(particle.weight);

  std::default_random_engine gen;
  discrete_distribution<int> distribution(temp_weights.begin(), temp_weights.end());

  // Loop over the distribution gathering the particles that were picked
  temp_particles.clear();
  temp_particles.reserve(num_particles);
  for (unsigned int i=0; i<particles.size(); i++)
  {
    int id = distribution(gen);
    temp_particles.push_back(particles[id]);
  }
  particles.swap(temp_particles);
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
  return particle;
}


string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
    vector<double> v = best.sense_x;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
    vector<double> v = best.sense_y;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

int Map::find_nearest_neighbor(double x, double y, double sensor_range) const
{
  // Brute-force search, however, ignore measurements that are
  // outside of the sensor range
  float minimum_distance = std::numeric_limits<double>::max();
  int landmark_id = -1;

  for (const auto& landmark : landmark_list)
  {
    if ((fabs(landmark.x_f - x) > sensor_range) ||
        (fabs(landmark.y_f - y) > sensor_range))
      continue;
    float current_distance = distance(x, y, landmark.x_f, landmark.y_f);
    if (current_distance < minimum_distance)
    {
      landmark_id = landmark.id_i;
      minimum_distance = current_distance;
    }
  }
  return landmark_id;
}

