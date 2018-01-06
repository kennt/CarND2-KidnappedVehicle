/*
 * map.h
 *
 *  Created on: Dec 12, 2016
 *      Author: mufferm
 */

#ifndef MAP_H_
#define MAP_H_

class Map {
public:
	
	struct single_landmark_s{

		int id_i ; // Landmark ID
		float x_f; // Landmark x-position in the map (global coordinates)
		float y_f; // Landmark y-position in the map (global coordinates)
	};


  std::vector<single_landmark_s> landmark_list ; // List of landmarks in the map

	/**
	 * find_nearest_neighbor Finds the landmark closest to the
	 * 	 point (x,y) using the L2-Norm (Euclidean distance)
	 * @param x x-coordinate in global (map) coordinates
	 * @param y y-coordinate in global (map) coordinates
	 * @param sensor_range The limits of the sensor, this (may) be
	 *	 used to exclude points outside of sensor range.
	 * @returns the id of the landmark that is closes to (x,y)
	 * or -1 if the landmark could not be identified.
	 */
	int find_nearest_neighbor(double x, double y, double sensor_range ) const;

};



#endif /* MAP_H_ */
