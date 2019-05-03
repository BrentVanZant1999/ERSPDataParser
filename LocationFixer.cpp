 // LocationFixer.cpp : This file contains the 'main' function. Program execution begins and ends there.
// ERSP Early Research Scholars Program
// Brent VanZant 
// Functionality: Take in an input csv with sporadic timestamps on data collection and to organize it into precise one minute intervals. 
#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <map> 
#include <math.h>
#define M_PI 3.14159265358979323846
#include <cmath> 



using namespace std;

//file constants
int NUM_FILES = 8; 
string FILENAMES[8] = { "2019-03-06", "2019-03-10", "2019-03-11", "2019-03-22", "2019-03-30", "2019-04-01", "2019-04-08", "2019-04-09" };
//splicing data constants
int NUM_COLS = 3; 
int CUT_LOCATION = 11;
int CUT_LOCATION_HR_INIT = 0; 
int CUT_LENGTH = 2; 
int CUT_LOCATION_MN_INIT = 3; 
int HOURS_IN_DAY = 24; 
double earthRadiusKm = 6371.0;
double LAT_SUPERCOMPUTER_SENSOR = 32.88432; 
double LONG_SUPERCOMPUTER_SENSOR = -117.23977;
//time editing constants
int CONVERSION_INCREMENT = 1020;
int MINS_IN_DAY = 1440;


// This function converts decimal degrees to radians
double deg2rad(double deg) {
	return (deg * M_PI / 180);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) {
	return (rad * 180 / M_PI);
}

/**
 * Returns the distance between two points on the Earth.
 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
 * @param lat1d Latitude of the first point in degrees
 * @param lon1d Longitude of the first point in degrees
 * @param lat2d Latitude of the second point in degrees
 * @param lon2d Longitude of the second point in degrees
 * @return The distance between the two points in kilometers
 */
double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
	double lat1r, lon1r, lat2r, lon2r, u, v;
	lat1r = deg2rad(lat1d);
	lon1r = deg2rad(lon1d);
	lat2r = deg2rad(lat2d);
	lon2r = deg2rad(lon2d);
	u = sin((lat2r - lat1r) / 2);
	v = sin((lon2r - lon1r) / 2);
	return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

int main() {
	for (int q = 0; q < NUM_FILES; q++) {
		string::size_type sz;
		string currentFileName = FILENAMES[q] + ".csv";
		ifstream data(currentFileName);


		string curLine;

		vector<string> collectedLong;
		vector<string> collectedLat;


		map<int, string> longMap;
		map<int, string> latMap;
		map<int, string> longDoubleMap;
		map<int, string> latDoubleMap;

		map<int, double> longValueMap;
		map<int, double> latValueMap;
		map<int, double> distFrom;

		map<int, bool> dataPointMap;


		double stationaryLat = 0.0f;
		double stationaryLong = 0.0f;


		vector<string> collectedZuluTimes; //store zulu times
		vector<string> collectedPDTTimes;  //hold PDT times 

		//handle if header has been read
		bool header = true;

		//handle file reading and array generation 
		while (getline(data, curLine))
		{
			//handle reading 
			if (header == false) {
				stringstream lineStream(curLine);

				//setup a string to hold individual cell data 
				string cell;
				int colReadCount = 0;
				//handle getting data
				while (getline(lineStream, cell, ','))
				{
					if (colReadCount == 0) {
						collectedZuluTimes.push_back(cell); //push back time 
					}
					else if (colReadCount == 1) {
						collectedLat.push_back(cell); //push back latitude
					}
					else if (colReadCount == 2) {
						collectedLong.push_back(cell);  //push back longitude
					}
					colReadCount++;
				}
			}
			else {
				header = false;
			}
		}
		data.close();

		int startTime = 0;
		int endTime = 0;
		//handle data string conversion and organizing
		for (int i = 0; i < collectedZuluTimes.size(); i++) {
			cout << "I VAL:" << i << endl;
			//updating zulu time string 
			cout << collectedZuluTimes[i] << endl;
			string toEdit = collectedZuluTimes[i];
			string newString = toEdit.substr(CUT_LOCATION);
			newString.pop_back();
			cout << newString << endl;
			//collectedZuluTimes[i] = newString;

			//get hour and min val of zulu time string. 
			string currHr = newString.substr(CUT_LOCATION_HR_INIT, CUT_LENGTH);
			string currMin = newString.substr(CUT_LOCATION_MN_INIT, CUT_LENGTH);
			int currHrVal = stoi(currHr, &sz);
			int currMinVal = stoi(currMin, &sz);
			int minVal = (currMinVal + (60 * currHrVal));
			cout << "ZULU TIME:" << minVal << endl;
			int newMinVal = (minVal + CONVERSION_INCREMENT) % MINS_IN_DAY;
			cout << "PDT TIME:" << newMinVal << endl;

			if (i == 0) {
				startTime = newMinVal;
			}
			if (i == collectedZuluTimes.size() - 1) {
				endTime = newMinVal;
			}
			dataPointMap.insert(pair<int, bool>(newMinVal, true));
			longMap.insert(pair<int, string>(newMinVal, collectedLong[i]));
			latMap.insert(pair<int, string>(newMinVal, collectedLat[i]));
		}
		bool hasFoundFirst = false;
		for (int n = 0; n < MINS_IN_DAY; n++) {
			auto searchData = dataPointMap.find(n);
			if (searchData == dataPointMap.end()) {
				dataPointMap.insert(pair<int, bool>(n, false));
			}
			else {
				string number;
				double lat_number;
				double long_number;
				number = latMap.find(n)->second;
				lat_number = atof(number.c_str());
				cout << "string: " << number << endl;
				cout << "double: " << setprecision(8) << lat_number << endl;

				latValueMap.insert(pair<int, double>(n, lat_number));

				number = longMap.find(n)->second;
				long_number = atof(number.c_str());
				longValueMap.insert(pair<int, double>(n, long_number));
				if (hasFoundFirst == false) {
					hasFoundFirst = true;
					stationaryLat = lat_number;
					stationaryLong = long_number;
				}


			}
		}

		std::ofstream outputFile;
		outputFile.open("UPDATED_" + currentFileName + ".csv");
		outputFile << "time" << "," << "latitude" << "," << "longitude" << "," << "isDataPoint" << "," << "within15(M)"  << "," << "distanceFromSensor(M)" << "," << "distanceFromSupercomputerSensor(M)" << endl;
		string lastLatString = "";
		string lastLongString = "";
		bool hasLast = false;
		int lastValIndex = 0;
		for (int z = startTime; z <= endTime; z++) {
			int hourToDisplay = z / 60;
			string display;
			display = currentFileName;
			display.pop_back();
			display.pop_back();
			display.pop_back();
			display.pop_back();
			display += "T";
			if (hourToDisplay < 10) {
				display += "0" + to_string(hourToDisplay);
			}
			else {
				display += to_string(hourToDisplay);
			}
			int minToDisplay = z % 60;
			string minDisplay;
			display += ":";
			if (minToDisplay < 10) {
				display += "0" + to_string(minToDisplay);
			}
			else {
				display += to_string(minToDisplay);
			}
			string newLongVal = "NULL";
			string newLatVal = "NULL";
			double kmDistToStationarySensor;
			double kmDistToUCSDSensor;
			bool is15M = false; 
			if (dataPointMap.find(z)->second == false) {
				if (hasLast) {
					kmDistToStationarySensor = distanceEarth(stationaryLat, stationaryLong, latValueMap.find(lastValIndex)->second, longValueMap.find(lastValIndex)->second);
					kmDistToStationarySensor = kmDistToStationarySensor * 1000;
					kmDistToUCSDSensor = distanceEarth(LAT_SUPERCOMPUTER_SENSOR, LONG_SUPERCOMPUTER_SENSOR, latValueMap.find(lastValIndex)->second, longValueMap.find(lastValIndex)->second);
					kmDistToUCSDSensor = kmDistToUCSDSensor * 1000;
					if (kmDistToStationarySensor < 15)
					{
						is15M = true; 
					}
					else {
						is15M = false;
					}
					outputFile << display << "," << lastLatString << "," << lastLongString << "," << "FALSE" << ","  << is15M << "," << setprecision(8) << kmDistToStationarySensor << "," << setprecision(8) << kmDistToUCSDSensor << endl;
				}
				else {
					outputFile << display << "," << newLatVal << "," << newLongVal << "," << "FALSE" << endl;
				}
			}
			else {
				hasLast = true;
				lastValIndex = z;
				lastLatString = latMap.find(z)->second;
				lastLongString = longMap.find(z)->second;
				kmDistToStationarySensor = distanceEarth(stationaryLat, stationaryLong, latValueMap.find(z)->second, longValueMap.find(z)->second);
				kmDistToStationarySensor = kmDistToStationarySensor * 1000;
				kmDistToUCSDSensor = distanceEarth(LAT_SUPERCOMPUTER_SENSOR, LONG_SUPERCOMPUTER_SENSOR, latValueMap.find(z)->second, longValueMap.find(z)->second);
				kmDistToUCSDSensor = kmDistToUCSDSensor * 1000;
				if (kmDistToStationarySensor < 15)
				{
					is15M = true;
				}
				else {
					is15M = false;
				}
				outputFile << display << "," << lastLatString << "," << lastLongString << "," << "TRUE" << "," << is15M << "," << setprecision(8) << kmDistToStationarySensor << "," << setprecision(8) << kmDistToUCSDSensor << endl;
			}

		}
		outputFile.close();
	}
}
		
