/*
 *   Copyright (C) 2015 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License
 *   version 3, as published by the Free Software Foundation.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with IsoSpec.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ISOMATH_HPP
#define ISOMATH_HPP

double logFactorial(int n);
double NormalCDFInverse(double p);
double NormalCDFInverse(double p, double mean, double stdev);
double NormalCDF(double x, double mean, double stdev);
double NormalPDF(double x, double mean = 0.0, double stdev = 1.0);



#endif /* ISOMATH_HPP */
