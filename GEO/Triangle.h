/*
 * \file Triangle.h
 *
 *  Created on: Mar 23, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include <vector>

// GeoLib
#include "Point.h"

namespace GEOLIB
{
/** \brief Class Triangle consists of a reference to a point vector and
 * a vector that stores the indices in the point vector.
 * A surface is composed by triangles. The class Surface stores the position
 * of pointers to the points of triangles in the m_sfc_pnt_ids vector.
 * */
class Triangle
{
public:
    /**
     * construction of object, initialization of reference to point vector
     */
    Triangle(std::vector<Point*> const& pnt_vec);

    /**
     * construction of object, initialization of reference to point vector,
     * saves the three indices describing a triangle
     */
    Triangle(std::vector<Point*> const& pnt_vec, size_t pnt_a, size_t pnt_b,
             size_t pnt_c);

    /**
     * saves three indices describing a triangle
     * */
    void setTriangle(size_t pnt_a, size_t pnt_b, size_t pnt_c);

    /** \brief const access operator to access the index
     * of the i-th triangle point
     */
    const size_t& operator[](size_t i) const
    {
        assert(i < 3);
        return _pnt_ids[i];
    }

    //	/** \brief access operator to access the index
    //	 * of the i-th triangle point
    //	 * */
    //	size_t& operator[] (size_t i) {
    //		assert (i < 3);
    //		return _pnt_ids[i];
    //	}

    /**
     * \brief const access operator to access the i-th triangle Point
     */
    const Point* getPoint(size_t i) const
    {
        assert(i < 3);
        return _pnts[_pnt_ids[i]];
    }

    /**
     * checks if point is in triangle
     * @param pnt
     * @return true, if point is in triangle, else false
     */
    bool containsPoint(const double* pnt, double eps = 0) const;

    bool containsPoint(const Point& pnt, double eps = 0) const
    {
        return containsPoint(pnt.getData(), eps);
    }
    /**
     * projects the triangle points to the x-y-plane and
     * checks if point pnt is contained into the triangle
     * @param pnt the point to test for
     * @return true, if the point is into the projected triangle
     */
    bool containsPoint2D(const double* pnt) const;
    bool containsPointAreaCheck(const double* pnt, double eps) const; //CMCD Check point in triangle not working for simple evaluation.
    void calculateNormal();
    double const* getNormal() const { return _normal_vector; }

protected:
    /** a vector of pointers to points */
    const std::vector<Point*>& _pnts;
    /** position of pointers to the geometric points */
    size_t _pnt_ids[3];
    bool _initialized;
    double _longest_edge;
    double _normal_vector[3];
};

void getPlaneCoefficients(Triangle const& tri, double c[3]);

std::ostream& operator<<(std::ostream& out, GEOLIB::Triangle const& tri);

}  // end namespace GEOLIB

#endif /* TRIANGLE_H_ */
