#pragma once

#include <iomanip>

#include "field.h"

struct Geometry {
    // mesh geometry
    const int ni;  // number of nodes
    const int nj;
    const int n_pml_xn;
    const int n_pml_xp;
    const int n_pml_yn;
    const int n_pml_yp;
    const Vec2i nn;  // another way to access node counts

    Vec2d dh;  // cell spacing
    Vec2d x0;  // mesh origin
    Vec2d xm;  // origin-diagonally opposite corner (max bound)
    Vec2d xc;  // dm centroid

    Field<double> node_area;  // node volumes

    Geometry();
    Geometry(int ni, int nj, int xn, int xp, int yn, int yp);

    // int GetNi() const;

    // int GetNj() const;

    // Vec2i GetNodeNumber() const;

    // Vec2d GetX0() const;

    // Vec2d GetXm() const;

    // Vec2d GetXc() const;

    // Vec2d GetDh() const;

    // Field<double> GetNodeAreas() const;

    /*functions to set mesh origin and spacing*/
    void SetExtents(const Vec2d x0, const Vec2d xm);

    bool InBounds(Vec2d pos);

    void ComputeNodeAreas();
};
