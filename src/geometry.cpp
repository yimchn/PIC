#include "geometry.h"

Geometry::Geometry()
    : ni{26},
      nj{26},
      nn{26, 26},
      n_pml_xn(5),
      n_pml_xp(5),
      n_pml_yn(5),
      n_pml_yp(5),
      node_area{26, 26} {}

Geometry::Geometry(int ni, int nj, int xn, int xp, int yn, int yp)
    : ni(ni),
      nj(nj),
      n_pml_xn(xn),
      n_pml_xp(xp),
      n_pml_yn(yn),
      n_pml_yp(yp),
      nn{ni, nj},
      node_area{ni, nj} {}

// int Geometry::GetNi() const { return ni; }

// int Geometry::GetNj() const { return nj; }

// Vec2i Geometry::GetNodeNumber() const { return nn; }

// Vec2d Geometry::GetX0() const { return Vec2d(x0); }

// Vec2d Geometry::GetXm() const { return Vec2d(xm); }

// Vec2d Geometry::GetXc() const { return Vec2d(xc); }

// Vec2d Geometry::GetDh() const { return Vec2d(dh); }

// Field<double> Geometry::GetNodeAreas() const { return node_area; }

bool Geometry::InBounds(Vec2d pos) {
    for (int i = 0; i < 2; i++)
        if (pos[i] < x0[i] || pos[i] >= xm[i]) return false;
    return true;
}

/*computes node volumes, dx*dy*dz on internal nodes and fractional
 * values on dm boundary faces*/
void Geometry::ComputeNodeAreas() {
    for (int i = 0; i < ni; i++)
        for (int j = 0; j < nj; j++) {
            double area = dh[0] * dh[1];  // default volume
            if (i == 0 || i == ni - 1)
                area *= 0.5;  // reduce by two for each boundary index
            if (j == 0 || j == nj - 1) area *= 0.5;
            node_area[i][j] = area;
        }
}

/*sets dm bounding box and computes mesh spacing*/
void Geometry::SetExtents(Vec2d _x0, Vec2d _xm) {
    /*set origin and the opposite corner*/
    x0 = _x0;
    xm = _xm;

    /*compute spacing by dividing length by the number of cells*/
    for (int i = 0; i < 2; i++) {
        dh[i] = (xm(i) - x0(i)) / (nn(i) - 1);
    }

    // compute centroid
    xc = 0.5 * (x0 + xm);

    /*recompute node volumes*/
    ComputeNodeAreas();
}
