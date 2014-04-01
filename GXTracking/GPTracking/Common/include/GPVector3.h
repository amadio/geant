#ifndef GPVector3_H
#define GPVector3_H

#include "GPConstants.h"

class GPVector3 {

public:

  FQUALIFIER GPVector3() : px(0.0), py(0.0), pz(0.0) {}
  FQUALIFIER GPVector3(G4double x, double y, double z) : px(x), py(y), pz(z) {}
  GPVector3(const GPVector3& r) : px(r.Px()), py(r.Py()), pz(r.Pz()) {}
  ~GPVector3() {}

  FQUALIFIER void Set(G4double x, double y, double z) { px = x; py = y; pz = z; }
  FQUALIFIER void SetX(G4double x) { px = x; }
  FQUALIFIER void SetY(G4double y) { py = y; }
  FQUALIFIER void SetZ(G4double z) { pz = z; }
  FQUALIFIER void SetPx(G4double x) { px = x; }
  FQUALIFIER void SetPy(G4double y) { py = y; }
  FQUALIFIER void SetPz(G4double z) { pz = z; }
  FQUALIFIER G4double X() const { return px; }
  FQUALIFIER G4double Y() const { return py; }
  FQUALIFIER G4double Z() const { return pz; }
  FQUALIFIER G4double Px() const { return px; }
  FQUALIFIER G4double Py() const { return py; }
  FQUALIFIER G4double Pz() const { return pz; }
  FQUALIFIER G4double Dot(const GPVector3& p) { return px*p.Px() + py*p.Py() + pz*p.Pz(); }
  FQUALIFIER GPVector3& operator=(const GPVector3& right) { px = right.Px(); py = right.Py(); pz = right.Pz(); return *this; }
  FQUALIFIER GPVector3 operator+(const GPVector3& right) { return GPVector3(px+right.Px(), py+right.Py(), pz+right.Pz()); }
  FQUALIFIER GPVector3 operator-(const GPVector3& right) { return GPVector3(px-right.Px(), py-right.Py(), pz-right.Pz()); }
  FQUALIFIER GPVector3 operator*(G4double a) { return GPVector3(a*Px(),a*Py(),a*Pz()); }
  //  FQUALIFIER GPVector3 operator*(G4double a, GPVector3& right) { return GPVector3(a*right.Px(),a*right.Py(),a*right.Pz()); }
  //  FQUALIFIER G4double operator*(GPVector3& a, GPVector3& b) { return a.Dot(b); }
  FQUALIFIER GPVector3& operator *= (G4double a) { px *= a; py *= a; pz *= a; return *this; }
  FQUALIFIER GPVector3& operator += (const GPVector3& p) { px += p.Px(); py += p.Py(); pz += p.Pz(); return *this; }
  FQUALIFIER GPVector3& operator -= (const GPVector3& p) { px -= p.Px(); py -= p.Py(); pz -= p.Pz(); return *this; }
  FQUALIFIER bool operator == (const GPVector3& p) { return (px == p.Px() && py == p.Py() && pz == p.Pz()); }
  FQUALIFIER bool operator != (const GPVector3& p) { return (px != p.Px() || py != p.Py() || pz != p.Pz()); }

  FQUALIFIER GPVector3 Unit() const {
    G4double tot = Mag2();
    GPVector3 p(px,py,pz);
    return tot > 0.0 ? p *= (1.0/sqrt(tot)) : p;
  }
  FQUALIFIER G4double Phi() const { return px==0.0 && py==0.0 ? 0.0 : atan2(px,py); }
  FQUALIFIER G4double Theta() const { return px==0.0 && py==0.0 && pz==0.0 ? 0.0 : atan2(Pt(),pz); }
  FQUALIFIER G4double Eta() const {
    G4double cosTheta = cos(Theta());
   if (cosTheta*cosTheta < 1.0) return -0.5* log( (1.0-cosTheta)/(1.0+cosTheta) );
   if (pz > 0.0) return 1.0e10;
   else        return -1.0e10;
  }

  FQUALIFIER G4double Mag() const { return sqrt(Mag2()); }
  FQUALIFIER G4double Mag2() const { return px*px + py*py + pz*pz; }
  FQUALIFIER G4double Pt() const { return sqrt(px*px + py*py); }

  FQUALIFIER void SetMag(G4double m) {
    G4double factor = Mag();
    if(factor > 0.0) {
      factor = m/factor;
      Set(px*factor,py*factor,pz*factor);
    }
  }

  FQUALIFIER void RotateUz(const GPVector3& UzVector) {
   // NewUzVector must be normalized !
    GPVector3 NewUzVector = UzVector.Unit();
    G4double u1 = NewUzVector.Px();
    G4double u2 = NewUzVector.Py();
    G4double u3 = NewUzVector.Pz();
    G4double up = u1*u1 + u2*u2;

    if (up) {
      up = sqrt(up);
      G4double x = px,  y = py,  z = pz;
      px = (u1*u3*x - u2*y + u1*up*z)/up;
      py = (u2*u3*x + u1*y + u2*up*z)/up;
      pz = (u3*u3*x -    x + u3*up*z)/up;
    } else if (u3 < 0.) { px = -px; pz = -pz; }      // phi=0  teta=pi
    else {};
  }


private:
  G4double px;
  G4double py;
  G4double pz;

};

#endif
