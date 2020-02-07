#include <iostream>
#include "SolGeom.h"
#include "SolGridCov.h"
#include "ObsTrk.h"

#include "TVector3.h"
#include "TLorentzVector.h"



using namespace std;


// needed input:


// - input Covariance Matrix in root file form
// 
// - provide vector position at origin, and vector momentum
//



int main() {

  // initialize Geometry
  SolGeom *G = new SolGeom();
  Double_t Bfield = G->B();

  cout<<Bfield<<endl;

  // initialize tracking resolution
  SolGridCov *GC = new SolGridCov();

  // reads covariance array
  GC->Read("CovIDEA-BASE.root");

  // apply track Resolution

  TVector3 tX(0., 0., 0.);
  TLorentzVector p4;
  Double_t pt = 10.;
  Double_t m = 0.1;
  Double_t e = TMath::Sqrt(pt*pt + m*m);



  p4.SetPtEtaPhiM(pt,0.,0.,e);

  TVector3 tP1 = p4.Vect();
  Double_t Q1 = 1; // charge

  ObsTrk *Tr = new ObsTrk(tX, tP1, Q1, Bfield, GC);

  TVector3 obsP1 = Tr->GetObsP();
  Double_t Eobs = TMath::Sqrt(m*m + obsP1.Mag2());

  TLorentzVector p4obs;
  p4obs.SetPxPyPzE(obsP1.Px(), obsP1.Py(), obsP1.Pz(), e);

  cout<<p4.Pt()<<","<<p4.Eta() <<","<<p4.Phi() <<","<<p4.E() <<endl;
  cout<<p4obs.Pt()<<","<<p4obs.Eta() <<","<<p4obs.Phi() <<","<<p4obs.E() <<endl;

}
