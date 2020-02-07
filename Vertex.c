#include <TMath.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TRandom.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include "SolGridCov.h"
#include "ObsTrk.h"
//
TVectorD Fillf(TVectorD par, TVectorD xin)
{
	//
	// Decode input arrays
	//
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	Double_t x = xin(0);
	Double_t y = xin(1);
	Double_t z = xin(2);
	Double_t R2 = x * x + y * y;
	//
	// Calculate constraints
	//
	TVectorD f(2);
	f(0) = R2 * C + D - y * TMath::Cos(p0) + x * TMath::Sin(p0);
	f(1) = z0 - z + ct * TMath::Sqrt(R2 - D * D);
	//
	return f;
}
//

TMatrixD FillD(TVectorD par, TVectorD xin)
{
	//
	// Decode input arrays
	//
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	Double_t x = xin(0);
	Double_t y = xin(1);
	Double_t z = xin(2);
	Double_t R2 = x * x + y * y;
	//
	// Calculate matrix elements
	//
	TMatrixD Do(2, 5);
	Do(0, 0) = 1.0;
	Do(0, 1) = x * TMath::Cos(p0) + y * TMath::Sin(p0);
	Do(0, 2) = R2;
	Do(0, 3) = 0.0;
	Do(0, 4) = 0.0;
	Do(1, 0) = -ct * D / TMath::Sqrt(R2 - D * D);
	Do(1, 1) = 0.0;
	Do(1, 2) = 0.0;
	Do(1, 3) = 1.0;
	Do(1, 4) = TMath::Sqrt(R2 - D * D);
	//
	return Do;
}
//
TMatrixD FillB(TVectorD par, TVectorD xin)
{
	//
	// Decode input arrays
	//
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	Double_t x = xin(0);
	Double_t y = xin(1);
	Double_t z = xin(2);
	Double_t R2 = x * x + y * y;
	//
	// Calculate constraints
	//
	TMatrixD B(2, 3);
	B(0, 0) = 2 * C*x + TMath::Sin(p0);
	B(0, 1) = 2 * C*y - TMath::Cos(p0);
	B(0, 2) = 0.0;
	B(1, 0) = ct * x / TMath::Sqrt(R2 - D * D);
	B(1, 1) = ct * y / TMath::Sqrt(R2 - D * D);
	B(1, 2) = -1;
	//
	return B;
}
//
TVectorD Vertex0(Int_t Ntr, ObsTrk **tracks)
{
	//
	// Preliminary estimate of the vertex position
	// based on transformation of track into points
	// and vertices into lines
	// No steering of track parameters
	// No error calculation
	//
	TVectorD xv(3);		// returned vertex position
	//
	TMatrixDSym H(2); 
	TVectorD xvt(2);
	TVectorD cxy(2); 
	//
	// Loop on tracks for transverse fit
	//
	TVectorD x0(2); x0.Zero();
	Double_t Rv = 0.0;	    // Radius of first iteration
	Int_t Ntry = 0;
	Int_t TryMax = 10;
	Double_t epsi = 1000.;	// Starting stability
	Double_t eps = 0.01;		// vertex stability required
	while (epsi > eps && Ntry < TryMax)
	{
		H.Zero(); cxy.Zero();
		for (Int_t i = 0; i < Ntr; i++)
		{
			// Get track helix parameters and their covariance matrix 
			ObsTrk *t = tracks[i];
			TVectorD par = t->GetObsPar();
			TMatrixDSym C = t->GetCov();
			//
			// Transverse fit
			Double_t D0i = par(0);
			Double_t phi = par(1);
			Double_t Ci = par(2);
			Double_t Di = (D0i*(1. + Ci*D0i) + Rv*Rv*Ci) / (1. + 2 * Ci*D0i);
			Double_t sDi2 = C(0, 0);
			Double_t sn = TMath::Sin(phi);
			Double_t cs = TMath::Cos(phi);
			H(0, 0) += sn*sn / sDi2;
			H(0, 1) += -cs*sn / sDi2;
			H(1, 0) += -cs*sn / sDi2;
			H(1, 1) += cs*cs / sDi2;
			cxy(0) += -Di*sn / sDi2;
			cxy(1) += Di*cs / sDi2;
		}
		//
		TMatrixDSym H0 = H;
		H.Invert();
		xvt = H*cxy;
		xv(0) = xvt(0); xv(1) = xvt(1);
		Rv = TMath::Sqrt(xv(0)*xv(0) + xv(1)*xv(1));
		TVectorD dx = xvt - x0;
		epsi = H0.Similarity(dx);
		x0 = xvt;
		Ntry++;
		//cout << "Vtx0: Iteration #" << Ntry << ", eps = " << epsi << ", x = " << xv(0) << ", y = " << xv(1) << endl;
	}
	//
	// Longitudinal fit
	Double_t Rv2 = Rv*Rv;
	//
	// Loop on tracks for longitudinal fit
	Double_t hz = 0.0;
	Double_t cz = 0.0;
	for (Int_t i = 0; i < Ntr; i++)
	{
		// Get track helix parameters and their covariance matrix 
		ObsTrk *t = tracks[i];
		TVectorD par = t->GetObsPar();
		TMatrixDSym C = t->GetCov();
		//
		// Longitudinal fit
		Double_t zi = par(3);
		Double_t cti = par(4);
		Double_t Di = par(0);
		Double_t sZi2 = C(3, 3);
		//
		hz += 1 / sZi2;
		cz += (cti*TMath::Sqrt(Rv2 - Di*Di) + zi)/sZi2;
	}
	xv(2) = cz / hz;
	//
	return xv;
}
//
Double_t Vertex(Int_t Ntr, ObsTrk **tracks, TVectorD &x, TMatrixDSym &covX)
{
	//
	// Get approximate vertex evaluation
	//
	TVectorD x0 = Vertex0(Ntr, tracks);
	//cout << "Preliminary vertex" << endl; x0.Print();
	TVectorD dx(3);		// Solution x variation
	//
	TVectorD f(2);		// Constraints
	TMatrixD D(2, 5);	// df/d alf (constraint over parameters)
	TMatrixD B(2, 3);	// df/dx    (constraint over x) 
	// Stored quantities
	TVectorD **fi = new TVectorD*[Ntr];
	TMatrixD **Bi = new TMatrixD*[Ntr];
	TMatrixDSym **Wi = new TMatrixDSym*[Ntr];
	//
	// Loop on tracks to calculate everything
	//
	Int_t Ntry = 0;
	Int_t TryMax = 10;
	Double_t eps = 0.01; // vertex stability
	Double_t epsi = 1000.;
	x = x0;
	while (epsi > eps && Ntry < TryMax)		// Iterate until found vertex is stable
	{
		TVectorD BtWf(3); BtWf.Zero();
		covX.Zero();		// Reset vertex covariance
		// 
		for (Int_t i = 0; i < Ntr; i++)
		{
			// Get track helix parameters and their covariance matrix 
			ObsTrk *t = tracks[i];
			TVectorD par = t->GetObsPar();
			TMatrixDSym C = t->GetCov();
			// Fill D
			D = FillD(par, x);
			// Fill B
			B = FillB(par, x);
			Bi[i] = new TMatrixD(B);
			//cout << "Bi" << endl; Bi[i]->Print();
			// Fill constraints
			f = Fillf(par, x);
			fi[i] = new TVectorD(f);
			//cout << "fi" << endl; fi[i]->Print();
			//
			TMatrixDSym W = C.Similarity(D);
			W.Invert();
			Wi[i] = new TMatrixDSym(W);
			//cout << "Wi" << endl; Wi[i]->Print();
			TMatrixD Bt(TMatrixD::kTransposed, B);
			TMatrixDSym W1(W);
			TMatrixDSym BtWB = W1.Similarity(Bt);
			covX += BtWB;
			BtWf += Bt * (W*f);
		}
		TMatrixDSym Hess = covX;
		covX.Invert();
		dx = (-1.0*covX) * BtWf;
		x += dx;
		epsi = Hess.Similarity(dx);
		Ntry++;
		//cout << "Vtx: Iteration #"<<Ntry<<", eps = "<<epsi<<", x = " << x(0) << ", " << x(1) << ", " << x(2) << endl;
	}
	//
	// Calculate Chi2
	//
	Double_t Chi2 = 0.0;
	for (Int_t i = 0; i < Ntr; i++)
	{
		TVectorD lambda = *fi[i] + (*Bi[i]) * dx;
		TMatrixDSym Wp = *Wi[i];
		Chi2 += Wp.Similarity(lambda);
	}
	//
	return Chi2;
}
Double_t fchi2(Double_t *x, Double_t *p)
{
	Double_t Ndof = p[0];
	Double_t Norm = p[1];
	Double_t z = x[0];
	//
	Double_t ex = Ndof / 2.0 - 1.0;
	Double_t chi2 = pow(z, ex)*TMath::Exp(-z / 2.0);
	Double_t bin = 5.*Ndof / 100.;
	chi2 = Norm * chi2 * bin / (pow(2, Ndof / 2.0)*TMath::Gamma(Ndof / 2.0));
	return chi2;
}
//
// Testing program
//
void TestVtx(Int_t Nvtx = 100, Int_t Ntr = 2)
{
	//
	// Init geometry
	Double_t Bfield = 2.0;	// Magnetic field
	SolGridCov *GC = new SolGridCov();
	GC->Read("CovIDEA-BASE.root");			// Read in covariance array
	//
	// Ranges
	//
	Double_t ThDegMin = 20.0;
	Double_t ThDegMax = 160.0;
	Double_t Lmin = 0.010;
	Double_t Lmax = 0.014;
	Double_t dTheta = 0.010;
	Double_t dPhi = 0.010;
	Double_t Pmin = 1.0;
	Double_t Pmax = 2.0;
	//
	// Histograms
	TH1D *hXpull = new TH1D("hXpull", "Pull X vertex component", 100, -5., 5.);
	TH1D *hYpull = new TH1D("hYpull", "Pull Y vertex component", 100, -5., 5.);
	TH1D *hZpull = new TH1D("hZpull", "Pull Z vertex component", 100, -5., 5.);
	Double_t Ndof = 2.0 * Ntr - 3.0;
	TH1D *hChi2 = new TH1D("hChi2", "Vertex #chi^{2}", 100, 0., 5 * Ndof);
	TF1 *fch = new TF1("fch", fchi2, 0., 5 * Ndof, 2);
	fch->SetParameter(0, Ndof);
	fch->SetParameter(1, Nvtx);
	//
	// Loop on # vertices
	//
	for (Int_t n = 0; n < Nvtx; n++)
	{
		//
		// Extract a direction
		Double_t ThMin = ThDegMin * TMath::Pi() / 180.;
		Double_t ThMax = ThDegMax * TMath::Pi() / 180.;
		Double_t rnTh = gRandom->Rndm();
		Double_t Th = ThMin + rnTh * (ThMax - ThMin);
		Double_t rnPh = gRandom->Rndm();
		Double_t Ph = 2 * TMath::Pi()*rnPh;
		//
		// Extract a flight distance (m)
		Double_t rnL = gRandom->Rndm();
		Double_t Lvtx = Lmin + rnL * (Lmax - Lmin);
		//
		TVector3 x;
		x(0) = Lvtx * TMath::Sin(Th)*TMath::Cos(Ph);
		x(1) = Lvtx * TMath::Sin(Th)*TMath::Sin(Ph);
		x(2) = Lvtx * TMath::Cos(Th);
		//
		//cout << "True vertex: x = " << x(0) << ", y = " << x(1) << ", z = " << x(2) << endl;
		//
		// Loop on tracks
		ObsTrk **tracks = new ObsTrk*[Ntr];
		for (Int_t i = 0; i < Ntr; i++)
		{
			Double_t rnP = gRandom->Rndm();
			Double_t Ptot = Pmin + rnP * (Pmax - Pmin);
			Double_t ThP = gRandom->Gaus(Th, dTheta);
			Double_t PhP = gRandom->Gaus(Ph, dPhi);
			//
			TVector3 P;
			P(0) = Ptot * TMath::Sin(ThP)*TMath::Cos(PhP);
			P(1) = Ptot * TMath::Sin(ThP)*TMath::Sin(PhP);
			P(2) = Ptot * TMath::Cos(ThP);
			Double_t Q = 1.0;
			if (gRandom->Rndm() > 0.5)Q = -1.0;
			tracks[i] = new ObsTrk(x, P, Q, Bfield, GC);
			//TVectorD par = tracks[i]->GetObsPar();
			//cout << "i = " << i << ", loading par = " << endl; par.Print();
		}
		//Double_t xa[3];
		//x.GetXYZ(xa);
		TVectorD xvtx(3); xvtx.Zero();
		//cout << "xvtx = " << endl; xvtx.Print();
		TMatrixDSym covX(3);
		Double_t Chi2 = Vertex(Ntr, tracks, xvtx, covX);
		//
		//cout << "Fit vertex: x = " << xvtx(0) << ", y = " << xvtx(1) << ", z = " << xvtx(2) << endl;
		Double_t PullX = (xvtx(0) - x(0)) / TMath::Sqrt(covX(0, 0));
		hXpull->Fill(PullX);
		Double_t PullY = (xvtx(1) - x(1)) / TMath::Sqrt(covX(1, 1));
		hYpull->Fill(PullY);
		Double_t PullZ = (xvtx(2) - x(2)) / TMath::Sqrt(covX(2, 2));
		hZpull->Fill(PullZ);
		hChi2->Fill(Chi2);
	}
	//
	// Plots
	TCanvas *cnv = new TCanvas("cnv", "Normalized residuals", 20, 20, 500, 500);
	cnv->Divide(2, 2);
	cnv->cd(1);
	hXpull->Fit("gaus");
	hXpull->Draw();
	cnv->cd(2);
	hYpull->Fit("gaus");
	hYpull->Draw();
	cnv->cd(3);
	hZpull->Fit("gaus");
	hZpull->Draw();
	cnv->cd(4);
	hChi2->Fit("fch");
	hChi2->Draw();
}

