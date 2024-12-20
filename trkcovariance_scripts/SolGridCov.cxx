﻿#include <TMath.h>
#include <TString.h>
#include <TFile.h>
#include <iostream>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TMatrixDSymEigen.h>
#include <TTree.h>
#include "SolGridCov.h"
#include "../geometry_scripts/SolGeom.h"
#include "SolTrack.h"
//
//
//
SolGridCov::SolGridCov()
{
	//
	// Define pt-polar angle grid
	//
	fNpt = 22;
	fPta.ResizeTo(fNpt);
	Double_t p[] = { 0.1, 0.2, 0.5, 0.7, 1., 2., 3., 4., 6., 8., 10., 15., 
		             20., 25., 30., 40., 50., 60., 80., 100., 150., 200. };
	for (Int_t ip = 0; ip < fNpt; ip++)fPta(ip) = p[ip];
	//
	fNang = 13;
	fAnga.ResizeTo(fNang);
	Double_t a[] = { 10., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90. };
	for (Int_t ia = 0; ia < fNang; ia++)fAnga(ia) = a[ia];
	//
	fCov = new TMatrixDSym **[fNpt];
	//
	for (Int_t ip = 0; ip < fNpt; ip++)
	{
		fCov[ip] = new TMatrixDSym *[fNang];
		for (Int_t ia = 0; ia < fNang; ia++ ) fCov[ip][ia] = 0;
	}
//
}
// Destructor
SolGridCov::~SolGridCov()
{
	for (Int_t ip = 0; ip < fNpt; ip++)
	{
		for (Int_t ia = 0; ia < fNang; ia ++ )fCov[ip][ia]->Clear();
	}
}

void SolGridCov::Write(TString fname, SolGeom *G)
{
	TVectorD pta = fPta;
	//cout << "Test fPta" << endl; pta.Print();
	TVectorD anga = fAnga;
	//cout << "Test fAnga" << endl; anga.Print();
	//
	// Write covariance matrices to file//
	//Open output file and initialize TTree
	//
	TFile *fout = new TFile(fname, "RECREATE");
	TTree *tree = new TTree("treeID", "Covariance matrix tree");
	const Int_t nPar = 5;
	TMatrixDSym Cov(nPar);
	TMatrixDSym *pointer = &Cov;
	//
	// Define and fill branches
	//
	// Store matrices
	Bool_t Res = kTRUE; Bool_t MS = kTRUE;	// Resolution and multiple scattering flags
	for (Int_t ip = 0; ip < fNpt; ip++)		// Loop on pt grid
	{
		Int_t ipt = TMath::Nint(10*pta(ip));
		for (Int_t ia = 0; ia < fNang; ia++)	// Loop on angle grid
		{
			Double_t th = TMath::Pi()*(anga(ia)) / 180.;
			Int_t iang = TMath::Nint(anga(ia));
			TString Bname;
			Bname.Form("C_%dG_%dd", ipt, iang);
			tree->Branch(Bname, "TMatrixDSym", &pointer, 64000, 0);
			//cout << "Branch name: " << Bname << endl;
			//cout << "ipt, iang: " << ipt << ", " << iang << endl;
			Double_t x[3]; Double_t p[3];
			x[0] = 0; x[1] = 0; x[2] = 0;			// Set origin
			p[0] = pta(ip); p[1] = 0; p[2] = pta(ip) / TMath::Tan(th);
			//
			SolTrack *tr = new 	SolTrack(x, p, G);	// Initialize track
			tr->CovCalc(Res, MS);				// Calculate covariance
			Cov = tr->Cov();						// Get covariance
			tree->Fill();							// Store it
		}
	}
	//
	fout->Write();
	fout->Close();
	delete fout;
}

void SolGridCov::Read(TString fname)
{
	//
	// Read covariance matrix grid from file
	//
	TFile *f = new TFile(fname, "READ");
	TTree *T = (TTree*)f->Get("treeID");
	//
	// Read branches
	//
	TVectorD pta = fPta;
	TVectorD anga = fAnga;
	for (Int_t ip = 0; ip < fNpt; ip++)		// Loop on pt grid
	{
		Int_t ipt = TMath::Nint(10*pta(ip));
		for (Int_t ia = 0; ia < fNang; ia++)	// Loop on angle grid
		{
			fCov[ip][ia] = 0;
			Int_t iang = TMath::Nint(anga(ia));
			TString Bname;
			Bname.Form("C_%dG_%dd", ipt, iang);
			T->SetBranchAddress(Bname, &fCov[ip][ia]);
			T->GetEntry(0);
			TMatrixDSym C = *fCov[ip][ia];
			//cout << "Branch name (Read): " << Bname << ", Sig(D) = "
			//	<< 1.e6*sqrt(C(0, 0)) << endl;
		}
	}
	//
	f->Close();
	delete f;
}

//
// Find bin in grid
Int_t SolGridCov::GetMinIndex(Double_t xval, Int_t N, TVectorD x)
{
	Int_t min = -1;	// default for xval below the lower limit
	if (xval < x(0))return min;
	if (xval>x(N - 1)){ min = N; return min; }
	for (Int_t i = 0; i < N; i++) if (xval>x(i))min = i;
	return min;
}
// Force positive definitness in normalized matrix
TMatrixDSym SolGridCov::MakePosDef(TMatrixDSym NormMat)
{
	//
	// Input: symmetric matrix with 1's on diagonal
	// Output: positive definite matrix with 1's on diagonal
	//
	// Default return value
	TMatrixDSym rMatN = NormMat;
	// Check the diagonal
	Bool_t Check = kFALSE;
	Int_t Size = NormMat.GetNcols();
	for (Int_t i = 0; i < Size; i++)if (TMath::Abs(NormMat(i, i) - 1.0)>1.0E-15)Check = kTRUE;
	if (Check)
	{
		cout << "SolGridCov::MakePosDef: input matrix doesn ot have 1 on diagonal. Abort." << endl;
		return rMatN;
	}
	//
	// Diagonalize matrix
	TMatrixDSymEigen Eign(NormMat);
	TMatrixD U = Eign.GetEigenVectors();
	TVectorD lambda = Eign.GetEigenValues();
	// Reset negative eigenvalues to small positive value
	TMatrixDSym D(Size); D.Zero(); Double_t eps = 1.0e-13;
	for (Int_t i = 0; i < Size; i++)
	{
		D(i, i) = lambda(i);
		if (lambda(i) <= 0) D(i, i) = eps;
	}
	//Rebuild matrix
	TMatrixD Ut(TMatrixD::kTransposed, U);
	TMatrixD rMat = (U*D)*Ut;				// Now it is positive defite
	// Restore all ones on diagonal
	for (Int_t i1 = 0; i1 < Size; i1++)
	{
		Double_t rn1 = TMath::Sqrt(rMat(i1, i1));
		for (Int_t i2 = 0; i2 <= i1; i2++)
		{
			Double_t rn2 = TMath::Sqrt(rMat(i2, i2));
			rMatN(i1, i2) = 0.5*(rMat(i1, i2) + rMat(i2, i1)) / (rn1*rn2);
			rMatN(i2, i1) = rMatN(i1, i2);
		}
	}
	return rMatN;
}
//
// Interpolate covariance matrix: Bi-linear interpolation
//
TMatrixDSym SolGridCov::GetCov(Double_t pt, Double_t ang)
{
	//
	// pt in GeV and ang in degrees
	Int_t minPt = GetMinIndex(pt, fNpt, fPta);
	if (minPt == -1)minPt = 0;
	if (minPt == fNpt - 1)minPt = fNpt - 2;
	Double_t dpt = fPta(minPt + 1) - fPta(minPt);
	//
	// Put ang in 0-90 range
	ang = TMath::Abs(ang);
	while (ang > 90.)ang -= 90.;	// Needs to be fixed
	Int_t minAng = GetMinIndex(ang, fNang, fAnga);
	if (minAng == -1)minAng = 0;
	if (minAng == fNang - 1)minAng = fNang - 2;
	Double_t dang = fAnga(minAng + 1) - fAnga(minAng);
	//
	Double_t tpt = (pt - fPta(minPt)) / dpt;
	Double_t tang = (ang - fAnga(minAng)) / dang;
	//
	TMatrixDSym C11 = *fCov[minPt][minAng];
	TMatrixDSym C12 = *fCov[minPt][minAng + 1];
	TMatrixDSym C21 = *fCov[minPt + 1][minAng];
	TMatrixDSym C22 = *fCov[minPt + 1][minAng + 1];
	TMatrixDSym Cv = ((1-tpt) * (1-tang)) *C11 +
		             ((1-tpt) *    tang ) *C12 +
		             (tpt     * (1-tang)) *C21 +
		             (tpt    *     tang ) *C22;
	//
	// Check for positive definiteness
	//
	TMatrixDSym CvN = Cv;
	TMatrixDSym DCvInv(5); DCvInv.Zero();
	for (Int_t id = 0; id < 5; id++) DCvInv(id, id) = 1.0 / TMath::Sqrt(Cv(id, id));
	CvN.Similarity(DCvInv);	// Normalize diagonal to 1
	TDecompChol Chl(CvN);
	if (!Chl.Decompose())
	{
		cout << "SolGridCov::GetCov: Interpolated matrix not positive definite. Recovering ...." << endl;
		TMatrixDSym rCv = MakePosDef(CvN); CvN = rCv;
		TMatrixDSym DCv(5); DCv.Zero();
		for (Int_t id = 0; id < 5; id++) DCv(id, id) = TMath::Sqrt(Cv(id, id));
		Cv = CvN.Similarity(DCv);	// Normalize diagonal to 1
	}
	//
	return Cv;
}
//
