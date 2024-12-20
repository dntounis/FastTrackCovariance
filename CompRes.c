﻿#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TString.h>
#include <THStack.h>
#include <TFile.h>
#include <iostream>
#include "geometry_scripts/SolGeom.h"
#include "trkcovariance_scripts/SolTrack.h"
#include "trkcovariance_scripts/SolGridCov.h"
//
void CompRes(Double_t Ang)
{
	//**********************
	// Manage input files  *
	//**********************
	//
	// Available angles:
	const Int_t nAng = 8;
	Double_t theta[nAng] = { 10., 20., 30., 40., 45., 60., 75., 90. };	// angles
	Bool_t found = kFALSE; Int_t nA = -1;
	for (Int_t i = 0; i < nAng; i++)	if (Ang == theta[i])
	{
		found = kTRUE; 
		nA = i;
	}
	if (!found)
	{
		cout << "Angle not available. Choose one among the following:" << endl;
		cout << "10, 20, 30, 40, 45, 60, 75, 90" << endl;
		return;
	}
	Double_t th = TMath::Pi()*Ang / 180.;
	const Int_t nPar = 5;
	TString sPar[nPar] = { "d0", "Phi", "Pt", "Theta", "Z" };
	TString sAng[nAng] = { "10deg", "20deg", "30deg", "40deg",
		"45deg", "60deg", "75deg", "90deg" };
	TString npt = "g" + sPar[2] + "-" + sAng[nA];
	TString nd0 = "g" + sPar[0] + "-" + sAng[nA];
	TString nz0 = "g" + sPar[4] + "-" + sAng[nA];
	TString nth = "g" + sPar[3] + "-" + sAng[nA];
	//
	// Input  file
	TFile *InFile = new TFile("./dump/graphs.root", "READ");
	TGraphErrors *grept = (TGraphErrors *)InFile->Get(npt);
	TGraphErrors *gred0 = (TGraphErrors *)InFile->Get(nd0);
	TGraphErrors *grez0 = (TGraphErrors *)InFile->Get(nz0);
	TGraphErrors *greth = (TGraphErrors *)InFile->Get(nth);
	//
	// Input graphs loaded
	//*****************************************************
	//
	//**************************
	//	Initialize geometry    *
	//**************************
	//
	Bool_t Res = kTRUE;		// Enable detector resolution effects
	Bool_t MS  = kTRUE;		// Enable multiple scattering
	SolGeom *G;				// Initialize geometry
	const Int_t nDet = 9;
	Bool_t OK[nDet] = {		// Enable selected parts of the detector 
		1,					// Beam pipe
		1,					// Inner VTX pixel layers
		1,					// Outer VTX layers
		1,					// Drift chamber
		1,					// Barrel Si wrapper
		1,					// Barrel pre-shower
		1,					// Forw. VTX pixel layers
		1,					// Forw. Si wrapper
		1};					// Forw. pre-shower
	//G = new SolGeom("GeoIDEA_GT.txt",2.0);	// Geometry with selected detectors
	G = new SolGeom("geometry_files/GeoIDEA_GT.txt",2.0);	// Geometry with selected detectors
	
	G->Draw();				// Draw R-z geometry
	//char* fname = "GeoIDEA_BASE.txt";
	//G->GeoPrint(fname);
	//
	// Geometry initialized and drawn
	//********************************************************
	//
	//*****************************************
	// Draw track on top of geometry display  *
	//*****************************************
	//	
	TCanvas *cc = G->cnv();					// Get canvas with geo display
	Double_t x[3] = { 0.0, 0.0, 0.0 };		// Track starting point
	Double_t ppt = 0.7;						// Track pt
	Double_t ppz = ppt / TMath::Tan(th);	// Track pz
	Double_t p[3] = { ppt, 0.0, ppz };		// Track momentum
	SolTrack *trk = new SolTrack(x, p, G);	// Initialize track
	TGraph *gr = trk->TrkPlot();			// graph intersection with layers
	gr->Draw("PLSAME");						// plot track
	cc -> SaveAs("IDEA_geometry_track.pdf");
	//
	// End track plot
	//*********************************************************
	//
	// *****************************
	// Make plot of material       *
	//******************************
	// 
	THStack *hMat = new THStack("hMat", "IDEA: Material vs. cos(#theta)");
	TH1D *hPipe = new TH1D("hPisa", "IDEA: Pipe material vs. Theta", 100, 0., 1.);
	hPipe->SetFillColor(kRed);
	hPipe->SetLineColor(kBlack);
	hPipe->GetYaxis()->SetRange(0.0, 100.);
	TH1D *hSilVtx = new TH1D("hSilVtx", "IDEA: Silicon VTX material vs. Theta", 100, 0., 1.);
	hSilVtx->SetFillColor(kGreen);
	hSilVtx->SetLineColor(kBlack);
	hSilVtx->GetYaxis()->SetRange(0.0, 100.);
	TH1D *hDCH = new TH1D("hDCH", "IDEA: Drift chamber material vs. Theta", 100, 0., 1.);
	hDCH->SetFillColor(kCyan);
	hDCH->SetLineColor(kBlack);
	hDCH->GetYaxis()->SetRange(0.0, 100.);
	TH1D *hSilWrp = new TH1D("hSilWrp", "IDEA: Silicon wrapper material vs. Theta", 100, 0., 1.);
	hSilWrp->SetFillColor(kMagenta);
	hSilWrp->SetLineColor(kBlack);
	hSilWrp->GetYaxis()->SetRange(0.0, 100.);
	hSilWrp->GetXaxis()->SetTitle("#h");
	hSilWrp->GetYaxis()->SetTitle("% radiation length");
	Int_t nStep = hPipe->GetNbinsX();
	for (Int_t i = 1; i <= nStep; i++)
	{
		Double_t CosTh = hPipe->GetBinCenter(i);
		Double_t th = TMath::ACos(CosTh);
		Double_t *mat = new Double_t[14];
		mat = G->FracX0(th);
		Double_t mPipe = 100.*mat[0]; hPipe->SetBinContent(i, mPipe);
		Double_t mSilVtx = 100 * (mat[1] + mat[2] + mat[9]); hSilVtx->SetBinContent(i, mSilVtx);
		Double_t mDCH = 100.*(mat[3] + mat[4] + mat[5] + mat[10]); hDCH->SetBinContent(i, mDCH);
		Double_t mSilWrp = 100.*(mat[6] + mat[11]); hSilWrp->SetBinContent(i, mSilWrp);
	}
	hMat->Add(hPipe);
	hMat->Add(hSilVtx);
	hMat->Add(hDCH);
	hMat->Add(hSilWrp);
	TCanvas *cmat = new TCanvas("cmat", "IDEA material in tracking", 100, 100, 500, 500);
	cmat->cd(1);
	//TLegend *lg = new TLegend(0.5, 0.6, 0.9, 0.9);
	TLegend *lg = new TLegend(0.1, 0.9, 0.5, 0.7);
	lg->AddEntry(hPipe, "Beam pipe", "f");
	lg->AddEntry(hSilVtx, "Vertex silicon", "f");
	lg->AddEntry(hDCH, "Drift chamber", "f");
	lg->AddEntry(hSilWrp, "Silicon wrapper", "f");
	hMat->SetMaximum(30.);
	hMat->Draw(); lg->Draw();
	cmat->SaveAs("IDEA_material_tracking.pdf");
	//
	//******************************************************
	// Compare track parameter resolutions vs pt and theta *
	//******************************************************
	//
	TCanvas *resol = new TCanvas("resol", "Resolutions", 100, 100, 500, 500);
	resol->Divide(2, 2);
	// Define graphs
	TGraph *grpt;				// pt resolution graphs
	TGraph *grd0;				// D resolution graphs
	TGraph *grz0;				// z0 resolution graphs
	TGraph *grth;				// theta resolution
	// Setup graph arrays
	Int_t Npt = 200;			// Nr. of points per graph
	Double_t * pt = new Double_t[Npt];
	Double_t * pp = new Double_t[Npt];
	Double_t *spt = new Double_t[Npt];
	Double_t *sd0 = new Double_t[Npt];
	Double_t *sz0 = new Double_t[Npt];
	Double_t *sth = new Double_t[Npt];
	// Fill graph arrays
	Double_t ptmin = 0.5;
	Double_t ptmax = 100;
	Double_t pts = (ptmax - ptmin) / (Double_t)(Npt-1);
	for (Int_t k = 0; k < Npt; k++)	// Loop on pt
	{
		Double_t x[3]; Double_t p[3];
		x[0] = 0; x[1] = 0; x[2] = 0;			// Set origin
		pt[k] = ptmin+k*pts;					// Set transverse momentum
		p[0] = pt[k]; p[1] = 0;	p[2] = pt[k] / TMath::Tan(th);
		pp[k] = pt[k]/TMath::Sin(th);			// Set momentum
		SolTrack *tr = new 	SolTrack(x, p, G);	// Initialize track
		Int_t nH = tr->nHit();
		//cout << "Pt = " << pt[k] << ", #hits = " << nH << endl;
		tr->CovCalc(Res,MS);					// Calculate covariance
		spt[k] = tr->s_pt();							// Dpt/pt
		sd0[k] = tr->s_D()*1e6;							// D  res. - change to microns
		sz0[k] = tr->s_z0()*1e6;						// z0 res. - change to microns
		sth[k] = tr->s_ct() / (1 + pow(tr->ct(), 2));	// theta resolution
	}
	//
	// Compare pt resolution
	resol->cd(1);
	grpt = new TGraph(Npt, pt, spt);			// Estimated resolution
	grpt->SetLineColor(kRed);
	grpt->SetTitle("#sigma_{pt}/pt");
	grept->SetTitle("#sigma_{pt}/pt");
	grept->SetMinimum(0.0);
	grpt->SetMinimum(0.0);
	grpt->GetXaxis()->SetTitle("pt (GeV)");
	grept->Draw("AP");							// Simulated resolution
	grpt->Draw("SAME");							// Estimated resolution
	// Compare d0 resolution
	resol->cd(2);
	grd0 = new TGraph(Npt, pp, sd0);			// Estimated resolution
	grd0->SetLineColor(kRed);
	grd0->SetTitle("D_{0} (#mum)");
	gred0->SetTitle("D_{0} (#mum)");
	gred0->SetMinimum(0.0);
	grd0->SetMinimum(0.0);
	grd0->GetXaxis()->SetTitle("p (GeV)");
	gred0->Draw("AP");							// Simulated resolution
	grd0->Draw("SAME");							// Estimated resolution
	// Compare z0 resolution
	resol->cd(3);
	grz0 = new TGraph(Npt, pp, sz0);			// Estimated resolution
	grz0->SetLineColor(kRed);
	grz0->SetTitle("Z_{0} (#mum)");
	grez0->SetTitle("Z_{0} (#mum)");
	grez0->SetTitle("Z_{0} (#mum)");
	grez0->SetMinimum(0.0);
	grz0->GetXaxis()->SetTitle("p (GeV)");
	grez0->Draw("AP");							// Simulated resolution
	grz0->SetMarkerColor(kRed);
	grz0->SetMarkerSize(0.5);
	grz0->SetMarkerStyle(kFullCircle);
	grz0->Draw("PLSAME");						// Estimated resolution
	// Compare theta resolution
	resol->cd(4);
	grth = new TGraph(Npt, pp, sth);			// Estimated resolution
	grth->SetLineColor(kRed);
	grth->SetTitle("#theta (rad)");
	greth->SetTitle("#theta (rad)");
	greth->SetMinimum(0.0);
	grth->SetMinimum(0.0);
	grth->GetXaxis()->SetTitle("p (GeV)");
	greth->Draw("AP");							// Simulated resolution
	grth->SetMarkerColor(kRed);
	grth->SetMarkerSize(0.5);
	grth->SetMarkerStyle(kFullCircle);
	grth->Draw("PLSAME");						// Estimated resolution
	resol->SaveAs("IDEA_resolution_tracking.pdf");
	//
	//***************************************************
	// Repeat using interpolation                       *
	//***************************************************
	//
	TCanvas *resol1 = new TCanvas("resol1", "Resolutions from grid", 10, 10, 500, 500);
	resol1->Divide(2, 2);
	// Define graphs
	TGraph *ggrpt;				// pt resolution graphs
	TGraph *ggrd0;				// D resolution graphs
	TGraph *ggrz0;				// z0 resolution graphs
	TGraph *ggrth;				// theta resolution
	// Setup graph arrays
	Npt = 200;			// Nr. of points per graph
	Double_t * pt1 = new Double_t[Npt];
	Double_t * pp1 = new Double_t[Npt];
	Double_t *spt1 = new Double_t[Npt];
	Double_t *sd01 = new Double_t[Npt];
	Double_t *sz01 = new Double_t[Npt];
	Double_t *sth1 = new Double_t[Npt];
	// Fill graph arrays
	//Double_t ptmin = 1.0;
	//Double_t ptmax = 100;
	//Double_t pts = (ptmax - ptmin) / (Double_t)(Npt - 1);
	//
	SolGridCov *GC = new SolGridCov();
	GC->Read("CovIDEA_GT.root");
	for (Int_t k = 0; k < Npt; k++)	// Loop on pt
	{
		Double_t x[3]; Double_t p[3];
		x[0] = 0; x[1] = 0; x[2] = 0;			// Set origin
		pt1[k] = ptmin + k*pts;					// Set transverse momentum
		p[0] = pt1[k]; p[1] = 0; p[2] = pt1[k] / TMath::Tan(th);
		pp1[k] = pt1[k] / TMath::Sin(th);			// Set momentum
		//
		TMatrixDSym Cv = GC->GetCov(pt1[k], Ang);
		Double_t dptopt = 2 * TMath::Sqrt(Cv(2,2))*pt1[k] / (0.2998*G->B());
		spt1[k] = dptopt;					// Dpt/pt
		sd01[k] = TMath::Sqrt(Cv(0,0))*1e6;			// D  res. - change to microns
		sz01[k] = TMath::Sqrt(Cv(3, 3))*1e6;			// z0 res. - change to microns
		Double_t sint = TMath::Sin(th); Double_t sint2 = sint*sint;
		sth1[k] = TMath::Sqrt(Cv(4, 4)) *sint2;	// theta resolution
	}
	//
	// Compare pt resolution
	resol1->cd(1);
	ggrpt = new TGraph(Npt, pt1, spt1);			// Estimated resolution
	ggrpt->SetLineColor(kRed);
	ggrpt->SetTitle("#sigma_{pt}/pt");
	grept->SetTitle("#sigma_{pt}/pt");
	grept->SetMinimum(0.0);
	ggrpt->SetMinimum(0.0);
	ggrpt->GetXaxis()->SetTitle("pt (GeV)");
	grept->Draw("AP");							// Simulated resolution
	ggrpt->Draw("SAME");							// Estimated resolution
	// Compare d0 resolution
	resol1->cd(2);
	ggrd0 = new TGraph(Npt, pp1, sd01);			// Estimated resolution
	ggrd0->SetLineColor(kRed);
	ggrd0->SetTitle("D_{0} (#mum)");
	gred0->SetTitle("D_{0} (#mum)");
	gred0->SetMinimum(0.0);
	ggrd0->SetMinimum(0.0);
	ggrd0->GetXaxis()->SetTitle("p (GeV)");
	gred0->Draw("AP");							// Simulated resolution
	ggrd0->Draw("SAME");							// Estimated resolution
	// Compare z0 resolution
	resol1->cd(3);
	ggrz0 = new TGraph(Npt, pp1, sz01);			// Estimated resolution
	ggrz0->SetLineColor(kRed);
	ggrz0->SetTitle("Z_{0} (#mum)");
	grez0->SetTitle("Z_{0} (#mum)");
	grez0->SetTitle("Z_{0} (#mum)");
	grez0->SetMinimum(0.0);
	ggrz0->GetXaxis()->SetTitle("p (GeV)");
	grez0->Draw("AP");							// Simulated resolution
	ggrz0->SetMarkerColor(kRed);
	ggrz0->SetMarkerSize(0.5);
	ggrz0->SetMarkerStyle(kFullCircle);
	ggrz0->Draw("LSAME");						// Estimated resolution
	// Compare theta resolution
	resol1->cd(4);
	ggrth = new TGraph(Npt, pp1, sth1);			// Estimated resolution
	ggrth->SetLineColor(kRed);
	ggrth->SetTitle("#theta (rad)");
	greth->SetTitle("#theta (rad)");
	greth->SetMinimum(0.0);
	ggrth->SetMinimum(0.0);
	ggrth->GetXaxis()->SetTitle("p (GeV)");
	greth->Draw("AP");							// Simulated resolution
	ggrth->SetMarkerColor(kRed);
	grth->SetMarkerSize(0.5);
	ggrth->SetMarkerStyle(kFullCircle);
	ggrth->Draw("LSAME");						// Estimated resolution
	resol1->SaveAs("IDEA_resolution_tracking_grid.pdf");
}

