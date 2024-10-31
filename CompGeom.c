#include <TMath.h>
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
#include <TSystem.h>
#include <iostream>
#include "geometry_scripts/SolGeom.h"
#include "trkcovariance_scripts/SolTrack.h"
#include "trkcovariance_scripts/SolGridCov.h"
//
void CompGeom(Double_t Ang)
{

	Int_t iang = TMath::Nint(Ang); // Angle to integer

	//
	// Track to draw
	Double_t x[3] = { 0.0, 0.0, 0.0 };		// Track starting point
	Double_t ppt = 5.0;						// Track pt
	Double_t th = Ang * TMath::Pi() / 180.;
	Double_t ppz = ppt / TMath::Tan(th);	// Track pz
	Double_t p[3] = { ppt, 0.0, ppz };		// Track momentum


	//
	// ***********************
	// Jim: magnetic fields!!!
	// ***********************
	double B_IDEA = 2.0;		// IDEA magnetic field
	double B_CLD = 2.0;			// CLD magnetic field
	double B_SiD = 5.0;			// SiD magnetic field



	//
	//**************************
	//	Initialize geometry    *
	//**************************
	//
	//
	SolGeom *Gidea;			// Initialize IDEA geometry
	Gidea = new SolGeom("geometry_files/GeoIDEA_BASE.txt",B_IDEA);	// Geometry IDEA
	Gidea->Draw();			// Draw IDEA geometry
	TCanvas *cc_id = Gidea->cnv();					// Get canvas with geo display
	cc_id->cd(1);
	SolTrack *trk_id = new SolTrack(x, p, Gidea);	// Initialize track
	TGraph *gr_id = trk_id->TrkPlot();			// graph intersection with layers
	gr_id->Draw("PLSAME");						// plot track
	TCanvas *cnv_id = new TCanvas("cnv_id", "IDEA Geometry");
	cnv_id->Divide(1, 1); cnv_id->cd(1); cc_id->DrawClonePad();
	//cnv_id->SaveAs("IDEA_geometry_track.pdf");
	TString filename_IDEA_track;
	filename_IDEA_track.Form("plot_dump/IDEA_geometry_track_%d.pdf", iang);
	cnv_id->SaveAs(filename_IDEA_track);



	//
	SolGeom *Gcld;			// Initialize CLD  geometry	
	Gcld  = new SolGeom("geometry_files/GeoCLD.txt",B_CLD);			// Geometry CLD
	Gcld->Draw();			// Draw CLD  geometry
	TCanvas *cc_cl = Gcld->cnv();
	cc_cl->cd(1);
	SolTrack *trk_cl = new SolTrack(x, p, Gcld);	// Initialize track
	TGraph *gr_cl = trk_cl->TrkPlot();			// graph intersection with layers
	gr_cl->Draw("PLSAME");						// plot track
	TCanvas *cnv_cl = new TCanvas("cnv_cl", "CLD Geometry");
	cnv_cl->Divide(1, 1); cnv_cl->cd(1); cc_cl->DrawClonePad();
	//cnv_cl->SaveAs("CLD_geometry_track.pdf");
	TString filename_CLD_track;
	filename_CLD_track.Form("plot_dump/CLD_geometry_track_%d.pdf", iang);
	cnv_cl->SaveAs(filename_CLD_track);
	cc_cl->Close(); gSystem->ProcessEvents(); delete cc_cl; cc_cl = 0;
	SolGeom *Gsid;			// Initialize SiD  geometry	
	Gsid  = new SolGeom("geometry_files/GeoSiD.txt",B_SiD);			// Geometry SiD
	Gsid->Draw();			// Draw SiD  geometry
	TCanvas *cc_sid = Gsid->cnv();
	cc_sid->cd(1);
	SolTrack *trk_sid = new SolTrack(x, p, Gsid);	// Initialize track
	TGraph *gr_sid = trk_sid->TrkPlot();			// graph intersection with layers
	gr_sid->Draw("PLSAME");						// plot track
	TCanvas *cnv_sid = new TCanvas("cnv_sid", "SiD Geometry");
	cnv_sid->Divide(1, 1); cnv_sid->cd(1); cc_sid->DrawClonePad();
	//cnv_sid->SaveAs("SiD_geometry_track.pdf");
	TString filename_SiD_track;
	filename_SiD_track.Form("plot_dump/SiD_geometry_track_%d.pdf", iang);
	cnv_sid->SaveAs(filename_SiD_track);
	cnv_sid->Close(); gSystem->ProcessEvents(); delete cc_sid; cc_sid = 0;

	//
	// IDEA-BASE without wrapper
	//
	SolGeom *GideaNW;			// Initialize IDEA geometry withe no wrapper
	const Int_t nDet = 9;
	Bool_t OK[nDet] = {		// Enable selected parts of the detector 
		1,					// Beam pipe
		1,					// Inner VTX pixel layers
		1,					// Outer VTX layers
		1,					// Drift chamber
		0,					// Barrel Si wrapper
		0,					// Barrel pre-shower
		1,					// Forw. VTX pixel layers
		0,					// Forw. Si wrapper
		0 };				// Forw. pre-shower
	GideaNW = new SolGeom(OK);
	//
	//char* fname = "GeoIDEA_BASE.txt";	// Write out IDEA geometry
	//Gidea->GeoPrint(fname);
	//
	// Geometry initialized and drawn
	//********************************************************
	//	
	//******************************************************
	// Compare track parameter resolutions vs pt and theta *
	//******************************************************
	//
	TCanvas *resol = new TCanvas("resol", "Comparison of resolutions", 100, 100, 500, 500);
	TString CnvTitle; 
	CnvTitle.Form("Comparison of resolutions - Track angle %d deg.",iang);
	resol->SetTitle(CnvTitle);
	
	resol->Divide(2, 2);
	// Define graphs for IDEA
	TGraph *grpt_id;				// pt resolution graphs
	TGraph *grptms_id;				// pt resolution graphs MS only
	TGraph *grd0_id;				// D resolution graphs
	TGraph *grz0_id;				// z0 resolution graphs
	TGraph *grth_id;				// theta resolution
	// Define graphs for CLD
	TGraph *grpt_cl;				// pt resolution graphs
	TGraph *grptms_cl;				// pt resolution graphs MS only
	TGraph *grd0_cl;				// D resolution graphs
	TGraph *grz0_cl;				// z0 resolution graphs
	TGraph *grth_cl;				// theta resolution
	// Define graphs for SiD
	TGraph *grpt_sid;				// pt resolution graphs
	TGraph *grptms_sid;				// pt resolution graphs MS only
	TGraph *grd0_sid;				// D resolution graphs
	TGraph *grz0_sid;				// z0 resolution graphs
	TGraph *grth_sid;				// theta resolution
	// Setup graph arrays
	Int_t Npt = 200;			// Nr. of points per graph
	Double_t * pt = new Double_t[Npt];
	Double_t * pp = new Double_t[Npt];
	Double_t *spt_id = new Double_t[Npt];
	Double_t *sptms_id = new Double_t[Npt];
	Double_t *spt_idnw = new Double_t[Npt];
	Double_t *sd0_id = new Double_t[Npt];
	Double_t *sz0_id = new Double_t[Npt];
	Double_t *sth_id = new Double_t[Npt];
	Double_t *spt_cl = new Double_t[Npt];
	Double_t *sptms_cl = new Double_t[Npt];
	Double_t *sd0_cl = new Double_t[Npt];
	Double_t *sz0_cl = new Double_t[Npt];
	Double_t *sth_cl = new Double_t[Npt];
	Double_t *spt_sid = new Double_t[Npt];
	Double_t *sptms_sid = new Double_t[Npt];
	Double_t *sd0_sid = new Double_t[Npt];
	Double_t *sz0_sid = new Double_t[Npt];
	Double_t *sth_sid = new Double_t[Npt];
	// Fill graph arrays
	Double_t ptmin = 2.0;
	Double_t ptmax = 100;
	Double_t pts = (ptmax - ptmin) / (Double_t)(Npt-1);
	for (Int_t k = 0; k < Npt; k++)	// Loop on pt
	{
		Double_t x[3]; Double_t p[3];
		x[0] = 0; x[1] = 0; x[2] = 0;			// Set origin
		pt[k] = ptmin+k*pts;					// Set transverse momentum
		p[0] = pt[k]; p[1] = 0;	p[2] = pt[k] / TMath::Tan(th);
		pp[k] = pt[k]/TMath::Sin(th);			// Set momentum
		
		
		// Fill IDEA arrays
		SolTrack *tr_id = new SolTrack(x, p, Gidea);	// Initialize track
		Bool_t Res = kTRUE;		// Enable detector resolution effects
		Bool_t MS  = kTRUE;		// Enable multiple scattering
		tr_id->CovCalc(Res,MS);					// Calculate covariance
		spt_id[k] = tr_id->s_pt();							// Dpt/pt
		//spt_id[k] = tr_id->s_pt()/pt[k];							// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2
		sd0_id[k] = tr_id->s_D()*1e6;							// D  res. - change to microns
		sz0_id[k] = tr_id->s_z0()*1e6;						// z0 res. - change to microns
		sth_id[k] = tr_id->s_ct() / (1 + pow(tr_id->ct(), 2));	// theta resolution
		//
		SolTrack *tr_idnw = new SolTrack(x, p, GideaNW);	// Initialize track
		Res = kTRUE;
		MS  = kTRUE;
		tr_idnw->CovCalc(Res, MS);					// Calculate covariance
		spt_idnw[k] = tr_idnw->s_pt();							// Dpt/pt
		//spt_idnw[k] = tr_idnw->s_pt()/pt[k];							// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2
		//
		Res = kFALSE;
		MS = kTRUE;
		tr_id->CovCalc(Res,MS);					// Calculate covariance with only MS
		sptms_id[k] = tr_id->s_pt();							// Dpt/pt
		//sptms_id[k] = tr_id->s_pt()/pt[k];							// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2

		
		// Fill CLD arrays
		SolTrack *tr_cl = new SolTrack(x, p, Gcld);	// Initialize track
		Res = kTRUE;		// Enable detector resolution effects
		MS  = kTRUE;		// Enable multiple scattering
		tr_cl->CovCalc(Res,MS);					// Calculate covariance
		spt_cl[k] = tr_cl->s_pt();							// Dpt/pt
		//spt_cl[k] = tr_cl->s_pt()/pt[k];					// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2
		sd0_cl[k] = tr_cl->s_D()*1e6;							// D  res. - change to microns
		sz0_cl[k] = tr_cl->s_z0()*1e6;						// z0 res. - change to microns
		sth_cl[k] = tr_cl->s_ct() / (1 + pow(tr_cl->ct(), 2));	// theta resolution
		Res = kFALSE;
		tr_cl->CovCalc(Res,MS);					// Calculate covariance with only MS
		sptms_cl[k] = tr_cl->s_pt();							// Dpt/pt
		//sptms_cl[k] = tr_cl->s_pt()/pt[k];					// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2


		// Fill SiD arrays
		SolTrack *tr_sid = new SolTrack(x, p, Gsid);	// Initialize track
		Res = kTRUE;		// Enable detector resolution effects
		MS  = kTRUE;		// Enable multiple scattering
		tr_sid->CovCalc(Res,MS);					// Calculate covariance
		spt_sid[k] = tr_sid->s_pt();							// Dpt/pt
		//spt_sid[k] = tr_sid->s_pt()/pt[k];					// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2
		sd0_sid[k] = tr_sid->s_D()*1e6;							// D  res. - change to microns
		sz0_sid[k] = tr_sid->s_z0()*1e6;						// z0 res. - change to microns
		sth_sid[k] = tr_sid->s_ct() / (1 + pow(tr_sid->ct(), 2));	// theta resolution
		Res = kFALSE;
		tr_sid->CovCalc(Res,MS);					// Calculate covariance with only MS
		sptms_sid[k] = tr_sid->s_pt();							// Dpt/pt
		//sptms_sid[k] = tr_sid->s_pt()/pt[k];					// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2

	}
	//
	// Compare pt resolution
	resol->cd(1);
	// CLD
	grpt_cl = new TGraph(Npt, pt, spt_cl);			// pt resolution
	grpt_cl->SetLineColor(kRed);
	grpt_cl->SetMarkerColor(kRed);
	grpt_cl->SetTitle("#sigma_{pt}/pt");
	//grpt_cl->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grpt_cl->SetMinimum(0.0);
	grpt_cl->GetXaxis()->SetTitle("pt (GeV)");
	grpt_cl->Draw("APL");
	grptms_cl = new TGraph(Npt, pt, sptms_cl);			// pt resolution MS only
	grptms_cl->SetLineColor(kRed);
	grptms_cl->SetMarkerColor(kRed);
	grptms_cl->SetLineStyle(7);
	grptms_cl->SetTitle("#sigma_{pt}/pt");
	//grptms_cl->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grptms_cl->SetMinimum(0.0);
	grptms_cl->GetXaxis()->SetTitle("pt (GeV)");
	grptms_cl->Draw("SAME");
	// IDEA
	grpt_id = new TGraph(Npt, pt, spt_id);			// pt resolution
	grpt_id->SetLineColor(kBlue);
	grpt_id->SetMarkerColor(kBlue);
	grpt_id->SetTitle("#sigma_{pt}/pt");
	//grpt_id->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grpt_id->SetMinimum(0.0);
	grpt_id->GetXaxis()->SetTitle("pt (GeV)");
	grpt_id->Draw("SAME");
	grptms_id = new TGraph(Npt, pt, sptms_id);			// pt resolution MS only
	grptms_id->SetLineColor(kBlue);
	grptms_id->SetMarkerColor(kBlue);
	grptms_id->SetLineStyle(7);
	grptms_id->SetTitle("#sigma_{pt}/pt");
	//grptms_id->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grptms_id->SetMinimum(0.0);
	grptms_id->GetXaxis()->SetTitle("pt (GeV)");
	grptms_id->Draw("SAME");
	//SiD
	grpt_sid = new TGraph(Npt, pt, spt_sid);			// pt resolution
	grpt_sid->SetLineColor(kGreen);
	grpt_sid->SetMarkerColor(kGreen);
	grpt_sid->SetTitle("#sigma_{pt}/pt");
	//grpt_sid->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grpt_sid->SetMinimum(0.0);
	grpt_sid->GetXaxis()->SetTitle("pt (GeV)");
	grpt_sid->Draw("SAME");
	grptms_sid = new TGraph(Npt, pt, sptms_sid);			// pt resolution MS only
	grptms_sid->SetLineColor(kGreen);
	grptms_sid->SetMarkerColor(kGreen);
	grptms_sid->SetLineStyle(7);
	grptms_sid->SetTitle("#sigma_{pt}/pt");
	//grptms_sid->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grptms_sid->SetMinimum(0.0);
	grptms_sid->GetXaxis()->SetTitle("pt (GeV)");
	grptms_sid->Draw("SAME");

	
	TLegend *lgpt = new TLegend(0.2, 0.9, 0.6, 0.70);
	TString LgTitle; 
	LgTitle.Form("Track angle %d deg.",iang);
	lgpt->SetHeader(LgTitle);
	lgpt->AddEntry(grpt_id, "IDEA", "L");
	lgpt->AddEntry(grptms_id, "IDEA MS only", "L");
	lgpt->AddEntry(grpt_cl, "CLD", "L");
	lgpt->AddEntry(grptms_cl, "CLD MS only", "L");
	lgpt->AddEntry(grpt_sid, "SiD", "L");
	lgpt->AddEntry(grptms_sid, "SiD MS only", "L");
	lgpt->Draw();

	// Compare d0 resolution
	resol->cd(2);
	grd0_id = new TGraph(Npt, pp, sd0_id);			// D resolution
	grd0_id->SetLineColor(kBlue);
	grd0_id->SetMarkerColor(kBlue);
	grd0_id->SetTitle("#sigma_{D_{0}} (#mum)");
	grd0_id->SetMinimum(0.0);
	grd0_id->GetXaxis()->SetTitle("p (GeV)");
	grd0_id->Draw("APL");
	grd0_cl = new TGraph(Npt, pp, sd0_cl);			// D resolution
	grd0_cl->SetLineColor(kRed);
	grd0_cl->SetMarkerColor(kRed);
	grd0_cl->SetTitle("#sigma_{D_{0}} (#mum)");
	grd0_cl->SetMinimum(0.0);
	grd0_cl->GetXaxis()->SetTitle("p (GeV)");
	grd0_cl->Draw("SAME"); 
	grd0_sid = new TGraph(Npt, pp, sd0_sid);			// D resolution
	grd0_sid->SetLineColor(kGreen);
	grd0_sid->SetMarkerColor(kGreen);
	grd0_sid->SetTitle("#sigma_{D_{0}} (#mum)");
	grd0_sid->SetMinimum(0.0);
	grd0_sid->GetXaxis()->SetTitle("p (GeV)");
	grd0_sid->Draw("SAME");
	TLegend *lgd0 = new TLegend(0.2, 0.9, 0.6, 0.70);
	lgd0->SetHeader(LgTitle);
	lgd0->AddEntry(grpt_id, "IDEA", "L");
	lgd0->AddEntry(grpt_cl, "CLD", "L");
	lgd0->AddEntry(grpt_sid, "SiD", "L");
	lgd0->Draw();


	// Compare z0 resolution
	resol->cd(3);
	//IDEA
	grz0_id = new TGraph(Npt, pp, sz0_id);			// z0 resolution
	grz0_id->SetLineColor(kBlue);
	grz0_id->SetMarkerColor(kBlue);
	grz0_id->SetTitle("#sigma_{Z_{0}} (#mum)");
	grz0_id->GetXaxis()->SetTitle("p (GeV)");
	grz0_id->Draw("APL");
	//CLD
	grz0_cl = new TGraph(Npt, pp, sz0_cl);			// z0 resolution
	grz0_cl->SetLineColor(kRed);
	grz0_cl->SetMarkerColor(kRed);
	grz0_cl->SetTitle("#sigma_{Z_{0}} (#mum)");
	grz0_cl->GetXaxis()->SetTitle("p (GeV)");
	grz0_cl->Draw("SAME");			
	//SiD
	grz0_sid = new TGraph(Npt, pp, sz0_sid);			// z0 resolution
	grz0_sid->SetLineColor(kGreen);
	grz0_sid->SetMarkerColor(kGreen);
	grz0_sid->SetTitle("#sigma_{Z_{0}} (#mum)");
	grz0_sid->GetXaxis()->SetTitle("p (GeV)");
	grz0_sid->Draw("SAME");

	TLegend *lgz0 = new TLegend(0.2, 0.9, 0.6, 0.70);
	lgz0->SetHeader(LgTitle);
	lgz0->AddEntry(grpt_id, "IDEA", "L");
	lgz0->AddEntry(grpt_cl, "CLD", "L");
	lgz0->AddEntry(grpt_sid, "SiD", "L");
	lgz0->Draw();
	resol->cd(4);

	// Compare theta resolution
	//IDEA
	grth_id = new TGraph(Npt, pp, sth_id);			// theta resolution
	grth_id->SetLineColor(kBlue);
	grth_id->SetMarkerColor(kBlue);
	grth_id->SetTitle("#sigma_{#theta} (rad)");
	grth_id->SetMinimum(0.0);
	grth_id->GetXaxis()->SetTitle("p (GeV)");
	grth_id->Draw("APL");
	//CLD
	grth_cl = new TGraph(Npt, pp, sth_cl);			// theta resolution
	grth_cl->SetLineColor(kRed);
	grth_cl->SetMarkerColor(kRed);
	grth_cl->SetTitle("#sigma_{#theta} (rad)");
	grth_cl->SetMinimum(0.0);
	grth_cl->GetXaxis()->SetTitle("p (GeV)");
	grth_cl->Draw("SAME");
	//SiD
	grth_sid = new TGraph(Npt, pp, sth_sid);			// theta resolution
	grth_sid->SetLineColor(kGreen);
	grth_sid->SetMarkerColor(kGreen);
	grth_sid->SetTitle("#sigma_{#theta} (rad)");
	grth_sid->SetMinimum(0.0);
	grth_sid->GetXaxis()->SetTitle("p (GeV)");
	grth_sid->Draw("SAME");

	TLegend *lgth = new TLegend(0.2, 0.9, 0.6, 0.70);
	lgth->SetHeader(LgTitle);
	lgth->AddEntry(grpt_id, "IDEA", "L");
	lgth->AddEntry(grpt_cl, "CLD", "L");
	lgth->AddEntry(grpt_sid, "SiD", "L");
	lgth->Draw();

	TString filename;
	filename.Form("plot_dump/IDEA_CLD_SiD_resolution_tracking_%d.pdf", iang);
	resol->SaveAs(filename);
	//


	TCanvas *resolp = new TCanvas("resolp", "Comparison of pt resolutions", 50, 50, 500, 500);
	resolp->Divide(1, 1);
	// Compare pt resolution
	// CLD
	resolp->cd(1); 
	grpt_cl->SetMaximum(0.005);
	grpt_cl->Draw("APL");
	grptms_cl->Draw("SAME");
	// SiD
	grpt_sid->Draw("SAME");
	// IDEA
	TGraph *grpt_idnw = new TGraph(Npt, pt, spt_idnw);			// pt resolution
	grpt_idnw->SetLineColor(kBlue);
	grpt_idnw->SetLineStyle(2);
	grpt_idnw->SetMarkerColor(kBlue);
	grpt_idnw->SetTitle("#sigma_{pt}/pt");
	grpt_idnw->SetMinimum(0.0);
	grpt_idnw->GetXaxis()->SetTitle("pt (GeV)");
	grpt_idnw->Draw("SAME");
	grpt_id->Draw("SAME");
	grptms_id->Draw("SAME");


	TLegend *lgpt1 = new TLegend(0.2, 0.9, 0.6, 0.70);
	lgpt1->SetHeader(LgTitle);
	lgpt1->AddEntry(grpt_id, "IDEA", "L");
	lgpt1->AddEntry(grptms_id, "IDEA", "L");
	lgpt1->AddEntry(grpt_idnw, "IDEA No Si wrapper", "L");
	lgpt1->AddEntry(grpt_cl, "CLD", "L");
	lgpt1->AddEntry(grptms_cl, "CLD MS only", "L");
	lgpt1->AddEntry(grpt_sid, "SiD", "L");
	lgpt1->Draw();

	TString filenamep;
	filenamep.Form("plot_dump/IDEA_CLD_SiD_resolution_tracking_pt_%d.pdf", iang);
	resolp->SaveAs(filenamep);



	/*
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
	GC->Read("CovCLD.root");
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
	*/
}

