#include <TMath.h>
#include <TGraph.h>
#include <iostream>
#include <TCanvas.h>
#include <TPave.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TF1.h>
#include <TString.h>
#include "SolGeom.h"
#include "SolTrack.h"
#include <TLegend.h>

SolGeom::SolGeom()
{
	SolGeoInit();
	for (Int_t i = 0; i < fNdet; i++)fEnable[i] = kTRUE;	// default is everything enabled
	SolGeoFill();
	TString OldLab = " "; Int_t k = 0;
	for (Int_t i = 0; i < fNlay; i++)
	{
		if (fLyLabl[i] != OldLab)
		{
			fDtype[k] = fLyLabl[i];
			fDfstLay[k] = i;
			OldLab = fLyLabl[i];
			k++;
			if (k > fNdty)
			{
				cout << "SolGeom::SolGeom : Too many detector types! Layer "<<i<<" reached" << endl;
				return;
			}
		}
	}
	for (Int_t i=0; i < k; i++)cout << "i = " << i << ", Detector = " << fDtype[i]<<endl;
}
//
SolGeom::SolGeom(Bool_t *OK)
{
	SolGeoInit();
	for (Int_t i = 0; i < fNdet; i++)fEnable[i] = OK[i];	// User defined list
	SolGeoFill();
	TString OldLab = " "; Int_t k = 0;
	for (Int_t i = 0; i < fNlay; i++)
	{
		if (fLyLabl[i] != OldLab)
		{
			fDtype[k] = fLyLabl[i];
			fDfstLay[k] = i;
			OldLab = fLyLabl[i];
			k++;
			if (k > fNdty)
			{
				cout << "SolGeom::SolGeom : Too many detector types! Layer " << i << " reached" << endl;
				return;
			}
		}
	}
	for (Int_t i=0; i < k; i++)cout << "i = " << i << ", Detector = " << fDtype[i]<<endl;
}
SolGeom::SolGeom(char *fname)
{
	SolGeoInit();
	for (Int_t i = 0; i < fNdet; i++)fEnable[i] = kTRUE;	// default is everything enabled
	GeoRead(fname);
	TString OldLab = " "; Int_t k = 0;
	for (Int_t i = 0; i < fNlay; i++)
	{
		if (fLyLabl[i] != OldLab)
		{
			fDtype[k] = fLyLabl[i];
			fDfstLay[k] = i;
			OldLab = fLyLabl[i];
			k++;
			if (k > fNdty)
			{
				cout << "SolGeom::SolGeom : Too many detector types! Layer " << i << " reached" << endl;
				return;
			}
		}
	}
	for (Int_t i = 0; i < k; i++)cout << "i = " << i << ", Detector = " << fDtype[i] << endl;
}

void SolGeom::SolGeoInit()
{
	//
	// Magnetic field
	//
	fB = 5.0;
	//
	// Create arrays
	//
	ftyLay = new Int_t[fNlMax];		// Layer type 1 = R (barrel) or 2 = z (forward/backward)
	fLyLabl = new TString[fNlMax];	// Layer label
	fxMin = new Double_t[fNlMax];	// Minimum dimension z for barrel  or R for forward
	fxMax = new Double_t[fNlMax];	// Maximum dimension z for barrel  or R for forward
	frPos = new Double_t[fNlMax];	// R/z location of layer
	fthLay = new Double_t[fNlMax];	// Thickness (meters)
	frlLay = new Double_t[fNlMax];	// Radiation length (meters)
	fnmLay = new Int_t[fNlMax];		// Number of measurements in layers (1D or 2D)
	fstLayU = new Double_t[fNlMax];	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
	fstLayL = new Double_t[fNlMax];	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
	fsgLayU = new Double_t[fNlMax];	// Resolution Upper side (meters) - 0 = no measurement
	fsgLayL = new Double_t[fNlMax];	// Resolution Lower side (meters) - 0 = no measurement
	fflLay = new Bool_t[fNlMax];	// measurement flag = T, scattering only = F
	fEnable = new Bool_t[fNdet];	// list of enabled detectors
	fDtype = new TString[fNdty];	// Array with layer labels 
	fDfstLay = new Int_t[fNdty];	// Array with start layer
	//
	// Load geometry info in SolGeom.h
	//
	fNlay = 0;	// Actual number of layers
	fBlay = 0;	// Nr. of barrel layers
	fFlay = 0;	// Nr. of forward/backward layers
	fNm = 0;	// Nr. of measuring layers
}
	//
void SolGeom::SolGeoFill()
{
	//===================================================================================
	//		BARREL REGION
	//===================================================================================
	//
	Double_t R12 = TMath::Sqrt(12);
	//
	// Beam pipe
	//
	if (fEnable[0])
	{
		ftyLay[fNlay] = 1;			// Layer type 1 = R (barrel) or 2 = z (forward/backward)
		fLyLabl[fNlay] = "PIPE";
		fxMin[fNlay] = -100.;		// Minimum dimension z for barrel  or R for forward
		fxMax[fNlay] = 100.;		// Maximum dimension z for barrel  or R for forward
		frPos[fNlay] = 0.012;		// R/z location of layer (in meters)
		fthLay[fNlay] = 0.0004;		// Thickness (meters)
		frlLay[fNlay] = 35.276e-2;	// Radiation length (meters) //Jim: for Be - taken from https://pdg.lbl.gov/2022/AtomicNuclearProperties/HTML/beryllium_Be.html
		fnmLay[fNlay] = 0;			// Number of measurements in layers (1D or 2D)
		fstLayU[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
		fstLayL[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
		fsgLayU[fNlay] = 0.;			// Resolution Upper side (meters) - 0 = no measurement
		fsgLayL[fNlay] = 0.;			// Resolution Lower side (meters) - 0 = no measurement
		fflLay[fNlay] = kFALSE;		// measurement flag = T, scattering only = F
		fNlay++; fBlay++;
	}




	//
	// Vertex  detector (inner)
	if (fEnable[1])
	{
		const Int_t NlVtx = 5;	// 5 vertex barrel pixel layers
		Double_t rVtx[NlVtx] = { 1.4, 2.2, 3.5, 4.8, 6.0 };		// Vertex layer radii in cm
		Double_t lVtx[NlVtx] = { 6.3,6.3,6.3,6.3,6.3 };		// Vertex layer half length in cm
		for (Int_t i = 0; i < NlVtx; i++)
		{
			ftyLay[fNlay] = 1;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "VTX";			// Layer label
			fxMin[fNlay] = -lVtx[i] * 1.e-2;		// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = lVtx[i] * 1.e-2;	// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = rVtx[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = 50.E-6;			// Thickness (meters) - Jim: assume 50 microns thickness (ARCADIA/ATLASPIX3)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters) for Si: https://pdg.lbl.gov/2023/AtomicNuclearProperties/HTML/silicon_Si.html
			fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0;					// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
			fstLayL[fNlay] = TMath::Pi() / 2.;		// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
			fsgLayU[fNlay] = 3.E-6;				// Resolution Upper side (meters) - 0 = no measurement - see TDR 2.2.1: "All of these technologies have the capability of delivering sensors [..] with 5 µm hit resolution"
			fsgLayL[fNlay] = 3.E-6;				// Resolution Lower side (meters) - 0 = no measurement - see TDR 2.2.1: "All of these technologies have the capability of delivering sensors [..] with 5 µm hit resolution"
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			fNlay++; fBlay++;
			fNm++;
		}


		// Describe associated material -- Jim: asume 0.1% X0 per layer, i.e. 93.7microns per layer
		const Int_t NlVtxM = 5;

		Double_t rVtxM[NlVtxM] = { 1.4+0.0050, 2.2+0.0050, 3.5+0.0050, 4.8+0.0050, 6.0+0.0050 };		// Vertex layer radii in cm -- Jim: assume r for vtx sensitive material above + sensor width of 50microns
		Double_t lVtxM[NlVtxM] = { 6.3,6.3,6.3,6.3,6.3 };		// Vertex layer half length in cm
		Double_t lThkM[NlVtxM] = { 43.7,43.7,43.7,43.7,43.7 };	// Layer thickness in um -- Jim: 93.7 - 50 = 43.7 microns for Si
		for (Int_t i = 0; i < NlVtxM; i++)
		{
			ftyLay[fNlay] = 1;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "VTX";				// Layer label
			fxMin[fNlay] = -lVtxM[i] * 1.e-2;	// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = lVtxM[i] * 1.e-2;	// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = rVtxM[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = lThkM[i] * 1.E-6;	// Thickness (meters)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 0;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0;					// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
			fstLayL[fNlay] = 0;	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
			fsgLayU[fNlay] = 0;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 0;				// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kFALSE;				// measurement flag = T, scattering only = F
			fNlay++; fBlay++;
		}




	}


	//
	// Tracker
	if (fEnable[2])
	{
		const Int_t NlTrki = 5;	// Assume 5 long pixel layers
		Double_t rTrki[NlTrki] = { 21.95,46.95,71.95,96.95,121.95};		// Tracker layer radii in cm
		Double_t lTrki[NlTrki] = { 111.6/2.0,147.3/2.0,200.1/2.0,251.8/2.0,304.5/2.0 };	// Tracker layer half length in cm -- Jim: careful, ILC TDR gives full length, not half length!
		for (Int_t i = 0; i < NlTrki; i++)
		{
			ftyLay[fNlay] = 1;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "TRK";				// Layer label
			fxMin[fNlay] = -lTrki[i] * 1.e-2;	// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = lTrki[i] * 1.e-2;	// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = rTrki[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = 100.E-6;			// Thickness (meters)-- Jim: assume 100microns
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters) -- Jim: Si
			fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0 = axial layer - Upper side
			fstLayL[fNlay] = TMath::Pi() / 2.;	// Stereo angle (rad) - pi/2 = z layer - Lower side
			fsgLayU[fNlay] = 7.E-6;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 7.E-6;			// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			fNlay++; fBlay++;
			fNm++;
		}
		// Describe associated material -- Jim: asume 0.3% X0 per layer, i.e. 281 microns per layer
		const Int_t NlTrkiM = 5;	// Assume 5 material layers
		Double_t rTrkiM[NlTrkiM] = { 21.95+0.01,46.95+0.01,71.95+0.01,96.95+0.01,121.95+0.01 };	// Tracker layer radii in cm -- Jim: assume r for trk sensitive material above + sensor width of 100microns
		Double_t lTrkiM[NlTrkiM] = { 111.6/2.0,147.3/2.0,200.1/2.0,251.8/2.0,304.5/2.0 };	// Tracker layer half length in cm
		Double_t lThkiM[NlTrkiM] = { 181., 181.,181., 181., 181.};	// Layer thickness in um of Si -- Jim: 281 - 100 = 181 microns for Si
		for (Int_t i = 0; i < NlTrkiM; i++)
		{
			ftyLay[fNlay] = 1;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "TRK";				// Layer label
			fxMin[fNlay] = -lTrkiM[i] * 1.e-2;	// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = lTrkiM[i] * 1.e-2;	// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = rTrkiM[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = lThkiM[i] * 1.E-6;	// Thickness (meters)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 0;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0 = axial layer - Upper side
			fstLayL[fNlay] = 0;				// Stereo angle (rad) - pi/2 = z layer - Lower side
			fsgLayU[fNlay] = 0;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 0;				// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kFALSE;				// measurement flag = T, scattering only = F
			fNlay++; fBlay++;
		}
	}
	
	//
	//================================================================================================
	//		FORWARD/BACKWARD
	//================================================================================================
	//
	// Vertex disks
	if (fEnable[3])
	{
		const Int_t NlVtxd = 8;	// 4 vertex disks on each side
		Double_t zVtxd[NlVtxd] = { -17.2, -12.3, -9.2, -7.2, 
			                        7.2, 9.2, 12.3, 17.2};		// Vertex layer z in cm
		Double_t riVtxd[NlVtxd] = { 2.0,1.8, 1.6, 1.4, 
			                        1.4, 1.6, 1.8, 2.0};	// Vertex layer R min in cm
		Double_t rοVtxd[NlVtxd] = { 7.1,7.1,7.1,7.1, 
			                        7.1,7.1,7.1,7.1 };	// Vertex layer R max in cm

		for (Int_t i = 0; i < NlVtxd; i++)
		{
			ftyLay[fNlay] = 2;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "VTXDSK";			// Layer label
			fxMin[fNlay] = riVtxd[i] * 1.e-2;	// Minimum dimension R for forward disk
			fxMax[fNlay] = rοVtxd[i] * 1.e-2;		// Maximum dimension R for forward disk
 			frPos[fNlay] = zVtxd[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = 50.E-6;				// Thickness (meters) - Jim: assume 50 microns thickness (ARCADIA/ATLASPIX3)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0 = axial layer - Upper side
			fstLayL[fNlay] = TMath::Pi() / 2.;	// Stereo angle (rad) - pi/2 = z layer - Lower side
			fsgLayU[fNlay] = 3.E-6;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 3.E-6;			// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			fNlay++; fFlay++;
			fNm++;
		}
		// Describe associated material -- Jim: asume 0.1% X0 per layer, i.e. 93.7microns per layer
		const Int_t NlVtxdM = 8;	// 4 vertex disks on each side
		Double_t zVtxdM[NlVtxdM] =  { -17.2-0.005, -12.3-0.005, -9.2-0.005, -7.2-0.005, 
			                           7.2+0.005, 9.2+0.005, 12.3+0.005, 17.2+0.005};	// Material layer z in cm -- Jim: same as z for sensitive material above + 50microns
		Double_t riVtxdM[NlVtxdM] = { 2.0,1.8, 1.6, 1.4, 
			                        1.4, 1.6, 1.8, 2.0};	// Material layer R min in cm
		Double_t rοVtxdM[NlVtxdM] = { 7.1,7.1,7.1,7.1, 
			                        7.1,7.1,7.1,7.1 };	// Material layer R max in cm
		Double_t lThVdM[NlVtxdM] =  { 43.7,43.7,43.7,43.7,
									43.7,43.7,43.7,43.7 };		// Layer thickness in um -- Jim: 93.7 - 50 = 43.7 microns for Si
		for (Int_t i = 0; i < NlVtxdM; i++)
		{
			ftyLay[fNlay] = 2;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "VTXDSK";			// Layer label
			fxMin[fNlay] = riVtxdM[i] * 1.e-2;	// Minimum dimension R for forward disk
			fxMax[fNlay] = rοVtxdM[i]* 1.e-2;		// Maximum dimension R for forward disk
			frPos[fNlay] = zVtxdM[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = lThVdM[i] * 1.E-6;	// Thickness (meters)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 0;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0 = axial layer - Upper side
			fstLayL[fNlay] = 0;				// Stereo angle (rad) - pi/2 = z layer - Lower side
			fsgLayU[fNlay] = 0;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 0;				// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kFALSE;				// measurement flag = T, scattering only = F
			fNlay++; fFlay++;
		}
	}


	// Vertex disks - forward
	if (fEnable[4])
	{
		const Int_t NlVtxfwdd = 6;	// 3 forward vertex disks on each side
		Double_t zVtxfwdd[NlVtxfwdd] = {-83.2, -54.1, -20.7, 
			                             20.7, 54.1, 83.2};		// Vertex layer z in cm
		Double_t riVtxfwdd[NlVtxfwdd] = { 11.7, 7.6, 2.8, 
			                              2.8, 7.6,11.7 };	// Vertex layer R min in cm
		Double_t rοVtxfwdd[NlVtxfwdd] = { 16.6,16.6,16.6, 
			                              16.6,16.6,16.6 };	// Vertex layer R max in cm
		for (Int_t i = 0; i < NlVtxfwdd; i++)
		{
			ftyLay[fNlay] = 2;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "VTXDSK_FWD";			// Layer label
			fxMin[fNlay] = riVtxfwdd[i] * 1.e-2;	// Minimum dimension R for forward disk
			fxMax[fNlay] = rοVtxfwdd[i] * 1.e-2;		// Maximum dimension R for forward disk
 			frPos[fNlay] = zVtxfwdd[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = 50.E-6;				// Thickness (meters) - Jim: assume 50 microns thickness (ARCADIA/ATLASPIX3)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0 = axial layer - Upper side
			fstLayL[fNlay] = TMath::Pi() / 2.;	// Stereo angle (rad) - pi/2 = z layer - Lower side
			fsgLayU[fNlay] = 3.E-6;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 3.E-6;			// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			fNlay++; fFlay++;
			fNm++;
		}
		// Describe associated material -- Jim: asume 0.1% X0 per layer, i.e. 93.7microns per layer
		const Int_t NlVtxdfwdM = 6;	// 3 forward vertex disks on each side
		Double_t zVtxdfwdM[NlVtxdfwdM] = {-83.2-0.005, -54.1-0.005, -20.7-0.005, 
			                             20.7+0.005, 54.1+0.005, 83.2+0.005};	// Material layer z in cm -- Jim: same as z for sensitive material above + 50microns
		Double_t riVtxdfwdM[NlVtxdfwdM] = { 11.7, 7.6, 2.8, 
			                              2.8, 7.6,11.7 };	// Material layer R min in cm
		Double_t rοVtxdfwdM[NlVtxdfwdM] = { 16.6,16.6,16.6, 
			                              16.6,16.6,16.6 };	// Material layer R max in cm
		Double_t lThVdfwdM[NlVtxdfwdM] =  { 43.7,43.7,43.7,
									        43.7,43.7,43.7 };		// Layer thickness in um -- Jim: 93.7 - 50 = 43.7 microns for Si
		for (Int_t i = 0; i < NlVtxdfwdM; i++)
		{
			ftyLay[fNlay] = 2;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "VTXDSK_FWD";			// Layer label
			fxMin[fNlay] = riVtxdfwdM[i] * 1.e-2;	// Minimum dimension R for forward disk
			fxMax[fNlay] = rοVtxdfwdM[i] * 1.e-2;		// Maximum dimension R for forward disk
			frPos[fNlay] = zVtxdfwdM[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = lThVdfwdM[i] * 1.E-6;	// Thickness (meters)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 0;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0 = axial layer - Upper side
			fstLayL[fNlay] = 0;				// Stereo angle (rad) - pi/2 = z layer - Lower side
			fsgLayU[fNlay] = 0;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 0;				// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kFALSE;				// measurement flag = T, scattering only = F
			fNlay++; fFlay++;
		}
	}



	//
	// Tracker disks	
	//
	if (fEnable[5])
	{
		const Int_t NlItkd = 8;	// 8 tracker disks on each side
		Double_t zItkd[NlItkd] = { -164.09, -135.55, -107.50, -78.89,
									78.89,107.50, 135.55, 164.09 };		// Vertex layer z in cm
		Double_t riItkd[NlItkd] = { 20.89, 20.89, 20.89, 20.89,
									20.89, 20.89, 20.89, 20.89  };	// Vertex layer R min in cm
		Double_t roItkd[NlItkd] = { 125.36, 100.31, 75.14, 49.80,
									49.80, 75.14, 100.31, 125.36 };		// Vertex layer R max in cm
		for (Int_t i = 0; i < NlItkd; i++)
		{
			ftyLay[fNlay] = 2;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "TRKDSK";			// Layer label
			fxMin[fNlay] = riItkd[i] * 1.e-2;	// Minimum dimension R for forward disk
			fxMax[fNlay] = roItkd[i] * 1.e-2;	// Maximum dimension R for forward disk
			frPos[fNlay] = zItkd[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] =100.E-6;			// Thickness (meters) -- Jim: assume 100microns
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0 = axial layer - Upper side
			fstLayL[fNlay] = TMath::Pi() / 2.;	// Stereo angle (rad) - pi/2 = z layer - Lower side
			fsgLayU[fNlay] = 7.E-6;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 7.E-6;			// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			fNlay++; fFlay++;
			fNm++;
		}
		// Describe associated material -- Jim: asume 0.3% X0 per layer, i.e. 281 microns per layer
		const Int_t NlItkdM = 8;	// 8 tracker disks on each side
		Double_t zItkdM[NlItkdM] = { -164.09-0.01, -135.55-0.01, -107.50-0.01, -78.89-0.01,
									78.89+0.01,107.50+0.01, 135.55+0.01, 164.09+0.01 };			// Material layer z in cm -- Jim: same as z for sensitive material above + 100microns
		Double_t riItkdM[NlItkdM] = { 20.89, 20.89, 20.89, 20.89,
									20.89, 20.89, 20.89, 20.89  };		// Vertex layer R min in cm	
		Double_t roItkdM[NlItkdM] = { 125.36, 100.31, 75.14, 49.80,
									49.80, 75.14, 100.31, 125.36 };		// Vertex layer R max in cm
		Double_t lThItdM[NlItkdM] = { 181., 181., 181., 181.,
									  181., 181., 181., 181. };		// Layer thickness in um of Si -- Jim: 281 - 100 = 181 microns for Si
		for (Int_t i = 0; i < NlItkdM; i++)
		{
			ftyLay[fNlay] = 2;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "TRKDSK";			// Layer label
			fxMin[fNlay] = riItkdM[i] * 1.e-2;	// Minimum dimension R for forward disk
			fxMax[fNlay] = roItkdM[i] * 1.e-2;		// Maximum dimension R for forward disk
			frPos[fNlay] = zItkdM[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = lThItdM[i] * 1.E-6;	// Thickness (meters)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 0;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0 = axial layer - Upper side
			fstLayL[fNlay] = 0;				// Stereo angle (rad) - pi/2 = z layer - Lower side
			fsgLayU[fNlay] = 0;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 0;				// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kFALSE;				// measurement flag = T, scattering only = F
			fNlay++; fFlay++;
		}
	}
	//
	//
	//
	// Magnet
	//
	ftyLay[fNlay] = 1;			// Layer type 1 = R (barrel) or 2 = z (forward/backward)
	fLyLabl[fNlay] = "MAG";		// Layer label
	fxMin[fNlay] = -2.793;		// Minimum dimension z for barrel  or R for forward -- Jim: ILC TDR Table II-6.1
	fxMax[fNlay] = 2.793;			// Maximum dimension z for barrel  or R for forward  -- Jim: ILC TDR Table II-6.1
	frPos[fNlay] = 2.731;		// R/z location of layer  -- Jim: ILC TDR Table II-6.1
	fthLay[fNlay] = 0.38;		// Thickness (meters)  -- Jim: ILC TDR Table II-6.1 3.112-2.731 = 0.38m
	frlLay[fNlay] = 6.58e-2;	// Radiation length (meters)
	fnmLay[fNlay] = 0;			// Number of measurements in layers (1D or 2D)
	fstLayU[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
	fstLayL[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
	fsgLayU[fNlay] = 0.;		// Resolution Upper side (meters) - 0 = no measurement
	fsgLayL[fNlay] = 0.;		// Resolution Lower side (meters) - 0 = no measurement
	fflLay[fNlay] = kFALSE;		// measurement flag = T, scattering only = F
	fNlay++; fBlay++;

	cout << "Geometry created with " << fNlay << "/" << fNm << " layers" << endl;
}
//
// Print geometry
void SolGeom::GeoPrint(char *fname)
{
	FILE *fdata = fopen(fname, "w");
	if (!fdata)
	{
		cout << "SolGeom::GeoPrint - can't open output file" << endl;
		return;
	}
	for (Int_t l = 0; l < fNlay; l++)
	{
		fprintf(fdata, "%d %s %g %g %g %g %g %d %g %g %g %g %d\n",
		ftyLay[l], fLyLabl[l].Data(), fxMin[l], fxMax[l], frPos[l], fthLay[l],
		frlLay[l], fnmLay[l], fstLayU[l], fstLayL[l], fsgLayU[l], fsgLayL[l], fflLay[l]);
		//cout << strng << endl<< endl;
	}	
	fclose(fdata);
}
//
// Material counter (Units are fraction of X0)
Double_t *SolGeom::FracX0(Double_t theta)
{
	//
	// Calculates amount of material crossed by a straight track at a polar angle theta
	// for each subdetector: Jim
	// 0: Pipe, 1: VTX, 2: TRK, 3: VTXDSK, 4: VTXDSK_FWD, 5: TRKDSK, 6: MAG,
	//
	Double_t *Mat;
	Mat = new Double_t[fNdty];
	for (Int_t i = 0; i < fNdty; i++)Mat[i] = 0;
	if (fNlay <= 0)
	{
		cout << "SolGeom::FracX0 : No geometry available. # layers = " << fNlay << endl;
		return Mat;
	}
	//
	// Loop over all layers
	Double_t lmb = 0.0;
	if (TMath::Abs(theta - TMath::PiOver2()) > 1.0e-10 && 
		TMath::Abs(theta)  > 1.0e-10) lmb = 1.0 / TMath::Tan(theta);	// Cot(theta)
	if (theta == 0.0) lmb = 1.e10;
	for (Int_t il = 0; il < fNlay; il++)
	{
		Int_t dNum;
		for (Int_t i = 0; i<fNdty; i++) if (fLyLabl[il] == fDtype[i])dNum = i;
		cout << "dnum = " << dNum << ", detector: "<<fDtype[dNum]<<endl;
		if (ftyLay[il] == 1)		// Cylinder at constant R
		{
			Double_t R = frPos[il];
			Double_t z = lmb*R;
			cout << "l num: " << il << ", R = " << R << ", z min: "<<fxMin[il]<<", z = " << z<<" , z max: "<<fxMax[il] << endl;
			if (z>fxMin[il] && z < fxMax[il])	// the layer is hit
			{
				Mat[dNum] += fthLay[il] / (TMath::Sin(theta)*frlLay[il]);
			}
		}
		else if (ftyLay[il] == 2) // disk at constant z
		{
			Double_t z = frPos[il];
			Double_t R = z / lmb;
			cout << "l num: " << il << ", z = " << z << ", R min: " << fxMin[il] << ", R = " << R << " , R max: " << fxMax[il] << endl;
			if (R>fxMin[il] && R < fxMax[il])	// the layer is hit
			{
				Mat[dNum] += fthLay[il] / (TMath::Cos(theta)*frlLay[il]);
			}
		}
	}
	//
	return Mat;
}
//
// Read geometry
void SolGeom::GeoRead(char *fname)
{
	char strng[200];
	int nbytes = 200;
	FILE *fdata = fopen(fname, "r");
	if (!fdata)
	{
		cout << "SolGeom::GeoRead - can't open input file" << endl;
		return;
	}
	Int_t tyLay;
	char LyLabl[20];
	float xMin;
	float xMax;
	float rPos;
	float thLay;
	float rlLay;
	Int_t nmLay;
	float stLayU;
	float stLayL;
	float sgLayU;
	float sgLayL;
	Int_t flLay;
	//
	while (fgets(strng, nbytes, fdata) != NULL)
	{
		cout << strng;
		int status = sscanf(strng, "%d %s %g %g %g %g %g %d %g %g %g %g %d",
			&tyLay, LyLabl, &xMin, &xMax, &rPos, &thLay,
			&rlLay, &nmLay, &stLayU, &stLayL, &sgLayU, &sgLayL, &flLay);
		ftyLay[fNlay] = tyLay; 
		fLyLabl[fNlay] = LyLabl; 
		fxMin[fNlay] = (Double_t) xMin; 
		fxMax[fNlay] = (Double_t) xMax; 
		frPos[fNlay] = (Double_t) rPos; 
		fthLay[fNlay] = (Double_t) thLay; 
		frlLay[fNlay] = (Double_t) rlLay; 
		fnmLay[fNlay] = nmLay; 
		fstLayU[fNlay] = (Double_t) stLayU; 
		fstLayL[fNlay] = (Double_t) stLayL; 
		fsgLayU[fNlay] = (Double_t) sgLayU; 
		fsgLayL[fNlay] = (Double_t) sgLayL; 
		fflLay[fNlay] = (Bool_t) flLay; 
		//cout << "Layer # " << fNlay << ": " << fLyLabl[fNlay] << ", Position: " << frPos[fNlay]
		//	<< ", Measurement: " << fflLay[fNlay] << endl;
		
		fNlay++;
		if (tyLay == 1)fBlay++;
		if (flLay == 1)fNm++;

	}
	fclose(fdata);
	cout << "SolGeom::GeoRead completed with " << fNlay << " layers input" << endl;
}
//
// Destructor
SolGeom::~SolGeom()
{
	fNlay = 0;
	fBlay = 0;
	fNm = 0;

	delete[] & ftyLay;
	delete[] & fxMin;
	delete[] & fxMax;
	delete[] & frPos;
	delete[] & fthLay;
	delete[] & frlLay;
	delete[] & fnmLay;
	delete[] & fstLayU;
	delete[] & fstLayL;
	delete[] & fsgLayU;
	delete[] & fsgLayL;
	delete[] & fflLay;
	delete[] & fEnable;
}
//
// Draw the geometry (just a sketch)
//
void SolGeom::Draw()
{
	//Double_t zMin = -2.9; Double_t zMax = 2.9;
	//Double_t rMax = 3.3;

	Double_t zMin = -3.5; Double_t zMax = 3.5;
	Double_t rMax = 3.5;


	fcnv = new TCanvas("cnv", "Geometry sketch", 10, 10, 950, 550);
	fcnv->Range(zMin, -1.0, zMax, rMax);

	// 
	// beam pipe
	//if (fEnable[0])
	//{
		TPave *pipe = new TPave(zMin, -frPos[0], zMax, frPos[0], 0, "");
		pipe->SetFillColor(kYellow);
		pipe->Draw();
	//}
	// Beamline
	TLine *beam = new TLine(zMin, 0.0, zMax, 0.0);
	beam->SetLineColor(kBlack);
	beam->SetLineWidth(1);
	beam->SetLineStyle(9);
	beam->Draw("SAME");
	// Magnet
	TPave *sol = new TPave(-2.793, 2.731, 2.793, 3.112, 0, ""); //Jim: change to SiD solenoid dimensions from https://pages.uoregon.edu/silicondetector/sid-dimensions.html
	sol->SetFillColor(30);
	sol->Draw("SAME");
	//
	// Draw Calorimeter
	// ECAL Barrel
	const Int_t nP = 4;
	Double_t brECALX[nP] = { -1.765, 1.765, 1.765, -1.765};
	Double_t brECALY[nP] = {  1.265, 1.265, 1.409,  1.409};
	TPolyLine *brECALor = new TPolyLine(nP, brECALX, brECALY,"F");
	brECALor->SetFillColor(38);
	brECALor->SetLineColor(kBlack);
	brECALor->Draw("FSAME");

	// HCAL Barrel
	Double_t brHCALX[nP] = { -3.018, 3.018, 3.018, -3.018};
	Double_t brHCALY[nP] = {  1.417, 1.417, 2.493,  2.493};
	TPolyLine *brHCALor = new TPolyLine(nP, brHCALX, brHCALY,"F");
	brHCALor->SetFillColor(39);
	brHCALor->SetLineColor(kBlack);
	brHCALor->Draw("FSAME");


	// ECAL Backward
	Double_t bkECALX[nP] = { -1.657, -1.800, -1.800, -1.657 };
	Double_t bkECALY[nP] = { 0.2195, 0.2195, 1.250, 1.250};
	TPolyLine *bkECALor = new TPolyLine(nP, bkECALX, bkECALY, "F");
	bkECALor->SetFillColor(38);
	bkECALor->SetLineColor(kBlack);
	bkECALor->Draw("FSAME");

	// ECAL Forward
	Double_t bfECALX[nP] = { 1.657, 1.800, 1.800, 1.657 };
	Double_t bfECALY[nP] = { 0.2195, 0.2195, 1.250, 1.250};
	TPolyLine *bfECALor = new TPolyLine(nP, bfECALX, bfECALY, "F");
	bfECALor->SetFillColor(38);
	bfECALor->SetLineColor(kBlack);
	bfECALor->Draw("FSAME");



	// HCAL Backward
	Double_t bkHCALX[nP] = { -1.805, -3.028, -3.028, -1.815 };
	Double_t bkHCALY[nP] = { 0.2195, 0.2195, 1.402, 1.402};
	TPolyLine *bkHCALor = new TPolyLine(nP, bkHCALX, bkHCALY, "F");
	bkHCALor->SetFillColor(39);
	bkHCALor->SetLineColor(kBlack);
	bkHCALor->Draw("FSAME");

	// HCAL Forward
	Double_t bfHCALX[nP] = { 1.805, 3.028, 3.028, 1.815 };
	Double_t bfHCALY[nP] = { 0.2195, 0.2195, 1.402, 1.402};
	TPolyLine *bfHCALor = new TPolyLine(nP, bfHCALX, bfHCALY, "F");
	bfHCALor->SetFillColor(39);
	bfHCALor->SetLineColor(kBlack);
	bfHCALor->Draw("FSAME");





	// All other layers
	// Measurement silicon (red), blue (DCH), scattering black
	//
	const Int_t lMax = 200;
	TLine *ln[lMax];
	TF1   *fn[lMax];
	Int_t il = 0; 
	Int_t ig = 0;
	for (Int_t i = 0; i < fNlay; i++)
	{
		if (fLyLabl[i] == "DCH")		// Drift chamber layers (hypeboloids)
		{
				char lab[10]; 
				Int_t stat;
				stat = sprintf(lab, "fun%d", ig);
				fn[ig] = new TF1(lab, this, &SolGeom::StereoHyp, lxMin(i), lxMax(i), 3, "SolGeom","StereoHyp");
				fn[ig]->SetParameter(0, lPos(i));
				fn[ig]->SetParameter(1, lStU(i));
				fn[ig]->SetParameter(2, (Double_t) i);
				fn[ig]->SetLineColor(kBlue);
				fn[ig]->Draw("SAME");
				ig++;
		}
		else
		{
			if(ftyLay[i] == 1)ln[il] = new TLine(lxMin(i), lPos(i), lxMax(i), lPos(i));
			else ln[il] = new TLine(lPos(i), lxMin(i), lPos(i), lxMax(i));
			ln[il]->SetLineColor(kBlack);
			if (isMeasure(i))
			{
				ln[il]->SetLineColor(kRed);
				ln[il]->SetLineWidth(2);
			}
			ln[il]->Draw("SAME");
			il++;
		}

	}
	//
	// Forward/backward geometry
	//
	for (Int_t i = fBlay; i < fNlay; i++)
	{
		ln[il] = new TLine(lPos(i), lxMin(i), lPos(i), lxMax(i));
		ln[il]->SetLineColor(kBlack);
			if (isMeasure(i))
			{
				ln[il]->SetLineColor(kRed);
				ln[il]->SetLineWidth(2);
			}
		ln[il]->Draw("SAME");
		il++;
	}


	// Legend
	TLegend *legend = new TLegend(0.1, 0.55, 0.3, 0.7);  // Adjust positioning as needed
	legend->SetTextSize(0.03);
	legend->SetBorderSize(0);
	legend->SetFillColorAlpha(kWhite, 0.0);  // Transparent background

	// Add entries to the legend
	legend->AddEntry(pipe, "Beam pipe", "f");
	legend->AddEntry(beam, "Beamline", "l");
	legend->AddEntry(sol, "Magnet", "f");
	legend->AddEntry(brECALor, "ECAL Barrel", "f");
	legend->AddEntry(brHCALor, "HCAL Barrel", "f");
	legend->Draw();

	// Add Axes
	TGaxis *axisZ = new TGaxis(zMin, 0, zMax, 0, zMin, zMax, 510, "");
	axisZ->SetTitle("z [m]");
	axisZ->SetTitleOffset(1.2);
	axisZ->SetLabelColor(kGray+2);
	axisZ->SetLineColorAlpha(kGray+2, 0.8);  // 50% transparent grey
	axisZ->SetLabelSize(0.03);
	axisZ->Draw("SAME");

	TGaxis *axisR = new TGaxis(0, 0, 0, rMax, 0, rMax, 510, "");
	axisR->SetTitle("r [m]");
	axisR->SetTitleOffset(1.3);
	axisR->SetLabelColor(kGray+2);
	axisR->SetLineColorAlpha(kGray+2, 0.8);  // 50% transparent grey
	axisR->SetLabelSize(0.03);
	// Suppress all labels initially
	axisR->SetLabelOffset(999);  // Move default labels off the canvas
	axisR->Draw("SAME");

	// Add custom labels at specific positions (e.g., z = 1, 2, and 3)
	TText *label1 = new TText(-0.05,1, "1");
	TText *label2 = new TText(-0.05,2, "2");
	TText *label3 = new TText(-0.05,3, "3");

	// Customize label appearance
	label1->SetTextColor(kGray+2);
	label2->SetTextColor(kGray+2);
	label3->SetTextColor(kGray+2);

	label1->SetTextSize(0.03);
	label2->SetTextSize(0.03);
	label3->SetTextSize(0.03);

	label1->SetTextAlign(22);  // Center alignment
	label2->SetTextAlign(22);
	label3->SetTextAlign(22);

	// Draw labels
	label1->Draw("SAME");
	label2->Draw("SAME");
	label3->Draw("SAME");






	fcnv->SaveAs("Geometry_SiD.png");
	fcnv->SaveAs("Geometry_SiD.pdf");

}
//
Double_t SolGeom::StereoHyp(Double_t *x, Double_t *p)
{
	Double_t R   = p[0];
	Double_t tg  = TMath::Tan(p[1]);
	Int_t i = (Int_t)p[2];
	Double_t r = TMath::Sqrt(R*R - lxMax(i)*lxMax(i)*tg*tg);
	Double_t z   = x[0];
	//
	return TMath::Sqrt(r*r + z*z*tg*tg);
}
