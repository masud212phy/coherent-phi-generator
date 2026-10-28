#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <random>  // place with other headers
#include <TROOT.h>
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TGFrame.h"
#include "TRootEmbeddedCanvas.h"


using namespace std;


const double Me = 0.00051;                      //Mass of electron in GeV/C^2
const double Mp = 0.938272;                    //Mass of proton in GeV/c^2
const double mp = 1.49;                        //magnetic moment of GeV-1
const double Mphi = 1.02;
const double Mk = 0.493677;                      //Mass of phi meson in GeV/c^2
const double Mh = 3.7274;                      //Mass of helium-4 nucleus




const double pi = TMath::Pi();                //value of pi
const double D2R = pi/180.0;
const double R2D = 180.0/pi;

//electron's energy scattering off proton as a target

double ese_p(double Q2P, double xp, double Ee){
	
	double M = Mp;
	return Ee - Q2P/(2*Mp*xp);    // Q^2 = -q^2
	
	
}




//electron's energy scattering off helium-4 nuclus as a target
double ese_H(double Q2H, double xh, double Ee){
	
	double M = Mp;
	return Ee - Q2H/(2*Mh*xh);    // Q^2 = -q^2
	
	
}

//electron's angle scattering off proton as a target
double esa_p(double Ee, double ese_p, double Q2P){
	
	return 2*asin(sqrt(Q2P/(4*Ee*ese_p)));
}






//electron's angle scattering off helium-4 nuclus as a target
double esa(double Ee, double ese, double Q2H){
	
	return 2*asin(sqrt(Q2H/(4*Ee*ese)));
}



//helum-4's scattered energy as a function of t channel variable 
double spe(double t_h, double Mh){
	//double Mp =M;
	return (2*Mh*Mh-t_h)/(2*Mh);
	
}




// Proton and helium form-factor functions



// Proton charge form factor: standard dipole
double protonChargeFormFactor(double tp){
	
    const double dipoleMass2 = 0.71;  // GeV^2

    return std::pow(1.0 - tp/dipoleMass2,-2.0);
	
}



// Helium-4 charge-like form factor with radius rescaling
//
// lambda = R_g/R_C
//
// F_lambda(t) = F_C(lambda^2 t)
//
double heliumFormFactor(double t, double lambda){
    double scaledT = lambda*lambda*t;

    return(1.0 - std::pow(2.5*scaledT, 6))*std::exp(11.7*scaledT);


}




// Smooth exponential amplitude with no zero
//
// The slope 11.7 GeV^-2 is chosen to retain the exponential
// envelope of the charge-like baseline while removing its node.
//
double heliumSmoothFormFactor(double t)
{
    return std::exp(11.7*t);
}



/*
// Proton exponential cross-section shape
double protonCrossSectionShape(double tp, double W){
    const double B0 = 2.2;         // GeV^-2
    const double alphaPrime = 0.24; // GeV^-2

    double B = B0 + 4.0*alphaPrime*std::log(W/1.0); // W measured in GeV

    return std::exp(B*tp);
}

*/



// Helium-4 cross-section shape for a chosen lambda
double heliumCrossSectionShape(double t, double lambda){
	
	
    const double A = 4.0;

    // Your current proton momentum-transfer mapping
    double tp = t/16.0;

    //double dsigmaProton = protonCrossSectionShape(tp, W);

    double FCp = protonChargeFormFactor(tp);

    double FHe = heliumFormFactor(t, lambda);

    double nuclearAmplitude = A*FHe/FCp;

    //return dsigmaProton*nuclearAmplitude*nuclearAmplitude;
	
	return nuclearAmplitude*nuclearAmplitude;
}




// Smooth exponential helium hypothesis
double heliumSmoothCrossSectionShape(double t)
	
	
{
    const double A = 4.0;

    double tp = t/16.0;

    //double dsigmaProton = protonCrossSectionShape(tp, W);

    double FCp = protonChargeFormFactor(tp);

    double FHeSmooth = heliumSmoothFormFactor(t);

    double nuclearAmplitude = A*FHeSmooth/FCp;

    //return dsigmaProton*nuclearAmplitude*nuclearAmplitude;
	return nuclearAmplitude*nuclearAmplitude;
}










void generator(int NEVENTS = 1000000, const char* outname = "dvmpH.dat", int seed = 1){
	
	
	std::mt19937 rng(seed);

	const double zCenter = -3.7315;               //cm
	const double zHalfLength = 22.5;               //cm
	
	std::uniform_real_distribution<double> zVertexDist(zCenter - zHalfLength, zCenter + zHalfLength);
	
    ofstream myfile;
    
   // myfile.open("dvmpH_5k.txt");
    //myfile.open("generator.txt");
    myfile.open(outname);

    if(!myfile.is_open()){
        cerr << "ERROR: could not open output file: " << outname << endl;
        return;
      }
	
	
	
    TRandom3 *rand3 = new TRandom3(seed);

    cout << "Generating " << NEVENTS << " events" << endl;
    cout << "Output LUND file = " << outname << endl;
    cout << "Seed = " << seed << endl;
	
	
	
	//double Ee =12.0;  
	//double Ee = 10.6;
	double Ee = 11.0;                           //beam energy in GeV
	double M = Mp;
	double m = Mphi;
	double Mh = 3.7274; 
	
	
    TH1F *Q2_P_generated = new TH1F("Q^{2}_P_generated for helium-4","Q^{2} distribution for Helium-4 as a target",100,0,15);
	TH1F *hQ2P_rec  = new TH1F("hQ2_rec",  "Q^{2}_{rec};Q^{2} [GeV^{2}];Counts", 100, 0.0, 15.0);
	TH1F *Q2_P_Weighted = new TH1F("Q^{2}_Pweighted for helium-4","Q^{2} distribution for Helium-4 as a target",100,0,15);
	
	
	
    TH1F *Q2_H_generated = new TH1F("Q^{2}_H_generated for helium-4","Q^{2} distribution for Helium-4 as a target",100,0,25);
	TH1F *Q2_H_Weighted = new TH1F("Q^{2}_Hweighted for helium-4","Q^{2} distribution for Helium-4 as a target",100,0,25);
	
	
	TH1F *wP_generated = new TH1F("wP_generated","W distribution for proton and ^{4}He;W [GeV];Counts",100, 1.5, 12.0);
	TH1F *wP_weighted = new TH1F("wP_weighted","W distribution for proton and ^{4}He;W [GeV];Weighted counts",100, 1.5, 12.0);

	TH1F *wH_generated = new TH1F("wH_generated","W distribution for proton and ^{4}He;W [GeV];Counts",100, 1.5, 12.0);
	TH1F *wH_weighted = new TH1F("wH_weighted","W distribution for proton and ^{4}He;W [GeV];Weighted counts",100, 1.5, 12.0);

	wP_weighted->Sumw2();
	wH_weighted->Sumw2();
	
	
	
	TH1F *hxp_generated = new TH1F("h_x_{B}","x_{B} distribution for Proton as a target",100,0,1);
	TH1F *hxh_generated = new TH1F("h_x_{B}","x_{B} distribution for Helium-4 as a target",100,0,1);
	
	TH1F *hxp_weighted = new TH1F("h_x_{B}","x_{B} distribution for Proton as a target",100,0,1);
	TH1F *hxh_weighted = new TH1F("h_x_{B}","x_{B} distribution for Helium-4 as a target",100,0,1);
	
	TH1F *htEff_generated = new TH1F("ht_{eff.}_{generated}","t_{eff}_generated distribution for Proton as a target",100,-20,10);
	TH1F *htEff_weighted = new TH1F("ht_{eff}_{weighted}","t_{eff}_weighted distribution for Helium-4 as a target",100,-20,10);
	
	TH1F *htH_generated = new TH1F("ht_{h.}_{generated}","t_{eff}_generated distribution for Proton as a target",100,-20,10);
	TH1F *htH_weighted = new TH1F("ht_{h}_{weighted}","t_{eff}_weighted distribution for Helium-4 as a target",100,-20,10);
	
	TH1F *tdiff_generated = new TH1F("tdiff_generated","t_{max}-t for ^{4}He;t_{max}-t [GeV^{2}];Counts", 100, -1, 2.0);
	TH1F *tdiff_weighted = new TH1F("tdiff_weighted","t_{max}-t for ^{4}He;t_{max}-t [GeV^{2}];Weighted counts", 100, -1, 2.0);
	
	
	//Phase-space in 2D histogram
	
	TH2F *Q2HW = new TH2F("W vs Q^{2}_{H}", "W vs. Q^{2}_{H}", 100,0,25,100,4,12);
	
	TH2F *Q2HXB = new TH2F("x_{B} vs Q^{2}_{H}", "x_{B} vs. Q^{2}_{H}", 100,0,0.8,100,0,30);
	
	//TH2F *Q2HtH = new TH2F("-t vs Q^{2}_{H}", "-t vs. Q^{2}_{H}", 100,0,10,100,0,30);
	
    TH2F *Q2HtH = new TH2F("-t vs Q^{2}_{H}", "-t vs. Q^{2}_{H}", 100,0,8,100,0,1);
	
	TH2F *Q2HW_rec  = new TH2F("Q2HW",  "Q^{2} vs W;Q^{2} [GeV^{2}];W [GeV]", 100, 0, 15, 100, 1.5, 5.0);
	TH2F *Q2HXB_rec = new TH2F("Q2HXB", "x_{B} vs Q^{2};x_{B};Q^{2} [GeV^{2}]", 100, 0, 1, 100, 0, 15);
	
	Q2HW_rec->Sumw2();
	Q2HXB_rec->Sumw2();
	
	
	//cross-section of sigma_L/sigma_R with Q^2
	
	TProfile *pRP = new TProfile("pRP","R vs Q^{2}_{H};Q^{2}_{P} [GeV^{2}];R = #sigma_{L}/#sigma_{T}",100, 0.8, 25.0);
	
	TProfile *pRH = new TProfile("pRH","R vs Q^{2}_{H};Q^{2}_{H} [GeV^{2}];R = #sigma_{L}/#sigma_{T}",100, 0.8, 25.0);
	
	

	//longitudinal-sigma_L vs Q^2 for measurement at different W: proton and helium-4 as target
	
	//TProfile *pSigmaT_Q2_W45_P = new TProfile("pSigmaT_Q2_W45_P","#sigma_{T} vs Q^{2} at W #approx 4.5 GeV;Q^{2} [GeV^{2}];#sigma_{T} [nb]",100, 0.8, 25.0);
	//TProfile *pSigmaT_Q2_W45_H = new TProfile("pSigmaT_Q2_W45_H","#sigma_{T} vs Q^{2}_{H} at W #approx 4.5 GeV;Q^{2}_{H} [GeV^{2}];#sigma_{T} [nb]",100, 0.8, 25.0);
	
	//TProfile *pSigmaT_Q2_W50_P = new TProfile( "pSigmaT_Q2_W50_P","#sigma_{T} vs Q^{2} at W #approx 5.0 GeV;Q^{2} [GeV^{2}];#sigma_{T} [nb]",100, 0.8, 25.0);
	//TProfile *pSigmaT_Q2_W50_H = new TProfile("pSigmaT_Q2_W50_H","#sigma_{T} vs Q^{2}_{H} at W #approx 5.0 GeV;Q^{2}_{H} [GeV^{2}];#sigma_{T} [nb]",100, 0.8, 25.0);
	
	//TProfile *pSigmaT_Q2_W75_P = new TProfile("pSigmaT_Q2_W75_P","#sigma_{T} vs Q^{2} at W #approx 7.5 GeV;Q^{2} [GeV^{2}];#sigma_{T} [nb]",100, 0.8, 25.0);
	//TProfile *pSigmaT_Q2_W75_H = new TProfile("pSigmaT_Q2_W75_H","#sigma_{T} vs Q^{2}_{H} at W #approx 7.5 GeV;Q^{2}_{H} [GeV^{2}];#sigma_{T} [nb]",100, 0.8, 25.0);
	
	TProfile *pSigmaT_Q2_P = new TProfile("pSigmaT_Q2_P",";Q^{2}_{H} [GeV^{2}];#sigma_{T} [nb]",100, 0.8, 25.0);
	
	
	TProfile *pFCsq_t_h = new TProfile("pFCsq_tdiff","F_{C,^{4}He}^{2} vs (t_{max}-t);t_{max}-t [GeV^{2}];F_{C,^{4}He}^{2}",100, 1e-3, 1.0);
	
	
	
	
	
	
	TH1F *ESE_P = new TH1F("E_{e'}","Distribution of E_{e'}",50,0,15);
	TH1F *ESE_H = new TH1F("E_{e'}","Distribution of E_{e'}",50,0,15);
	
	
	TH1F *THE = new TH1F("Theta for Electron","Distriubtion of #Theta_{e}",50,0,50);
	
	//t-distribution
	TH1F *TH = new TH1F("t for Helium", "Distribution of t_{H}",50,-50,5);
	TH1F *thmin = new TH1F("t for Helium", "Distribution of t_{H}_{min}",50,-50,5);
	TH1F *thmax = new TH1F("t for Helium", "Distribution of t_{H}_{max}",50,-50,5);
	
	
	
	TH1F *PH1 = new TH1F("Total Momentum for Helium","Distribution of P_{H}",100,0,1.0);
	TH1F *PT = new TH1F("Transverse Momentum for Helium","Distribution of P_{t}",100,0,1);
	TH2F *PZ = new TH2F("Longitudinal Momentum for Helium","Distribution of P_{z}",100,0,90,100,0,1);
	
	TH2F *PHT = new TH2F("Total Momentum for Helium","tmin-t vs. P_{H}",100,0,1,100,0,1.0);
	TH2F *PTT = new TH2F("Transverse Momentum for Helium","tmin-t vs. P_{t}",100,0,1,100,0,1);
	
	TH1F *TE = new TH1F("t_eff","Distribution of effective t", 100,-10,1.0);
	
	
    //CM frame analysis
	
	TH2F *PTheta = new TH2F("PvsTheta", "Distribution of Total Mom. vs. Theta", 100,-10,200,100,0,2.5);
		
		
	TH1F *P = new TH1F("P", "Distribution of P for scattered proton in CMS for scattered angle 90", 100,-5,5);
	
	TH1F *Theta = new TH1F("Theta", "Distribution of Theta for scattered proton in CMS ", 100,0,180);
	
	
	
	
    TH2F *ECMTheta = new TH2F("EcmvsTheta", "Distribution of Ecm vs. Theta", 100,0,50,100,0,10);
	TH2F *PTA = new TH2F("Angle vs. Pt", "Distribution of Ecm vs. Theta", 100,0,80,100,0,1.0);
	
	
	
	TH1F *ECM = new TH1F("EHCM", "Distribution of Total CM Energy of Initial State", 100,0,10);
	
	TH1F *EHCM = new TH1F("EHCM", "Distribution of Total CM Energy for scattered helium-4 nucleus", 100,0,10);
	
	TH1F *EHCM1 = new TH1F("EHCM", "Distribution of CM Energy for scattered helium-4 nucleus", 100,0,10);
	
	TH1F *ECMtot1 = new TH1F("Ecm1", "E_{cm} = E_{q}_{cm}+E_{H}_{cm}", 100,0,10);
	TH1F *ECMtot2 = new TH1F("Ecm2", "E_{cm} = E'_{H}_{cm}+E_{#phi}", 100,0,10);
	
	TH1F *SHCM = new TH1F("SHCM", "CM Energy for scattered helium-4 after boosting", 100,0,10);
	TH1F *SPCM = new TH1F("SPCM", "CM Energy for phi meson after boosting ", 100,0,10);
	
	TH1F *EHPhi = new TH1F("EHPhi", "CM Energy for phi meson after boosting ", 100,0,15);
	
	TH1F *PqCM = new TH1F("Pq_{cm}", "Before Collision: P_{q}_{cm} ", 100,0,11);
	TH1F *PHCM1 = new TH1F("P_{H}_{cm}", "Before Collision: P_{H}_{cm} ", 100,0,11);
	
	TH1F *PHPhi = new TH1F("PHPhi", "CM Momentum for phi meson after boosting ", 100,0,15);
	
	TH1F *PHCM2 = new TH1F("PH", "After Collision: P'_{H}_{cm} ", 100,0,10);
	TH1F *PPhiCM = new TH1F("PPhiCM", "After Collision: P_{#phi}", 100,0,10);
	
	
	
	//TH1F *Q2HMEq = new TH1F("E_{q} vs. Q^{2}_{H}", "ME_{q} vs. Q^{2}_{H}", 100,0,12,100,0,50);
	
	TH2F *ThetaPe = new TH2F("Theta_{e} vs. P_{e}", "Theta vs. Pe", 100,0,10,100,0,60);
	
	
	//TH1F *PHL = new TH1F("PHL", "LAB Momentum for Helium-4 nucleus ", 100,0.1,0.4);
	
	//form factors
	
	TH1F *CFFH = new TH1F("CFF for Helium","Distribution of CFFH",50,0,1);
	
	TH1F *CFFPt = new TH1F("CFFPt for effective proton","Distribution of CFFPt",50,0,1);
	
	
	//cross-section
	TH1F *R = new TH1F("R", "Distribution of R",50,0.8,100);
	
	
	TH1F *EKH = new TH1F("EKH", "Distribution of Kinetic Energy of Helium-4",50,0,0.5);
	TH1F *CTW = new TH1F("#sigma_{T}", "#sigma_{L} vs. W",50,1,100);
	TH1F *CTQ2 = new TH1F("#sigma_{T}", "#sigma_{L} vs. Q^{2}",50,1,100);
	TH1F *thtmin = new TH1F("t-t_{min}","Distribution of abs(t-t_{min})",100,0,-1.0);
	TH2F *Q2thtmin = new TH2F("Q^{2} vs. t-t_{min}","Distribution of abs(t-t_{min})",100,0.5,8.0,100,0,1.0);
	
	TH2F *ththeta = new TH2F("th vs. theta","Theta for scattered helium-4 in CM frame without cut with beamline axis",100,0.5,-20,100,0,180);
	TH2F *ththeta1 = new TH2F("th vs. theta1","Theta for scattered helium in lab frame without cut with beamline axis",100,0.0,-5.0,100,0,180);
	
	
	//recoil momentum and theta
	
	
	TH1F *pHe_generated = new TH1F("pHe_generated", "Recoil ^{4}He momentum;p_{He} [GeV/c];Counts", 100, 0.0, 10.0);

	TH1F *pHe_weighted = new TH1F("pHe_weighted","Recoil ^{4}He momentum;p_{He} [GeV/c];Weighted counts",100, 0.0, 10.0);

	TH1F *thetaHe_generated = new TH1F("thetaHe_generated","Recoil ^{4}He polar angle;#theta_{He} [degrees];Counts",100, 0.0, 180.0);

	TH1F *thetaHe_weighted = new TH1F("thetaHe_weighted","Recoil ^{4}He polar angle;#theta_{He} [degrees];Weighted counts",100, 0.0, 180.0);
	
	
	//2D histogram for tdiff vs theta_He and t_h vs theta_He
	
	TH2D *tdiff_vs_thetaHe = new TH2D("tdiff_vs_thetaHe","Generated t-t_{min} vs recoil ^{4}He angle;""#theta_{He}^{lab} [deg];t-t_{min} [GeV^{2}]",100, 0.0, 100.0,100, 0.0, 1.0);

	TH2D *tH_vs_thetaHe = new TH2D("tH_vs_thetaHe","Generated t_{h} vs recoil ^{4}He angle;""#theta_{He}^{lab} [deg];t_{h} [GeV^{2}]",100, 0.0, 100.0,500, -5.0, 0.0);
	
	
	// histogram for recoil p vs theta, and vz vs p or theta
	
	TH2D *h_p_vs_theta = new TH2D("h_p_vs_theta","Generated recoil ^{4}He;#theta_{He} [deg];p_{He} [GeV/c]",90, 0.0, 180.0,100, 0.0, 5.0);

	TH2D *h_theta_vs_vz = new TH2D("h_theta_vs_vz","Generated recoil ^{4}He;v_{z} [cm];#theta_{He} [deg]",100, -30.0, 25.0,90, 0.0, 180.0);

	TH2D *h_p_vs_vz = new TH2D("h_p_vs_vz","Generated recoil ^{4}He;v_{z} [cm];p_{He} [GeV/c]",100, -30.0, 25.0,100, 0.0, 1.0);
	
	
	
	TH1D *pKplus_generated = new TH1D("pKplus_generated","Generated K^{+} momentum;p_{K^{+}} [GeV/c];Counts",150, 0.0, 10.0);

	TH1D *pKminus_generated = new TH1D("pKminus_generated","Generated K^{-} momentum;p_{K^{-}} [GeV/c];Counts",150, 0.0, 10.0);

	TH1D *thetaKplus_generated = new TH1D("thetaKplus_generated","Generated K^{+} polar angle;#theta_{K^{+}} [deg];Counts",180, 0.0, 180.0);

	TH1D *thetaKminus_generated = new TH1D("thetaKminus_generated","Generated K^{-} polar angle;#theta_{K^{-}} [deg];Counts",180, 0.0, 180.0);
	
	
	TH2D *pVsThetaElec = new TH2D("pVsThetaElec","Generated K^{-};#theta_{e'} [deg];p_{e'} [GeV/c]",60, 0.0, 60.0,150, 0.0, 10.0);
	
	TH2D *pVsThetaKplus = new TH2D("pVsThetaKplus","Generated K^{+};#theta_{K^{+}} [deg];p_{K^{+}} [GeV/c]",180, 0.0, 180.0,150, 0.0, 10.0);

	TH2D *pVsThetaKminus = new TH2D("pVsThetaKminus","Generated K^{-};#theta_{K^{-}} [deg];p_{K^{-}} [GeV/c]",180, 0.0, 180.0,150, 0.0, 10.0);
	
	
	
	// Four distributions requested for the before/after t' comparison.
	TH1F *hMinusTBefore = new TH1F("hMinusTBefore", "-t before t' cut;-t [GeV^{2}];Events",200, 0.0, 10.0);
	TH1F *hTprimeBefore = new TH1F("hTprimeBefore","t'=t_{min}-t before cut;t' [GeV^{2}];Events",200, 0.0, 10.0);
	TH1F *hMinusTAfter = new TH1F("hMinusTAfter","-t after 0<t'<1.0 GeV^{2};-t [GeV^{2}];Events",200, 0.0, 10.0);
	TH1F *hTprimeAfter = new TH1F("hTprimeAfter","t'=t_{min}-t after 0<t'<1.0 GeV^{2};t' [GeV^{2}];Events",100, 0.0, 1.0);
	
	TH2F *Q2Hvst_h = new TH2F("Q2Pvst_h","Q2H and -t phase space; -t [GeV^{2}];Q^{2} [GeV^{2}]", 100, 0.0, 10.0, 100, 0.0, 25.0);
	
	
	// Define before the t' cut
	
	TH1D *h_minus_tmin = new TH1D("h_minus_tmin","Generated -t_{min} Before Cuts;" "-t_{min} [GeV^{2}];" "Events", 200, 0.0, 5.0);
	
	TH2D *h_tmin_vs_th = new TH2D( "h_tmin_vs_th","Generated -t_{min} vs -t_{h} Before Cuts;" "-t_{min} [GeV^{2}];" "-t_{h} [GeV^{2}]", 200, 0.0, 5.0, 400, 0.0, 10.0);
	
	TH2D *h_tmin_vs_recoil_theta = new TH2D( "h_tmin_vs_recoil_theta","Generated -t_{min} vs Recoil ^{4}He Angle; -t_{min} [GeV^{2}]; #theta_{^{4}He'}^{lab} [deg]",200, 0.0, 5.0,180, 0.0, 180.0);
	
	
	TH2D *h_t_vs_recoil_theta = new TH2D("h_t_vs_recoil_theta","Generated -t vs Recoil ^{4}He Angle;""-t [GeV^{2}];""#theta_{^{4}He'}^{lab} [deg]",200, 0.0, 5.0,90, 0.0, 90.0);
	
	
	TH2D *h_tdiff_vs_theta_after = new TH2D("h_tdiff_vs_theta_after", "Selected Events: t' vs Recoil ^{4}He Lab Angle;" "t' = t_{min}-t_{h} [GeV^{2}];" "#theta_{^{4}He'}^{lab} [deg]", 100, 0.0, 1.0, 90, 0.0, 90.0);
	
	TH2D *h_tmin_vs_th_after = new TH2D("h_tmin_vs_th_after", "Generated -t_{min,h} vs -t_{h} After t' Cut;" "-t_{min,h} [GeV^{2}];" "-t_{h} [GeV^{2}]", 200, 0.0, 5.0, 200, 0.0, 5.0);
	
	
	//invariant phi-mass check:
	
	TH1D *h_MKK = new TH1D( "h_MKK", "Generated K^{+}K^{-} Invariant Mass;" "M(K^{+}K^{-}) [GeV/c^{2}];" "Events",200, 0.98, 1.06);
	
	TH1D *h_deltaE = new TH1D("h_deltaE", "Generator Energy-Conservation Check;" "#Delta E = E_{initial}-E_{final} [GeV];" "Events", 200, -1.0e-6, 1.0e-6);
	
	TH1D *h_deltaP = new TH1D( "h_deltaP", "Generator Momentum-Conservation Check;" "|#Delta#vec{p}| [GeV/c];" "Events",200, 0.0, 0.1e-6);
	
	
	//TH1D *h_delta_t_recoil_phi = new TH1D("h_delta_t_recoil_phi", "Consistency of Momentum-Transfer Definitions;" "t_{recoil}-t_{#phi} [GeV^{2}];" "Events",200, -1e-6, 1.0e-6);
	
	TH1D *h_delta_t_recoil_phi = new TH1D("h_delta_t_recoil_phi","Consistency of Momentum-Transfer Definitions;" "t_{recoil}-t_{#phi} [GeV^{2}];" "Events",400, -10.0, 10.0);
	
	
	
	
	
	
	
	
	
	
	const int nTBins = 120;
	const double minusTLow  = 0.0;
	const double minusTHigh = 1.2;

	TH1D *hMinusTCompactModel = new TH1D("hMinusTCompactModel","Compact model;" "-t [GeV^{2}];" "d#sigma/dt [nb/GeV^{2}]", nTBins, minusTLow, minusTHigh);

	TH1D *hMinusTChargeModel = new TH1D( "hMinusTChargeModel","Charge-like model;" "-t [GeV^{2}];" "d#sigma/dt [nb/GeV^{2}]", nTBins, minusTLow, minusTHigh);

	TH1D *hMinusTBroadModel = new TH1D("hMinusTBroadModel", "Broad model;" "-t [GeV^{2}];" "d#sigma/dt [nb/GeV^{2}]", nTBins, minusTLow, minusTHigh);

	TH1D *hMinusTSmoothModel = new TH1D("hMinusTSmoothModel","Smooth exponential model;" "-t [GeV^{2}];" "d#sigma/dt [nb/GeV^{2}]",nTBins, minusTLow, minusTHigh);

	// Required for correctly propagated weighted uncertainties
	hMinusTCompactModel->Sumw2();
	hMinusTChargeModel->Sumw2();
	hMinusTBroadModel->Sumw2();
	hMinusTSmoothModel->Sumw2();
	
	
	
	
	
	
	
	
	
	long conservationDebugCount = 0;
	long conservationFailureCount = 0;

	const double conservationTolerance = 1.0e-8;  // GeV or GeV/c
	
	
	long tCheckCount = 0;
	long tFailureCount = 0;

	const double tTolerance =  1.0e-10;  // GeV^2
	
	
	TLorentzVector beam4v, scatelec4vP, scatelec4vH, q4v_P, q4v_H, target4vH, com4v_P, com4v_H, scat4vH, phi4v1, kaon4v1, kaon4v2,check4v1,check4v2,check4v3,target4vP,scat4vP,phi4vP,phi4vH;
	
	//for (int k = 0; k<NEVENTS; k++){
	
	
		long long nThrown = 0;
		long long nPassedTprime = 0;
		long long nSelected = 0;
		
		// Keep generating until exactly NEVENTS complete events are written.
		while (nSelected < static_cast<long long>(NEVENTS)) {
		++nThrown;
	
	//Smearing parameters
	 //const double sigmaE_rel = 0.005;   // 0.5% energy resolution
	 //const double sigmaTh    = 0.001;   // 1 mrad
	 //const double sigmaPhi   = 0.001;   // 1 mrad
	
	 
	double Q2P = rand3->Uniform(0.75,12.0);
	 
	 //double Q2P = 2.5;
	 Q2_P_generated->Fill(Q2P);
	
	 double Q2H = rand3->Uniform(0.81,25.52);
	 Q2_H_generated->Fill(Q2H);
	
	
    //if (Q2H < 0 || Q2H > 8) continue;
	
	//cout<<"Q2H:   "<<Q2H<<endl;
	
	//double Q2H = 10.0;
	
	double Eq_Pmin = (m*m+2*Mp*m+Q2P)/(2*Mp);
	double xminp = Q2P/(2*Mp*(Ee-0.5));
	double xmaxp = Q2P/(2*Mp*Eq_Pmin);
	double xp = rand3->Uniform(xminp,xmaxp);
	hxp_generated->Fill(xp);
	
	double Eq_Hmin = (m*m+2*Mh*m+Q2H)/(2*Mh);
	double xminh = Q2H/(2*Mh*(Ee-0.5));
    double xmaxh = Q2H/(2*Mh*Eq_Hmin);
	double xh = rand3->Uniform(xminh,xmaxh);
	hxh_generated->Fill(xh);
	//cout<<"xminh:    "<<xminh<<"xh:   "<<xh<<"xmaxh:  "<<xmaxh<<endl;
	
	
	double Wp = sqrt(Mp*Mp + Q2P/xp-Q2P);
	double Wh = sqrt(Mh*Mh + Q2H/xh-Q2H);
	
	// generated W distributions can go earlier
	wP_generated->Fill(Wp);
	wH_generated->Fill(Wh);
	
	
	
	beam4v.SetPxPyPzE(0,0,Ee,Ee);
	
	target4vP.SetPxPyPzE(0,0,0,Mp);
	target4vH.SetPxPyPzE(0,0,0,Mh);
	
	
	double eseP = ese_p(Q2P,xp,Ee);
	
    //double eseP_sm   = eseP * (1.0 + rand3->Gaus(0.0, sigmaE_rel));
	double eseH =ese_H(Q2H,xh,Ee);
	//cout<<"ese1:   "<<ese1<<endl;
	
	
	
	double esaP =esa_p(Ee,eseP,Q2P);
	//double esaP_sm = esaP*(1.0 + rand3->Gaus(0.0, sigmaTh));
	
	//double Q2P_rec = 4.0 * Ee * eseP_sm * pow(sin(esaP_sm/2.0), 2);
	//cout<<"Q2P_rec:   "<<Q2P_rec<<endl;
	//double nu_rec = Ee - eseP_sm;
	//double xB_rec = Q2P_rec / (2.0 * M * nu_rec);
	//cout<<"xB_rec:    "<<xB_rec<<endl;
	
	//double W_rec = M*M + Q2P_rec/xB_rec - Q2P_rec;
	//cout<<"Wrec:     "<<W_rec<<endl;
	
	
	//Helium-4
	double esaH =esa(Ee,eseH,Q2H);
	//cout<<"esa1:   "<<esa1<<endl;
	
	
	
	
	
	
	
    //if(esa1<=0.087 || esa1 >0.611 ) continue;                                                                                                                                                                    
         
	//cout<<"esa1:   "<<esa1*R2D<<endl;	 
		//double esa2 = esa1;
			                                                                                                                                                                                                                                                                                                                                                
        //THE->Fill(esa2*R2D); 
		
		//cout<<"esa2:   "<<esa1*R2D<<endl;                                                                                                                                                     
	
		//cout<<"ese1:   "<<ese1<<endl;
	//declaration of 4-vector of beam particle electron
		
	
	
	//double esa2 = esa1;
	





	//declaration of 4-vector of target as helium-4 nuclus
	
	//target4vH.SetPxPyPzE(0,0,0,Mh);
	
	
	
	//declaration of 4-vector of scattered electron 
	//scatelec4v.SetPxPyPzE(0,ese1*sin(esa1),ese1*cos(esa1),ese1);
	
	
	scatelec4vP.SetPxPyPzE(0,eseP*sin(esaP),eseP*cos(esaP),eseP);
	scatelec4vH.SetPxPyPzE(0,eseH*sin(esaH),eseH*cos(esaH),eseH);
	
	//double sephi = rand3->Uniform(0,2*pi);
	//scatelec4v.RotateZ(sephi);
	
	double Eq_P = Ee - eseP;
	q4v_P = beam4v - scatelec4vP;
	
	double Eq_H = Ee - eseH;
	q4v_H = beam4v - scatelec4vH;
	
	
	//double W = sqrt(Mh*Mh -Q2H*(1.0-1.0/xh));
	//double xB1 = 0.3;
	//double Q2H1 = 2*xB1*Mh*Eq;
	

	
	
	
	//TLorentzVector q4v1; 
	//q4v1 = beam4v - scatelec4v;
	
	//cout<<"q4vPxL:   "<<q4v.Px()<< "q4v1PxL1:   "<<q4v1.Px()<< "q4vPyL:   "<<q4v.Py()<< "q4v1PyL:   " <<q4v1.Py()<< "<q4vPzL:   "<<q4v.Pz()<< "q4v1PzL:  "<<q4v1.Pz()<< "mass diff virtual photon:  " <<(q4v.M()-q4v1.M())<<endl;
	
	//double phi = rand3->Uniform(0,2*pi);
	
	//declaration of 4-vector of 
	q4v_P.SetPxPyPzE(0,0,q4v_P.Vect().Mag(),q4v_P.E());
	q4v_H.SetPxPyPzE(0,0,q4v_H.Vect().Mag(),q4v_H.E());
	
	//cout<<"Eq:  "<<Eq<<"q4v.E():  "<<q4v.E()<< "mom tot diff:  " <<(q4v.P()-q4v1.P())<<endl;
	
	double ptot_P = q4v_P.Vect().Mag();
	double ptot_H = q4v_H.Vect().Mag();
	
	//total momentum magnitude of virtual photon proton as a target
	double qptot_P = sqrt(eseP*eseP*sin(esaP)*sin(esaP)+(Ee-eseP*cos(esaP))*(Ee-eseP*cos(esaP))); // using 3-momentum conservation of electron 
	
	
	//total momentum magnitude of virtual photon helium-4 as a target
	double qptot_H = sqrt(eseH*eseH*sin(esaH)*sin(esaH)+(Ee-eseH*cos(esaH))*(Ee-eseH*cos(esaH))); //using 3-momentum conservation of eletron
	
	//cout<<"ptot:   "<<ptot<<"qptot:  "<<qptot<<endl;
	
	//total energy of target-proton and virtual photon
	double Etot_P = Mp + q4v_P.E();
	
	//total energy of target-helium and virtual photon
	double Etot_H = Mh + q4v_H.E();
	
	
	//energy of virtual photon helium-4 as target
	double qE_P = (Mp+(Ee-eseP));
	
	//energy of virtual photon helium-4 as target
	double qE_H = (Mh+(Ee-eseH));
	//cout<<"Etot:  "<<Etot<<"qE:  "<<qE<<endl;
	
	
	double beta_P = ptot_P/Etot_P; // beta for target proton
	
	double beta_H = qptot_H/(Mh+(Ee-eseH)); // beta for target helium-4
	
	
	
	//cout<<"beta:    "<<beta<<"betaa:   "<<betaa<<endl;
	
	
	//relativistic factor gamma for target proton
	double rgamma_P = 1/sqrt(1-beta_P*beta_P);
	
	//relativistic factor gamma for target helium-4
	double rgamma_H = 1/sqrt(1-beta_H*beta_H);
	
	
	//total energy of relativistically moving particle proton as target
	double EPcm0 = rgamma_P*Mp;
	
	//total energy of relativistically moving particle helium-4 as target
	double EHcm1 = rgamma_H*Mh;
	
	
	//momentum of virtual photon in z-direction for proton as target
	double qzcmP0 = -rgamma_P*beta_P*Eq_P +(1+(rgamma_P-1))*q4v_P.Vect().Mag();
	double qzcmP1 = -rgamma_P*beta_P*Eq_P +(1+(rgamma_P-1))*qptot_P;
	
	
	
	//momentum of virtual photon in z-direction for helium-4 as target
	double qzcmH0 = -rgamma_H*beta_H*Eq_H +(1+(rgamma_H-1))*q4v_H.Vect().Mag();
	double qzcmH1 = -rgamma_H*beta_H*Eq_H +(1+(rgamma_H-1))*qptot_H;
		
	//cout<<"qzcm1:   "<<qzcm1<<"qzcm2:   "<<qzcm2<<endl;
		
		
		
		
	//transforming four-vectors from lab to CM frame using ROOT-BoostVecot() operation
	
	com4v_P = q4v_P +target4vP;
	com4v_H = q4v_H + target4vH;
	
	TVector3 boostVectorP = com4v_P.BoostVector();
	TVector3 boostVectorH = com4v_H.BoostVector(); //calculation of beta
	
	
	double boostBeta_P = boostVectorP.Mag();
	double boostBeta_H = boostVectorH.Mag();
	
	//cout<<"beta_ROOT:   "<<beta<<"betaa_manual:  "<<betaa<<"boostBeta_ROOT:   "<<boostBeta<<endl;
	com4v_P.Boost(-boostVectorP);
	com4v_H.Boost(-boostVectorH);
	
	
	q4v_P.Boost(-boostVectorP);
	q4v_H.Boost(-boostVectorH);
	
	
	target4vP.Boost(-boostVectorP);
	target4vH.Boost(-boostVectorH);
	
    // energy of helium in the center-mass frame  before collision     
	
	double PHcm = std::abs(-rgamma_H*beta_H*Mh);
	
	//cout<<"CM q4v.P():   "<<q4v.P()<<"CM target4v.P():   "<<target4vH.P()<<"CM PHcm:  "<<PHcm<<endl;

    double EHcm = target4vH.E(); 
	double Pqcm = q4v_H.Vect().Mag(); 
	//cout<<"qzcm1:   "<<qzcm1<< "qzcm2:  "<<qzcm2<<"Pqcm:   "<<Pqcm<<endl;
	
	double qEcm = rgamma_H*Eq_H -rgamma_H*beta_H*qptot_H;
	double Ecmtot1 = qEcm+EHcm;
	 
	double EcmtotP = q4v_P.E()+target4vP.E();                                                                                         
	double EcmtotH = q4v_H.E()+target4vH.E();
	
	double EcmP = EcmtotP;
	double EcmH = EcmtotH;
	//cout<<"Ecmtot:   "<<Ecmtot<<"Ecm:   "<<Ecm<<endl;
	
    // scattered helium-4 energy after collision
	double speP = (EcmP*EcmP+Mp*Mp-m*m)/(2*EcmP);                                                                                                    
	double speH = (EcmH*EcmH+Mh*Mh-m*m)/(2*EcmH);
	
	
	
	
	/*
	
	//calculation of effective t
    double tmin_eff = ((Mh*Mh) -(target4vH.E()*speH) - sqrt((target4vH.E()*target4vH.E()-Mh*Mh)*(speH*speH-Mh*Mh)))/8.0;
    double tmax_eff = ((Mh*Mh) -(target4vH.E()*speH) + sqrt((target4vH.E()*target4vH.E()-Mh*Mh)*(speH*speH-Mh*Mh)))/8.0;
	
	
	double t_eff = rand3->Uniform(tmin_eff,tmax_eff);
	
	cout<<"tmin_eff:  "<<tmin_eff<<"t_eff:   "<<t_eff<<"tmax_eff:  "<<tmax_eff<<endl;
	
	double theta_hh = acos((t_eff-2*Mh*Mh +2*target4vH.E()*speH)/(2*sqrt((target4vH.E()*target4vH.E()-Mh*Mh)*(speH*speH-Mh*Mh))));
    double theta_qhr = 180*D2R-theta_hh;

	*/





	//calculation of effective t and t values 

    double tmin_eff = ((Mh*Mh) -(target4vH.E()*speH) - sqrt((target4vH.E()*target4vH.E()-Mh*Mh)*(speH*speH-Mh*Mh)))/8.0;
    double tmax_eff = ((Mh*Mh) -(target4vH.E()*speH) + sqrt((target4vH.E()*target4vH.E()-Mh*Mh)*(speH*speH-Mh*Mh)))/8.0;

	
	double t_eff = rand3->Uniform(tmin_eff,tmax_eff);
	htEff_generated->Fill(t_eff);
	
	
	//cout<<"tmin_eff:  "<<tmin_eff<<"t_eff:   "<<t_eff<<"tmax_eff:  "<<tmax_eff<<endl;
	

    //calculation of t for helium-4
	double tmax_h = 2*Mh*Mh -2*target4vH.E()*speH - 2*sqrt((target4vH.E()*target4vH.E()-Mh*Mh)*(speH*speH-Mh*Mh));
	double tmin_h = 2*Mh*Mh -2*target4vH.E()*speH + 2*sqrt((target4vH.E()*target4vH.E()-Mh*Mh)*(speH*speH-Mh*Mh)); // t is minimum in terms of absolute of t
	
	
	 
	 double t_eff_check = rand3->Uniform(tmin_h/16,tmax_h/16);
	 
	 
	//double t = rand3->Uniform(tmin_H,tmin_H);
	
	double t_h =rand3->Uniform(tmin_h,tmax_h);
	//htH_generated->Fill(t_h);
	
	//cout<<"tmin_H:  "<<tmin_H<<"t_h:   "<<t_h<<"tmin_H:  "<<tmin_H<<endl;
	
	h_tmin_vs_th->Fill(-tmin_h,-t_h);
	
	h_minus_tmin->Fill(-tmin_h);
	
	double t_h_check =t_h/16;
	//double tdiff = std::abs(t_h-tmin_H);
	 double tdiff = (tmin_h-t_h);
	 
	// if (!std::isfinite(tdiff)) continue;
	 
	 if (!std::isfinite(tdiff) || !std::isfinite(t_h)) continue;
	 
	 // Optional before-cut histograms
	 hMinusTBefore->Fill(-t_h);
	 hTprimeBefore->Fill(tdiff);
	 
	 if (tdiff < 0.0 || tdiff > 1.0) continue;
	  ++nPassedTprime;
	  
	  h_tmin_vs_th_after->Fill(-tmin_h,-t_h);
	 
	 //h_tdiff_vs_theta_after->(tdiff, scat4vH.Theta()*R2D);
	 
	 // Fill the selected distributions only after the event passes the cut.
	 hMinusTAfter->Fill(-t_h);
	 hTprimeAfter->Fill(tdiff);
	 
	 
	//cout<<"tdiff:   "<<tdiff<<endl;
	 
	 //if (tdiff < 0.0 || tdiff > 1.0) continue;


    //if (std::isnan(tdiff)) continue;
	 
	//cout<<"tmin_H-t_h:    "<<tdiff<<"Q2H:   "<<Q2H<<endl;
	 
	//cout<<"tmin_H_after_cut:  "<<tmin_H<<"t_h_after_cut:   "<<t_h<<"tmin_H_after_cut:  "<<tmin_H<<endl;
	
	//double theta_hh = acos((t_eff-2*Mh*Mh +2*target4vH.E()*speH)/(2*sqrt((target4vH.E()*target4vH.E()-Mh*Mh)*(speH*speH-Mh*Mh))));
	
	double theta_pp = acos((t_eff-2*Mp*Mp +2*target4vP.E()*speP)/(2*sqrt((target4vP.E()*target4vP.E()-Mp*Mp)*(speP*speP-Mp*Mp))));
	double theta_hh = acos((t_h-2*Mh*Mh +2*target4vH.E()*speH)/(2*sqrt((target4vH.E()*target4vH.E()-Mh*Mh)*(speH*speH-Mh*Mh))));
	
	//cout<<"theta_hh:    "<<theta_hh<<endl;
	
	
	double theta_qpr = 180*D2R-theta_pp;
    double theta_qhr = 180*D2R-theta_hh;
	
    double p0 = sqrt(speP*speP-Mp*Mp);
    double p4 = sqrt(speH*speH-Mh*Mh);
	
    double spphi = rand3->Uniform(0,2*pi);
	
	
	
	
	scat4vP.SetXYZM(0,p0*sin(theta_qpr),p0*cos(theta_qpr),Mp);
    //scat4vH.SetXYZM(0,p4*sin(theta_qhr),p4*cos(theta_qhr),Mh);
	//scat4vH.SetXYZM(0,p4*sin(theta_qhr),p4*cos(theta_qhr),speH);
	scat4vH.SetXYZM(p4 * sin(theta_qhr) * cos(spphi),p4 * sin(theta_qhr) * sin(spphi),p4 * cos(theta_qhr),Mh);
	
	 //scat4vH.RotateZ(spphi);
	 
	double PH = scat4vH.P();

	 
	 //phi4v = beam4v + target4vH - scat4vH - scatelec4v;
	 
	 phi4vH = q4v_H + target4vH - scat4vH;
	  //cout<<"scat4vH.P():   "<<scat4vH.P()<<"phi4v.P():   "<<phi4v.P()<<endl;
	 
	 double phiE = phi4vH.E();
	 double PPhi = phi4vH.P();
	 
	 double EHP = scat4vH.E()+phi4vH.E();
	 double PHP = scat4vH.P()+phi4vH.P();
	 
	 double Ecmtot3 = scat4vH.E()+phi4vH.E();
	 
	 //cout<<"Ecm:  "<<Ecm<<"Ecm1:   "<<Ecm1<<endl;
	 
	  //TLorentzVector check4v;
	 
	 check4v1 = q4v_H+target4vH-scat4vH-phi4vH;
	 //cout<<"check4v1PxCM:   "<<check4v1.Px()<<"check4v1PyCM:     "<<check4v1.Py()<< "check4v1PzCM:   "<<check4v1.Pz()<< "check4v1MCM:  "<<check4v1.M()<<endl;
	 
	 
	 
	 //back to the lab frame from CM frame
	 
	 
	 //q4v.Boost(boostVector1);
	 
	 //cout<<"q4v1PxL:   "<<q4v1.Px()<< "q4vPxL:  "<<q4v.Px()<< "q4v1PyL:  "<<q4v1.Py()<<"q4vPyL:  "<<q4v.Py()<<endl;
	 
 	 q4v_P.Boost(boostVectorP);
 	 q4v_H.Boost(boostVectorH);
	 
	 
	 target4vP.Boost(boostVectorP);
     scat4vP.Boost(boostVectorP);
	 phi4vP.Boost(boostVectorP);
	 
	 target4vH.Boost(boostVectorH);
     scat4vH.Boost(boostVectorH);
	 phi4vH.Boost(boostVectorH);
	 
	 //cout<<"scat4vH.P():   "<<scat4vH.P()<<endl;
           
    //if(scat4vH.P()>=0.23 && scat4vH.P()<=0.4){                                                  

		//TLorentzVector scat4vH1 = scat4vH;
		
		//PH->Fill(scat4vH.P());
	
       //cout<<"target:  "<<target4v.Px() << "   "<<target4v.Py()<<"  "<<target4v.Pz()<<"  "<< target4v.E()<<" "<<target4v.M()<<endl;                                                     



     TVector3 q3v_nw0;
     q3v_nw0.SetXYZ(0,-eseP*sin(esaP), Ee-eseP*cos(esaP));

     TVector3 direction0 = q3v_nw0.Unit();
     scat4vP.RotateUz(direction0);

     phi4vP = beam4v + target4vP - scat4vP - scatelec4vP;
     //cout<<"phiLX: "<<phi4v0.Px()<< "phi4vLY:   "<<phi4v0.Py()<<"phi4vLZ:  "<<phi4v0.Pz()<<"PHiLE:  "<< phi4v0.E()<<"PhiLM:   "<<phi4v0.M()<<endl; 
     //TLorentzVector scat4vH1 = scat4vH;

     TVector3 q3v_nw;
     q3v_nw.SetXYZ(0,-eseH*sin(esaH), Ee-eseH*cos(esaH));

     TVector3 direction = q3v_nw.Unit();

      //direction.Print();                                                                                                                                                               
     scat4vH.RotateUz(direction);
	 
	 //q4v.RotateUz(direction);
	 
	 //cout<<"q4v1PxL:   "<<q4v1.Px()<< "q4vPxL:  "<<q4v.Px()<< "q4v1PyL:  "<<q4v1.Py()<<"q4vPyL:  "<<q4v.Py()<< "q4v1PzL:  " <<q4v1.Pz()<<"q4vPzL:  " <<q4v.Pz()<<endl;
	 
	 //cout<<"scat4vH.P():   "<<scat4vH.P()<<endl;

     //cout<<"scatprotX2:  "<<scatprot4v.Px()<<"scatprotY:  "<<scatprot4v.Py()<<" "<<scatprot4v.Pz()<<" "<<scatprot4v.M()<<endl;                                                        

	 //Momentum cut for recoiled helium-4 nucleus
	 //if(scat4vH.P() <= 0.23 || scat4vH.P() > 0.4 ) continue; 
	 
	  //if(scat4vH.Theta() <= TMath::Pi()/4 || scat4vH.Theta() > TMath::Pi()/2 ) continue; 

      //cout<<"tdiff:   "<<tdiff<<endl;


	 phi4vH = beam4v + target4vH - scat4vH - scatelec4vH;
	 
	 
	
	 
	 //cout<<"phi4vML:  "<<phi4v.M()<<endl;
	 //phi4v = q4v + target4vH - scat4vH;
	 //cout<<"phiLX: "<<phi4v.Px()<< "phi4vLY:   "<<phi4v.Py()<<"phi4vLZ:  "<<phi4v.Pz()<<"PHiLE:  "<< phi4v.E()<<"PhiLM:   "<<phi4v.M()<<endl;  
	  
	 //cout<<"scat4vH.P():   "<<scat4vH.P()<<"phi4v.P():   "<<phi4v.P()<<endl;
	 
     //TLorentzVector q4v2 = beam4v -scatelec4v;
   
    //cout<<"q4v2Px:   "<<q4v2.Px()<< "q4v2Py:   "<<q4v2.Py()<< "q4v2Pz:   "<<q4v2.Pz()<< "q4v2M:  "<<q4v2.M()<<endl;
   
    //cout<<"q4v2Px:   "<<q4v2.Px()<< "q4v2Py:   "<<q4v2.Py()<< "q4v2Pz:   "<<q4v2.Pz()<< "q4v2M:  "<<q4v2.M()<<endl;
	                                         
	  double sephi = rand3->Uniform(0,2*pi);
	  
	  scatelec4vP.RotateZ(sephi);
	  q4v_P.RotateZ(sephi);
      scat4vP.RotateZ(sephi);
	  phi4vP.RotateZ(sephi);
	  
	  q4v_H.RotateZ(sephi);
	  scatelec4vH.RotateZ(sephi);
      scat4vH.RotateZ(sephi);
	  phi4vH.RotateZ(sephi);
	  
	  
	  q4v_H  = beam4v - scatelec4vH;
	  phi4vH = q4v_H + target4vH - scat4vH;
	  
	  double pHe = scat4vH.P();  // recoil He-4 momentum in GeV/c
	  
	  // Check generated helium t
	  double t_check = (target4vH - scat4vH).M2();

	  double delta_t = t_check - t_h;
	  
	  if (!std::isfinite(t_check) ||
	      !std::isfinite(t_h) ||
	      !std::isfinite(delta_t)) {

	      std::cerr << "WARNING: non-finite t value"
	                << "  t_h=" << t_h
	                << "  t_check=" << t_check
	                << "  delta_t=" << delta_t
	                << std::endl;

	      tFailureCount++;
	  }
	  
	  
	  
	  
	  if (tCheckCount < 20) {

	      std::cout << std::scientific
	                << std::setprecision(10)

	                << "T-CHECK event="
	                << tCheckCount + 1

	                << "  t_h=" << t_h
	                << "  t_check=" << t_check
	                << "  delta_t=" << delta_t

	                << "  -t_h=" << -t_h
	                << "  pHe=" << scat4vH.P()
	                << "  thetaHe="
	                << scat4vH.Theta()*180.0/TMath::Pi()

	                << std::defaultfloat
	                << std::endl;
	  }
	  
	  
	  if (std::isfinite(delta_t) &&
	      std::abs(delta_t) > tTolerance) {

	      tFailureCount++;

	      if (tFailureCount <= 20) {

	          std::cerr << std::scientific
	                    << std::setprecision(10)

	                    << "WARNING: t inconsistency"

	                    << "  t_h=" << t_h
	                    << "  t_check=" << t_check
	                    << "  delta_t=" << delta_t

	                    << "  recoilE=" << scat4vH.E()
	                    << "  recoilP=" << scat4vH.P()
	                    << "  recoilM=" << scat4vH.M()

	                    << std::defaultfloat
	                    << std::endl;
	      }
	  }

	  tCheckCount++;
	  
	  
	  
	  
	  
	  
	 

	  //if (pHe < 0.2 || pHe > 2.0) continue;
	  
	  
	  
	  //phi4v1 = beam4v + target4vH - scat4vH - scatelec4v;
	   //cout<<"q4v1PxL:   "<<q4v1.Px()<< "q4vPxL:  "<<q4v.Px()<< "q4v1PyL:  "<<q4v1.Py()<<"q4vPyL:  "<<q4v.Py()<< "q4v1PzL:  " <<q4v1.Pz()<<"q4vPzL:  " <<q4v.Pz()<<endl;
	  
	  //cout<<"phiX_L: "<<phi4v.Px() <<" PhivY_L:   "<<phi4v.Py()<<" PhivZ_L:  "<<phi4v.Pz()<<"phiE1:  "<< phi4v.E()<<"PhivM1: "<<phi4v.M()<<endl;
	  
	  TLorentzVector q4v0 = beam4v -scatelec4vP;
	  TLorentzVector q4v1 = beam4v -scatelec4vH;
	  
	  //cout<<"q4v1PxL:   "<<q4v1.Px()<< "q4v2PxL:  "<<q4v2.Px()<< "q4v1PyL:  "<<q4v1.Py()<<"q4v2PyL:  "<<q4v2.Py()<< "q4v1PzL:  " <<q4v1.Pz()<<"q4v2PzL:  " <<q4v2.Pz()<<endl;
	   
	  //cout<<"q4v2Px:   "<<q4v2.Px()<< "q4v2Py:   "<<q4v2.Py()<< "q4v2Pz:   "<<q4v2.Pz()<< "q4v2M:  "<<q4v2.M()<<endl;
	   
	   //cout<<"q4v2Px:   "<<q4v2.Px()<< "q4v2Py:   "<<q4v2.Py()<< "q4v2Pz:   "<<q4v2.Pz()<< "q4v2M:  "<<q4v2.M()<<endl;
	  
	  check4v2 = beam4v+target4vH-scatelec4vH-scat4vH-phi4vH;
	 //check4v2 = q4v2+target4vH-scat4vH-phi4v;
	  //cout<<"check4v2PxL1:   "<<check4v2.Px()<<"check4v2PyL1:     "<<check4v2.Py()<< "check4v2PzL1:   "<<check4v2.Pz()<< "check4v2ML1:  "<<check4v2.M()<<endl;
	  
	  
	  //if(scat4vH.P() <= 0.24 || scat4vH.P()>0.4) continue;
	  
	  
	  //cout<<"scat4vH.P():   "<<scat4vH.P()<<endl;

	  //phi4v1= beam4v + target4vH - scat4vH - scatelec4v;
	  
	  
	   //cout<<"phi4vML:  "<<phi4v.M()<<endl;
	  
	  
	  
	  //agnels
	  
	  //double sha = acos(((beam4v.Vect()).(scat4vH.Vect()))/((beam4v.Vect().Mag()).(scat4vH.Vect())));
	  //cout<<"sha:    "<<sha/D2R<<endl;
	  
	  double sha = beam4v.Vect().Angle(scat4vH.Vect());
	  //cout << "sha: " << sha / D2R << endl;  // Converts from radians to degrees
	  
	  //kinetic energy of scattered helium-4
	  
	  double ekh = scat4vH.E()-Mh;
	  
	  //transverse momentum of scattered helium-4
	  
	  double pt = scat4vH.Pt(); // ROOT method
	  //cout<<"pt:   "<<pt<<endl;
	  
	  
	  
	 //boosting back phi meson to CM frame
	  
	  
      TVector3 boostVector2 = phi4vP.BoostVector();
      phi4vP.Boost(-boostVector2); //meson in the phi rest frame 
	  
     // TVector3 boostVector3 = phi4vH.BoostVector();
      //phi4vH.Boost(-boostVector3); //meson in the CM frame 
	  
	  //cout<<"phi4vPxCM:    " <<phi4v.Px()<<"phi4vPyCM:  " <<phi4v.Py()<<"phi4vPzCM:   " <<phi4v.Pz()<<"phi4vMCM:  " <<phi4v.M()<<endl;
	  
	  //double phiE = phi4v.E();
	  

      TLorentzVector  kaon4vP1, kaon4vP2;

      //phi4v = kaon4v1+kaon4v2;                                                                                                                                                 


      //cout<<"mag. of momentum of phi-meson(CM frame):  "<<phi4v.Vect().Mag()<<endl;                                                                                            
      //double thetaK = rand3->Uniform(0,pi);                                                                                                                                    

      //double sinetheta = rand3->Uniform(-1,1);
	  //double theta = asin(sinetheta);
      //double costheta = cos(theta);
	  
	  //isotropic phi generation
	  
	  double costhetaP = rand3->Uniform(-1.0, 1.0);
	  double sinthetaP = sqrt(1.0 - costhetaP*costhetaP);

    //double costheta = rand3->Uniform(-1,1);                                                                                                                                  


     double phiK_P = rand3->Uniform(0,2*pi);

     double pk_P = 0.1269;//in the unit of GeV/c^2                                                                                                                                

     double Ek_P = 0.51;//energy of kaons in the unit of GeV/c                                                                                                                    

      //kaon4v1 = kaon4v1.SetPxPyPzE(pk*sin(thetaK)*cos(phiK), pk*sin(thetaK)*sin(phiK),pk*cos(thetaK),0.51);                                                                    


     kaon4vP1.SetPxPyPzE(pk_P*sinthetaP*cos(phiK_P), pk_P*sinthetaP*sin(phiK_P),pk_P*costhetaP,Ek_P);

     kaon4vP2.SetPxPyPzE(-kaon4vP1.Px(), -kaon4vP1.Py(),-kaon4vP1.Pz(),kaon4vP1.E());
	 
	 phi4vP.Boost(boostVector2);
	 //phi4vH.Boost(boostVector3);
	 
	 //cout<<"phi4vML:  "<<phi4v.M()<<endl;
	 
	 kaon4vP1.Boost(boostVector2);
	 kaon4vP2.Boost(boostVector2);
	 
	 //kaon4v1.Boost(boostVector2);
	 //kaon4v2.Boost(boostVector2);
	 
	 
	 
	 
	 
	 TLorentzVector kaon4vH1, kaon4vH2;

	 // check parent phi before decay
	 if (!std::isfinite(phi4vH.Px()) || !std::isfinite(phi4vH.Py()) ||
	     !std::isfinite(phi4vH.Pz()) || !std::isfinite(phi4vH.E())) continue;

	 if (phi4vH.E() <= 0) continue;
	 if (phi4vH.M2() <= 0) continue;
	 if (phi4vH.E() < phi4vH.P()) continue;

	 // boost vector of helium phi
	 TVector3 boostVector3 = phi4vH.BoostVector();
	 if (!std::isfinite(boostVector3.X()) ||!std::isfinite(boostVector3.Y()) ||!std::isfinite(boostVector3.Z())) continue;

	 // go to phi rest frame
	 phi4vH.Boost(-boostVector3);

	 // kaon masses and two-body decay momentum
	 double Mphi_evt = phi4vH.M();
	 double MK = 0.493677;
	 if (Mphi_evt < 2.0*MK) continue;
	 
	 //double MK = 0.493677;
	 double Ek = Mphi_evt / 2.0;
	 double pk = sqrt(Ek*Ek - MK*MK);
	 //cout<<"pk:  "<<pk<<endl;
	 
	 
	 // isotropic decay
	 double costheta = rand3->Uniform(-1.0, 1.0);
	 double sintheta = sqrt(1.0 - costheta*costheta);
	 double phiK = rand3->Uniform(0, 2*pi);

	 kaon4vH1.SetPxPyPzE(pk*sintheta*cos(phiK),pk*sintheta*sin(phiK),pk*costheta,Ek);

	 kaon4vH2.SetPxPyPzE(-kaon4vH1.Px(),-kaon4vH1.Py(),-kaon4vH1.Pz(),kaon4vH1.E());

	 // boost back to lab
	 phi4vH.Boost(boostVector3);
	 kaon4vH1.Boost(boostVector3);
	 kaon4vH2.Boost(boostVector3);
	 
	 
	 
	 
	 // ============================================================
	 // Four-momentum conservation check in the laboratory frame
	 //
	 // Initial:
	 //     beam electron + target He-4
	 //
	 // Final:
	 //     scattered electron + recoil He-4 + K+ + K-
	 // ============================================================

	 TLorentzVector initial4v = beam4v + target4vH;

	 TLorentzVector final4v = scatelec4vH + scat4vH + kaon4vH1 + kaon4vH2;

	 TLorentzVector delta4v = initial4v - final4v;

	 // Component residuals
	 double deltaPx = delta4v.Px();
	 double deltaPy = delta4v.Py();
	 double deltaPz = delta4v.Pz();
	 double deltaE  = delta4v.E();
	 
	 if(std::isfinite(deltaE)) {
		 h_deltaE->Fill(deltaE);
	 }

	 double deltaP = delta4v.Vect().Mag();
	 
	 if(std::isfinite(deltaP)) {
	     h_deltaP->Fill(deltaP);
	 }

	 // Verify that the kaon pair reproduces the parent phi
	 TLorentzVector kaonPair4v = kaon4vH1 + kaon4vH2;

	 TLorentzVector phiDecayResidual4v = phi4vH - kaonPair4v;

	 double phiDecayDeltaPx = phiDecayResidual4v.Px();

	 double phiDecayDeltaPy = phiDecayResidual4v.Py();

	 double phiDecayDeltaPz = phiDecayResidual4v.Pz();

	 double phiDecayDeltaE = phiDecayResidual4v.E();

	 // On-shell mass checks
	 double recoilHeMass = scat4vH.M();

	 double kaon1Mass = kaon4vH1.M();

	 double kaon2Mass = kaon4vH2.M();

	 double reconstructedPhiMass = kaonPair4v.M();
	 
	 if(std::isfinite(reconstructedPhiMass)) {
	     h_MKK->Fill(reconstructedPhiMass);
	 }
	 
	 
	 
	 
	 
	 if (conservationDebugCount < 20) {

	     std::cout << std::scientific
	               << std::setprecision(6)

	               << "CONSERVATION event="
	               << conservationDebugCount + 1

	               << "  dPx=" << deltaPx
	               << "  dPy=" << deltaPy
	               << "  dPz=" << deltaPz
	               << "  dE="  << deltaE
	               << "  |dP|=" << deltaP

	               << "  phiDecay_dPx=" << phiDecayDeltaPx
	               << "  phiDecay_dPy=" << phiDecayDeltaPy
	               << "  phiDecay_dPz=" << phiDecayDeltaPz
	               << "  phiDecay_dE="  << phiDecayDeltaE

	               << "  MHe=" << recoilHeMass
	               << "  MK1=" << kaon1Mass
	               << "  MK2=" << kaon2Mass
	               << "  MKK=" << reconstructedPhiMass

	               << std::endl;
	 }

	 conservationDebugCount++;
	 
	 
	 
	 
	 //const double conservationTolerance = 1.0e-8;  // GeV or GeV/c
	 
	 
	 bool conservationOK =
	         std::abs(deltaPx) < conservationTolerance &&
	         std::abs(deltaPy) < conservationTolerance &&
	         std::abs(deltaPz) < conservationTolerance &&
	         std::abs(deltaE)  < conservationTolerance;

	 bool phiDecayOK =
	         std::abs(phiDecayDeltaPx) < conservationTolerance &&
	         std::abs(phiDecayDeltaPy) < conservationTolerance &&
	         std::abs(phiDecayDeltaPz) < conservationTolerance &&
	         std::abs(phiDecayDeltaE)  < conservationTolerance;

	 if (!conservationOK || !phiDecayOK) {

	     conservationFailureCount++;

	     if (conservationFailureCount <= 20) {
	         std::cerr << "WARNING: four-momentum failure"
	                   << " dPx=" << deltaPx
	                   << " dPy=" << deltaPy
	                   << " dPz=" << deltaPz
	                   << " dE="  << deltaE
	                   << std::endl;
	     }

	     // For debugging, do not reject the event yet.
	     // continue;
	 }
	 
	 
	 
	 //std::cout << "Conservation checks performed: "<< conservationDebugCount << std::endl;

	 //std::cout << "Conservation failures: " << conservationFailureCount << std::endl;
	 
	 
	 
	 // t from the recoil helium
	 double t_recoil = (scat4vH - target4vH).M2();

	 // t from the virtual photon and phi meson
	 double t_phi =(q4v_H - phi4vH).M2();

	 // Difference between the two definitions
	 double delta_t_recoil_phi = t_recoil - t_phi;

	 if(std::isfinite(delta_t_recoil_phi)) {
	     h_delta_t_recoil_phi->Fill(delta_t_recoil_phi);
		 
		 if(NEVENTS < 10) {
		        std::cout
		            << "t_recoil = " << t_recoil
		            << ", t_phi = " << t_phi
		            << ", delta_t = "
		            << delta_t_recoil_phi
		            << std::endl;
		    }
	 }
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 // Final lab-frame kaon kinematics exactly as written to LUND
	 double pKplus     = kaon4vH1.P();
	 double thetaKplus = kaon4vH1.Theta() * R2D;

	 double pKminus     = kaon4vH2.P();
	 double thetaKminus = kaon4vH2.Theta() * R2D;
	 
	 pKplus_generated->Fill(pKplus);
	 thetaKplus_generated->Fill(thetaKplus);

	 pKminus_generated->Fill(pKminus);
	 thetaKminus_generated->Fill(thetaKminus);
	 
	 pVsThetaElec->Fill(scatelec4vH.Theta()*R2D,scatelec4vH.P());
	 pVsThetaKplus->Fill(thetaKplus, pKplus);
	 pVsThetaKminus->Fill(thetaKminus, pKminus);
	 
	 
	 h_tmin_vs_recoil_theta->Fill(-tmin_h,scat4vH.Theta()*R2D);
	 
	 // Before the cuts
	 h_t_vs_recoil_theta->Fill(-t_h, scat4vH.Theta()*R2D);
	 
	 h_tdiff_vs_theta_after->Fill(tdiff, scat4vH.Theta()*R2D);

	 // Fill histograms here

	 // === Write LUND header ===
	 //myfile << "4 0 0 0 0 0 0" << std::endl;
	 
	 
	 
	 
	 
	
	 //lund file
	 
	 //double vx =0.0;
	 //double vy = 0.0;
	 //double vz = 0.0;
	 //double T = 0.0;
	/* 
	 const double vx = 0.0;
	 const double vy = 0.0;
	 const double vz = zVertexDist(rng);
	 
     // === Write LUND header ===
    // myfile << "3 0 0 0 0 0 0" << endl;
	 // event header: 4 final-state particles
	 myfile << "4 0 0 0 0 0 0" << std::endl;
	 
	 
	 //myfile << " 1 " << 321 << " 0 0 "<< kaon4vH1.Px() << " " << kaon4vH1.Py() << " " << kaon4vH1.Pz() << " "<< kaon4vH1.E()  << " " << kaon4vH1.M()  << " "<< vx << " " << vy << " " << vz << " " << T << endl;

	 //myfile << " 1 " << -321 << " 0 0 "<< kaon4vH2.Px() << " " << kaon4vH2.Py() << " " << kaon4vH2.Pz() << " "<< kaon4vH2.E()  << " " << kaon4vH2.M()  << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	 
	 
	 
      // === Write initial state ===
    //myfile << "-1 " << 11.0 << " 0 0 "<< beam4v.Px() << " " << beam4v.Py() << " " << beam4v.Pz() << " " << beam4v.E() << " "<< beam4v.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	 //myfile << "-1 " << 1000020040 << " 0 0 "<< target4vH.Px() << " " << target4vH.Py() << " " << target4vH.Pz() << " " << target4vH.E() << " "<< target4vH.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	

     // === Write final state ===
     //myfile << " 1 " << 11.0 << " 0 0 "<< scatelec4vH.Px() << " " << scatelec4vH.Py() << " " << scatelec4vH.Pz() << " " << scatelec4vH.E() << " "<< scatelec4vH.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 // scattered electron
	 myfile << "1 -1 1 11 0 0 "<< scatelec4vH.Px() << " "<< scatelec4vH.Py() << " "<< scatelec4vH.Pz() << " "<< scatelec4vH.E()  << " "<< Me  << " "<< vx << " " << vy << " " << vz<< std::endl;
	 

     //myfile << " 1 " << 321 << " 0 0 "<< kaon4vH1.Px() << " " << kaon4vH1.Py() << " " << kaon4vH1.Pz() << " " << kaon4vH1.E() << " "<< kaon4vH1.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	 // K+
	 myfile << "2 1 1 321 0 0 "<< kaon4vH1.Px() << " "<< kaon4vH1.Py() << " "<< kaon4vH1.Pz() << " "<< kaon4vH1.E()  << " "<< Mk  << " "<< vx << " " << vy << " " << vz<< std::endl;

     //myfile << " 1 " << -321 << " 0 0 "<< kaon4vH2.Px() << " " << kaon4vH2.Py() << " " << kaon4vH2.Pz() << " " << kaon4vH2.E() << " "<< kaon4vH2.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	 // K-
	 myfile << "3 -1 1 -321 0 0 "<< kaon4vH2.Px() << " "<< kaon4vH2.Py() << " "<< kaon4vH2.Pz() << " "<< kaon4vH2.E()  << " "<< Mk  << " "<< vx << " " << vy << " " << vz<< std::endl;

     //myfile << " 1 " << 1000020040 << " 0 0 "<< scat4vH.Px() << " " << scat4vH.Py() << " " << scat4vH.Pz() << " " << scat4vH.E() << " "<< scat4vH.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	 // recoil 4He
	 myfile << "4 2 1 1000020040 0 0 "<< scat4vH.Px() << " "<< scat4vH.Py() << " "<< scat4vH.Pz() << " "<< scat4vH.E()  << " "<< Mh  << " "<< vx << " " << vy << " " << vz<< std::endl;

    


	 ++nSelected;

	 if (nSelected % 100000 == 0 || nSelected == NEVENTS) {
	     std::cout << "Accepted " << nSelected
	               << " / " << NEVENTS
	               << " events after " << nThrown
	               << " throws"
	               << std::endl;
	 }


*/
	 

	 double eetheta = acos((beam4v.Vect().Dot(scatelec4vH.Vect()))/(beam4v.Vect().Mag()*scatelec4vH.Vect().Mag()));
	 
	 double setheta = scatelec4vH.Theta();
	 
	 //cout<<"eetheta:   "<<eetheta/D2R<<"e4.Theta():    "<<setheta/D2R<<endl;
	 
	 
	 double ehtheta = acos((beam4v.Vect().Dot(scat4vH.Vect()))/(beam4v.Vect().Mag()*scat4vH.Vect().Mag()));
	 
	 double shtheta = scat4vH.Theta();
	 //cout<<"ehtheta:   "<<ehtheta/D2R<<"h4.Theta():    "<<shtheta/D2R<<endl;
	 
	 double kktheta = acos((kaon4v1.Vect().Dot(kaon4v2.Vect()))/(kaon4v1.Vect().Mag()*kaon4v2.Vect().Mag()));
	 
	 
	 double k1theta = acos((beam4v.Vect().Dot(kaon4v1.Vect()))/(beam4v.Vect().Mag()*kaon4v1.Vect().Mag()));
	 
	 double k2theta = acos((beam4v.Vect().Dot(kaon4v2.Vect()))/(beam4v.Vect().Mag()*kaon4v2.Vect().Mag()));
	 
	 double kktheta1 = acos((beam4v.Vect().Dot(kaon4v1.Vect()))/(beam4v.Vect().Mag()*kaon4v1.Vect().Mag())) - acos((beam4v.Vect().Dot(kaon4v2.Vect()))/(beam4v.Vect().Mag()*kaon4v2.Vect().Mag()));
	 
	 
	 //cout<<"k1theta:   " <<k1theta/D2R<<"k2theta:   "<<k2theta/D2R<<endl;
	 //cout<<"kktheta:   "<<kktheta/D2R<< "kktheta1:    "<<kktheta1/D2R<<endl;




		 
     //check4v3 = beam4v+target4vH-scatelec4v-scat4vH-kaon4v1-kaon4v2;
	 check4v3 = beam4v+target4vH-scatelec4vH-scat4vH-kaon4vH1-kaon4vH2;
	 //cout<<"check4v3PxL2:   "<<check4v3.Px()<<"check4v3PyL2:     "<<check4v3.Py()<< "check4v3PzL2:   "<<check4v3.Pz()<< "check4v3ML2:  "<<check4v3.M()<<endl;
	 if (check4v3.P() > 1e-3 || fabs(check4v3.E()) > 1e-3) {
	     cout << "Warning: non-conserved event: "<< " dE=" <<check4v3.E()<< " dP=" << check4v3.P() << endl;
	 }
	 
	 //form factors
	 
	 
     //double cffh = (1-pow(2.5*t_h,6))*exp(11.7*t_h);
	 double cffh = (1-pow(2.5*(t_h-tmin_h),6))*exp(11.7*(t_h-tmin_h));
	 
	 double cffh_sq = cffh*cffh;
     double cffh_min = (1-pow(2.5*tmin_eff,6))*exp(11.7*tmin_eff);
	 
	 
     //double cffp = (1/(1-Q2/(0.84*0.84)))*(1/(1-Q2/(0.84*0.84)));
	 
	 
     //for t-channel
	 
	 double cffpt = (1/(1-t_eff/(0.84*0.84)));
     //double cffpt_sq = (1/(1-t_eff/(0.84*0.84)))*(1/(1-t_eff/(0.84*0.84)));
	 double cffpt_sq = cffpt*cffpt;
 
 
	 //Longitudinal and transverse cross-section
	 
	 
     double alpha1 = 400.0;  //in the unit of nanobarns
     double alpha2 = 1.0;
     double alpha3 = 0.32;





	 /*

	 //cross-section for proton as a target


	 double wt_p = Mp + m;
	 double wP = sqrt(Mp*Mp + 2.0*Mp*(Ee-eseP) - Q2P);
	 double wr_p = wt_p / wP;

	 double ct_p  = alpha1*(1.0-wr_p*wr_p)*pow(wP,0.32);
	 double cst_p = ct_p/pow((1.0+Q2P/(m*m)),3.0);
	 
	 //double fixed_Q2P = 2.5; 
	 //double cst_p = ct_p/pow((1.0+fixed_Q2P/(m*m)),3.0);
	 double csl_p = ((0.4*Q2P)/(m*m))*cst_p;
	 double R_p = csl_p/cst_p;

	 double yP = 1.0 - eseP/Ee;
	 double eps_p = (1.0 - yP - Q2P/(4.0*Ee*Ee)) / (1.0 - yP + 0.5*yP*yP + Q2P/(4.0*Ee*Ee));
	 
     double B_initial =  2.2; // in unit of 1/GeV^2
     double alpha4 = 0.24;
     double Bp = B_initial + 4*alpha4*(TMath::Log(wP));
	 
	 
	 // Effective proton momentum transfer
	 double tp = t_h/16.0;
	 
	 
	 
	 
	 //exponential dependence of dsigma/dt 
	 
     double Ft_p = exp(Bp*tp);
     //double Fint_p = exp(Bp*tmin_eff)/Bp;
	 
	 double tpNear = tmin_h/16.0;
	 double tpFar = tmax_h/16.0;
	 
	 double Fint_p =(std::exp(Bp*tpNear) - std::exp(Bp*tpFar))/Bp;

	 //double weightP = (cst_p + eps_p*csl_p) * exp(Bp * tp);                 //this is dsigmaPdtp -diff. for proton 
	 const double alpha = 1.0 / 137.035999084;
	 
	 double Gamma = (alpha / (8.0 * M_PI)) *(Q2P / (Mp * Mp * 11.0 * 11.0)) *((1.0 - xp) / std::pow(xp, 3)) * (1.0 / (1.0 - eps_p));
	 
	 double weightP = (cst_p + eps_p*csl_p) * Ft_p/Fint_p;                 //this is dsigmaPdtp -diff. for proton 
	 
	 */
	 
	 
	 
	 
	 
	 // ------------------------------------------------------------
	 // Proton cross section
	 // ------------------------------------------------------------

	 double wt_p = Mp + m;

	 double wP = sqrt(Mp*Mp + 2.0*Mp*(Ee - eseP) - Q2P);

	 double wr_p = wt_p / wP;

	 // W dependence
	 double ct_p = alpha1 * pow(1.0 - wr_p*wr_p, alpha2) * pow(wP, alpha3);

	 // sigma_T
	 double cst_p = ct_p / pow(1.0 + Q2P/(m*m), 3.0);

	 // sigma_L
	 double csl_p = (0.4 * Q2P/(m*m)) * cst_p;
	 
	  double R_p = csl_p/cst_p;


	 // ------------------------------------------------------------
	 // Virtual-photon polarization epsilon
	 // ------------------------------------------------------------

	 double yP = 1.0 - eseP/Ee;

	 double eps_p = (1.0 - yP - Q2P/(4.0*Ee*Ee)) /(1.0 - yP + 0.5*yP*yP + Q2P/(4.0*Ee*Ee));


	 // ------------------------------------------------------------
	 // t dependence
	 // ------------------------------------------------------------

	 double B_initial = 2.2;
	 double alpha4 = 0.24;

	 double Bp = B_initial + 4.0*alpha4*TMath::Log(wP);

	 double tp = t_h / 16.0;

	 double Ft_p = exp(Bp*tp);

	 double tpNear = tmin_h / 16.0;
	 double tpFar  = tmax_h / 16.0;

	 double Fint_p = (exp(Bp*tpNear) - exp(Bp*tpFar))/Bp;


	 // ------------------------------------------------------------
	 // Virtual-photon flux
	 // ------------------------------------------------------------

	 const double alpha = 1.0 / 137.035999084;

	 double Gamma = (alpha / (8.0*M_PI)) * (Q2P / (Mp*Mp*Ee*Ee)) * ((1.0 - xp) / pow(xp,3.0)) * (1.0 / (1.0 - eps_p));


	 // ------------------------------------------------------------
	 // d^3 sigma / (dQ^2 dx_B dt)
	 // ------------------------------------------------------------

	 double weightP = Gamma * (cst_p + eps_p*csl_p) * Ft_p/Fint_p;
	 
	 	 
	 
	 if(!std::isfinite(weightP) || weightP < 0.0) {
	     continue;
	 }
	 
	 hxp_weighted->Fill(xp,weightP);
	 htEff_weighted->Fill(t_eff,weightP);
	 
	 
	 //cross-section for helium-4 as a target
     //calculation of W and Wt for proton as a target considered
     //double wt = M+m;
     //for helium
	 
	 //tdiff = tmin_H - t_h;

	 //if (!std::isfinite(tdiff)) continue;
	 //if (tdiff < 0.0 || tdiff > 1.0) continue;
	 
	 //tdiff_generated->Fill(tdiff);
	 //tdiff_weighted->Fill(tdiff, weightH);
	 
 	 //cout<<"tmin_H-t_h:    "<<tdiff<<"tdiff:   "<<tdiff<<endl;
	 //htH_generated->Fill(t_h);
	 //htH_weighted->Fill(t_h,weightH);
	 
	 
     double wt_h = Mh+m;
     //double w = sqrt(M*M+2.0*M*(Ee-ese1)-Q2);
     //calculation of W and Wt for helium-4 as target consided
     double wH = sqrt(Mh*Mh+2.0*Mh*(Ee-eseH)-Q2H);

     //cout<<"W: "<<w<<"Wh:  "<< wH<<endl;
     //double wr = wt/w;
    //for helium
     double wr_h = wt_h/wH;
     //double term1 = M*M+2.0*M*(q4v.E());
     double q_e = Ee - eseH;
     //cout<<"q4v.E:  "<<q4v.E()<< "q_e:  "<<q_e<<endl;
     
     
     //cout<<" term1:  "<<term1<<" Q2:  "<<Q2<<endl;
     //cout<<" w: "<<w<<endl;
     //double ct1 = alpha1*(1.0-wr*wr)*pow(w,0.32);
     //cout<<"ct1: "<<ct1<<endl;

	 //for proton as a target

    //double ct = alpha1*(1 - wr*wr)*pow(w,0.32);
    //double cst = ct1/pow((1+Q2/(m*m)),3.0);
    //double csl = ((0.4*Q2)/(m*m))*cst;

    //for helium as a target
     double ct_h = alpha1*(1.0-wr_h*wr_h)*pow(wH,0.32);  //in the unit of nanobars
   //transverse and longitudinal cross-section for helium-4
     double cst_h = ct_h/pow((1.0+Q2H/(m*m)),3.0);
     double csl_h =((0.4*Q2H)/(m*m))*cst_h;
	 
	 
    //transverse and longitudinal cross-section for proton as a target
	 //double cst_p = ct;

	 //for fixed Q2_H

    //double wH_f = sqrt(Mh*Mh +2.0*);

     double gm = 1.0;// square of mass of gluon
     //cout<<"cst:  "<<cst<< "csl: "<<csl<<endl;
     
     //double dif_cs = 3*gamma_nw*(cst + esp*csl)*(pow((gm-tmax_nw),3.0)/(pow((gm-t_nw),4.0)));
     //double fact1 = gamma_nw*(cst + esp*csl);
     
     //cout<<"fact1:  "<<fact1<<endl;
     
    // double den2 = pow(gm-abs(t_nw),4.0);
     //double num2 = pow(gm-abs(tmax_nw),3.0);
     
     //cout<<"den2:  "<<den2<<"num2:  "<<num2<<endl;
     
     //cout<<"abs_tmax: "<<abs(tmax)<<"abs_t: "<<abs(t)<<endl;
     
     //cout<<"tmax: "<<abs(tmax_nw)<<endl;
     
     //cout<<"x1:  "<<x1<< "esp: "<<esp<<endl;
     //double B_initial =  2.2; // in unit of 1/GeV^2 and defined in line number 861
     //double alpha4 = 0.24; // it is declared in line number 862
     double Bh = B_initial + 4*alpha4*(TMath::Log(wH));
     double Ft_h = exp(Bh*t_h);
     double Fint_h = exp(Bh*tmin_h)/Bh;
    // double R = csl/cst;
     //R for helium
     double R_h = csl_h/cst_h;
	 
	 double R_htest = 0.4*(Q2H/(m*m));
	 
	 
	 //calculation polarization parameter
     //cross-section calulation
     double fsc = 1.0/137.0;
     double y = 1.0 - eseH/Ee;
     double y1 = 1.0 - eseH/beam4v.E();
     //cout<<"y:  "<<y<< "y1:  "<<y1<<endl;
     //double esp = (1.0-y-Q2/(4.0*Ee*Ee))/(1.0-y+(y*y)/2.0+(Q2/(4.0*Ee*Ee)));

     //for helium-4 as a target
     //double y2 = 1.0 - ese1/Ee;
     //double esp_h = (1.0 - y1-Q2H/(4.0*Ee*Ee))/(1.0-y1+(y1*y1)/2.0+(Q2H/(4.0*Ee*Ee)));
	 
	 //double y = y1;       // inelasticity
	 //double Q2 = Q2H;     // GeV^2
	 //double E = Ee;       // Beam energy in GeV

	 double esp_h = (1.0 - y1 - Q2H / (4.0 * Ee * Ee)) / (1.0 - y1 + (y1 * y1) / 2.0 + Q2H / (4.0 * Ee * Ee));
	 
   
   
    //double cross_tot = cst_h +esp_h*csl_h;	
    double cr = 0.4;  //independent of W
    double R1 = cr*Q2H/(m*m);
    double cross_tot = cst_h*(1+esp_h*R1); // in the units of nb
   
    //cout <<"tot.cross:    "<<cross_tot<<endl;
	
	
	//double FC_He4_sq = 16*(pow((1.0 - pow(2.5*t_h, 6)) * exp(11.7*t_h), 2));
	//pFCsq_t_h->Fill(-t_h,FC_He4_sq);
	
	
	double FC_He4_sq = 16*(pow((1.0 - pow(2.5*tdiff, 6)) * exp(-11.7*tdiff), 2));
	//double FC_He4_sq = 16*(pow((1.0 - pow(2.5*(t_h), 6)) * exp(-11.7*(t_h)), 2));
	pFCsq_t_h->Fill(tdiff,FC_He4_sq);
 
 
     double weightH = (cst_h + esp_h * csl_h) * exp(Bh * t_h);
	 
	 
	 
	 //dsigma/dt calculation for Helium-4
	 
	 
	 // ------------------------------------------------------------
	 // Model-dependent form factors and cross-section weights
	 // ------------------------------------------------------------

	 const double A = 4.0;

	 // Radius hypotheses
	 const double lambdaCompact = 0.93;
	 const double lambdaCharge  = 1.00;
	 const double lambdaBroad   = 1.10;

	 // Effective proton momentum transfer
	 tp = t_h/(A*A);
	 
	 
	 //double protonFt = protonCrossSectionShape(tp,Wp);
	 
	 //double sigmaU = sigmaT + epsilonP*sigmaL;
	 
	 
	 //double B = B0 + 4.0*alphaPrime*std::log(Wp/1.0);
	 
	 //double tpNear = tmin_h/16.0;
	 //double tpFar = tmax_h/16.0;
	 
	 //double Fint =(std::exp(B*tpNear) - std::exp(B*tpFar))/B;
	 
	 //double dsigmaPdtp = sigmaU*Ft/Fint;
	 
	 

	 // Proton charge form factor
	 double FCp = protonChargeFormFactor(tp);
	 
	 
	 // Helium-4 form-factor hypotheses
	 double FCHeCompact = heliumFormFactor(t_h, lambdaCompact);

	 double FCHeCharge = heliumFormFactor(t_h, lambdaCharge);

	 double FCHeBroad = heliumFormFactor(t_h, lambdaBroad);

	 double FCHeSmooth = heliumSmoothFormFactor(t_h);

	 // Complete cross-section shapes
	 double crossSectionCompact = weightP*heliumCrossSectionShape(t_h,lambdaCompact);

	 double crossSectionCharge = weightP*heliumCrossSectionShape(t_h,lambdaCharge);

	 double crossSectionBroad = weightP*heliumCrossSectionShape(t_h,lambdaBroad);
	 
	 double crossSectionSmooth = weightP*heliumSmoothCrossSectionShape(t_h);
	 
	 
	 
	 
	 
	 double minusT = -t_h;

	 bool validModelWeights =
	     std::isfinite(minusT) &&
	     std::isfinite(crossSectionCompact) &&
	     std::isfinite(crossSectionCharge) &&
	     std::isfinite(crossSectionBroad) &&
	     std::isfinite(crossSectionSmooth) &&
	     crossSectionCompact >= 0.0 &&
	     crossSectionCharge >= 0.0 &&
	     crossSectionBroad >= 0.0 &&
	     crossSectionSmooth >= 0.0;

	 if(validModelWeights) {

	     hMinusTCompactModel->Fill(minusT,crossSectionCompact);

	     hMinusTChargeModel->Fill(minusT, crossSectionCharge);

	     hMinusTBroadModel->Fill(minusT,crossSectionBroad);

	     hMinusTSmoothModel->Fill(minusT,crossSectionSmooth);
		 
	 }
	 
	 





	 const double vx = 0.0;
	 const double vy = 0.0;
	 const double vz = zVertexDist(rng);
	 
	 // event-by-event physics weight
	 //const double eventWeight = weightP;
	 
	 const double WEIGHT_SCALE = 1.0e40;
	 
	  const double eventWeight = crossSectionCharge*WEIGHT_SCALE;

	 // ============================================================
	 // LUND EVENT HEADER
	 // ============================================================

	 //myfile << 2  << " " << 4  << " " << 2  << " " << 0.0 << " "<< 0.0 << " "<< 11 << " "<< Ee << " "<< 0  << " " << 1  << " " << eventWeight << std::endl;
	  
	  const int nParticles = 4;
	  const int targetA = 4;
	  const int targetZ = 2;

	 // ============================================================
	 // LUND EVENT HEADER
	 // ============================================================

	 //myfile << 2  << " " << 4  << " " << 2  << " " << 0.0 << " "<< 0.0 << " "<< 11 << " "<< Ee << " "<< 0  << " " << 1  << " " << eventWeight << std::endl;
	  
	   myfile << nParticles << " "<< targetA << " " << targetZ << " "<< 0.0 << " "<< 0.0 << " "<< 11 << " "<< Ee << " "<< 0 << " "<< 1 << " "<< eventWeight<< std::endl;
	 
	 
	 //myfile << " 1 " << 321 << " 0 0 "<< kaon4vH1.Px() << " " << kaon4vH1.Py() << " " << kaon4vH1.Pz() << " "<< kaon4vH1.E()  << " " << kaon4vH1.M()  << " "<< vx << " " << vy << " " << vz << " " << T << endl;

	 //myfile << " 1 " << -321 << " 0 0 "<< kaon4vH2.Px() << " " << kaon4vH2.Py() << " " << kaon4vH2.Pz() << " "<< kaon4vH2.E()  << " " << kaon4vH2.M()  << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	 
	 
	 
      // === Write initial state ===
    //myfile << "-1 " << 11.0 << " 0 0 "<< beam4v.Px() << " " << beam4v.Py() << " " << beam4v.Pz() << " " << beam4v.E() << " "<< beam4v.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	 //myfile << "-1 " << 1000020040 << " 0 0 "<< target4vH.Px() << " " << target4vH.Py() << " " << target4vH.Pz() << " " << target4vH.E() << " "<< target4vH.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	

     // === Write final state ===
     //myfile << " 1 " << 11.0 << " 0 0 "<< scatelec4vH.Px() << " " << scatelec4vH.Py() << " " << scatelec4vH.Pz() << " " << scatelec4vH.E() << " "<< scatelec4vH.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 // scattered electron
	 myfile << "1 -1 1 11 0 0 "<< scatelec4vH.Px() << " "<< scatelec4vH.Py() << " "<< scatelec4vH.Pz() << " "<< scatelec4vH.E()  << " "<< Me  << " "<< vx << " " << vy << " " << vz<< std::endl;
	 

     //myfile << " 1 " << 321 << " 0 0 "<< kaon4vH1.Px() << " " << kaon4vH1.Py() << " " << kaon4vH1.Pz() << " " << kaon4vH1.E() << " "<< kaon4vH1.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	 // K+
	 myfile << "2 1 1 321 0 0 "<< kaon4vH1.Px() << " "<< kaon4vH1.Py() << " "<< kaon4vH1.Pz() << " "<< kaon4vH1.E()  << " "<< Mk  << " "<< vx << " " << vy << " " << vz<< std::endl;

     //myfile << " 1 " << -321 << " 0 0 "<< kaon4vH2.Px() << " " << kaon4vH2.Py() << " " << kaon4vH2.Pz() << " " << kaon4vH2.E() << " "<< kaon4vH2.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	 // K-
	 myfile << "3 -1 1 -321 0 0 "<< kaon4vH2.Px() << " "<< kaon4vH2.Py() << " "<< kaon4vH2.Pz() << " "<< kaon4vH2.E()  << " "<< Mk  << " "<< vx << " " << vy << " " << vz<< std::endl;

     //myfile << " 1 " << 1000020040 << " 0 0 "<< scat4vH.Px() << " " << scat4vH.Py() << " " << scat4vH.Pz() << " " << scat4vH.E() << " "<< scat4vH.M() << " "<< vx << " " << vy << " " << vz << " " << T << endl;
	 
	 // recoil 4He
	 myfile << "4 2 1 1000020040 0 0 "<< scat4vH.Px() << " "<< scat4vH.Py() << " "<< scat4vH.Pz() << " "<< scat4vH.E()  << " "<< Mh  << " "<< vx << " " << vy << " " << vz<< std::endl;

    


	 ++nSelected;

	 if (nSelected % 100000 == 0 || nSelected == NEVENTS) {
	     std::cout << "Accepted " << nSelected
	               << " / " << NEVENTS
	               << " events after " << nThrown
	               << " throws"
	               << std::endl;
	 }

















	 
	 
		 
	 
	 
 
	 hxh_weighted->Fill(xh,weightH);
	 htH_generated->Fill(t_h);
	 htH_weighted->Fill(t_h,weightH);
	 
	 tdiff_generated->Fill(tdiff);
	 //tdiff_weighted->Fill(tdiff, weightH);
	 tdiff_weighted->Fill(tdiff, FC_He4_sq);
	 
	 
	 pHe = scat4vH.P();
	 double thetaHe = scat4vH.Theta() * R2D;
	 
	 // Unweighted histograms
	 tdiff_vs_thetaHe->Fill(thetaHe, tdiff);
	 tH_vs_thetaHe->Fill(thetaHe, t_h);

	 if (std::isfinite(pHe) &&
	     std::isfinite(thetaHe) &&
	     std::isfinite(FC_He4_sq) &&
	     FC_He4_sq >= 0.0) {

	     pHe_generated->Fill(pHe);
	     pHe_weighted->Fill(pHe, FC_He4_sq);

	     thetaHe_generated->Fill(thetaHe);
	     thetaHe_weighted->Fill(thetaHe, FC_He4_sq);
	 }
	 
	 
	 h_p_vs_theta->Fill(thetaHe, pHe); // x=theta, y=p
	 h_theta_vs_vz->Fill(vz, thetaHe);     // x=vz, y=theta
	 h_p_vs_vz->Fill(vz, pHe);             // x=vz, y=p
	 
	 
	 
	 
	 // weighted W distributions must go here
	 wP_weighted->Fill(Wp, weightP);
	 wH_weighted->Fill(Wh, weightH);
 




 //Q2_P_generated->Fill(Q2P); // it is filled right after the generation
 Q2_P_Weighted->Fill(Q2P,weightP);
	
 //Q2_H_generated->Fill(Q2H);// it is filled right after the generation
 Q2_H_Weighted->Fill(Q2H,weightH);

 //xB->Fill(xh,weightH);
 THE->Fill(esaH*R2D);

 ESE_P->Fill(eseP);
 ESE_H->Fill(eseH);
 
 TE->Fill(t_eff);
 TH->Fill(t_h,weightH);
 thmin->Fill(tmin_h);
 thmax->Fill(tmax_h);
 
 PH1->Fill(scat4vH.P());
 PT->Fill(pt);
 ECMTheta->Fill(esaH/D2R,EHcm);
 
 ECM->Fill(EcmP);
 EHCM->Fill(EHcm);
 EHCM1->Fill(EHcm1);
 ECMtot1->Fill(Ecmtot1);
 ECMtot2->Fill(Ecmtot3);
 
 SHCM->Fill(speH);
 SPCM->Fill(phiE);
 
 PqCM->Fill(qzcmH1);
 
 EHPhi->Fill(EHP);
 PHPhi->Fill(PHP);
 
 PHCM1->Fill(PHcm);
 PHCM2->Fill(PH);
 PPhiCM->Fill(PPhi);
 
 //PHL->Fill(scat4vH.P());
 
 
 
 Q2HW->Fill(Q2H,Wh,FC_He4_sq);
 Q2HXB->Fill(xh,Q2H,FC_He4_sq);
 //Q2HtH->Fill(-t_h,Q2H);
 Q2HtH->Fill(Q2H,tdiff,FC_He4_sq);
 //Q2HtH->Fill(Q2H,-t_h,FC_He4_sq);
 
 //Q2HW_rec->Fill(Q2P_rec, W_rec, weightP);
 //Q2HXB_rec->Fill(xB_rec, Q2P_rec, weightP);
 
 
 
 ThetaPe->Fill(scatelec4vH.P(),scatelec4vH.Theta()/D2R);
 
 
 //cross-section
// R->Fill(Q2H,R_h);
 pRP->Fill(Q2P,R_p);
 pRH->Fill(Q2H, R_h);
 Q2Hvst_h->Fill(-t_h,Q2H);
 pSigmaT_Q2_P->Fill(Q2P, cst_p);
 CTW->Fill(Wp,cst_p);
 //CTQ2->Fill(Q2H,cst_h);
 
 //if (fabs(Wp - 2.0) < 0.25) pSigmaT_Q2_W45_P->Fill(Q2P, cst_p);
 //if (fabs(Wh - 4.5) < 0.25) pSigmaT_Q2_W45_H->Fill(Q2H, cst_h);
 
 //if (fabs(Wp - 3.0) < 0.25) pSigmaT_Q2_W50_P->Fill(Q2P, cst_p);
 //if (fabs(Wh - 5.0) < 0.25) pSigmaT_Q2_W50_H->Fill(Q2H, cst_h);
 
 //if (fabs(Wp - 4.0) < 0.25) pSigmaT_Q2_W75_P->Fill(Q2P, cst_p);
 //if (fabs(Wh - 7.5) < 0.25) pSigmaT_Q2_W75_H->Fill(Q2H, cst_h);
 
 
 
 thtmin->Fill(tdiff);
 
 Q2thtmin->Fill(Q2H,tdiff,cffh_sq);
 //Q2thtmin->Fill(Q2H,tdiff,cross_tot);
 
 //form factors
 
 CFFH->Fill(tdiff,cffh_sq);
 CFFPt->Fill(-t_eff,cffpt_sq);
 ththeta->Fill(tdiff,theta_hh/D2R);
 ththeta1->Fill(tdiff,sha/D2R);
 EKH->Fill(ekh);
 
 PHT->Fill(tdiff,scat4vH.P());
 PTT->Fill(tdiff,pt);
 PTA->Fill(sha/D2R,pt);
 PZ->Fill(sha/D2R,scat4vH.Pz());
 //f1->Draw();
 
 
 
	
}



std::cout << "\nDelta-t histogram diagnostics\n";
std::cout << "Entries   = "
          << h_delta_t_recoil_phi->GetEntries()
          << std::endl;

std::cout << "Underflow = "
          << h_delta_t_recoil_phi->GetBinContent(0)
          << std::endl;

std::cout << "Overflow  = "
          << h_delta_t_recoil_phi->GetBinContent(
                 h_delta_t_recoil_phi->GetNbinsX() + 1
             )
          << std::endl;

std::cout << "Mean      = "
          << h_delta_t_recoil_phi->GetMean()
          << " GeV^2"
          << std::endl;

std::cout << "RMS       = "
          << h_delta_t_recoil_phi->GetRMS()
          << " GeV^2"
          << std::endl;








const double fRetained =
    (nThrown > 0)
        ? static_cast<double>(nSelected) /
              static_cast<double>(nThrown)
        : 0.0;

const double fTprime =
    (nThrown > 0)
        ? static_cast<double>(nPassedTprime) /
              static_cast<double>(nThrown)
        : 0.0;

std::cout << "\n====================================\n";
std::cout << "Total thrown events:       "
          << nThrown << std::endl;

std::cout << "Passed 0 < t' < 1 GeV^2:  "
          << nPassedTprime << std::endl;

std::cout << "Selected/written events:   "
          << nSelected << std::endl;

std::cout << std::setprecision(8)
          << "t' cut retention fraction: "
          << fTprime << std::endl;

std::cout << std::setprecision(8)
          << "f_retained = Nselected/Nthrown = "
          << fRetained << std::endl;

if (nSelected > 0) {
    std::cout << "Throws per written event:  "
              << static_cast<double>(nThrown) /
                     static_cast<double>(nSelected)
              << std::endl;
}

std::cout << "====================================\n";









// Print once, after all events
std::cout << "\n====================================\n";

std::cout << "Conservation checks performed: "
          << conservationDebugCount
          << std::endl;

std::cout << "Conservation failures: "
          << conservationFailureCount
          << std::endl;

std::cout << "====================================\n";



std::cout << "\n====================================\n";

std::cout << "t checks performed: "
          << tCheckCount
          << std::endl;

std::cout << "t-check failures: "
          << tFailureCount
          << std::endl;

std::cout << "t tolerance: "
          << tTolerance
          << " GeV^2"
          << std::endl;

std::cout << "====================================\n";



//plotting section





 myfile.close();



 // =========================
 // Make clones for raw plot
 // =========================
 TH1F *raw_P_gen = (TH1F*)Q2_P_generated->Clone("raw_P_gen");
 TH1F *raw_P_wgt = (TH1F*)Q2_P_Weighted->Clone("raw_P_wgt");
 TH1F *raw_H_gen = (TH1F*)Q2_H_generated->Clone("raw_H_gen");
 TH1F *raw_H_wgt = (TH1F*)Q2_H_Weighted->Clone("raw_H_wgt");

 // ================================
 // Make separate clones to normalize
 // ================================
 TH1F *norm_P_gen = (TH1F*)Q2_P_generated->Clone("norm_P_gen");
 TH1F *norm_P_wgt = (TH1F*)Q2_P_Weighted->Clone("norm_P_wgt");
 TH1F *norm_H_gen = (TH1F*)Q2_H_generated->Clone("norm_H_gen");
 TH1F *norm_H_wgt = (TH1F*)Q2_H_Weighted->Clone("norm_H_wgt");

 // normalize only the normalized copies
 if (norm_P_gen->Integral() > 0) norm_P_gen->Scale(1.0 / norm_P_gen->Integral());
 if (norm_P_wgt->Integral() > 0) norm_P_wgt->Scale(1.0 / norm_P_wgt->Integral());
 if (norm_H_gen->Integral() > 0) norm_H_gen->Scale(1.0 / norm_H_gen->Integral());
 if (norm_H_wgt->Integral() > 0) norm_H_wgt->Scale(1.0 / norm_H_wgt->Integral());


 TCanvas *c_raw = new TCanvas("c_raw", "Q2 comparison raw", 1000, 800);

 raw_P_gen->SetLineColor(kBlack);
 raw_P_wgt->SetLineColor(kRed);
 raw_H_gen->SetLineColor(kBlue);
 raw_H_wgt->SetLineColor(kGreen+2);

 raw_P_gen->SetLineWidth(2);
 raw_P_wgt->SetLineWidth(2);
 raw_H_gen->SetLineWidth(2);
 raw_H_wgt->SetLineWidth(2);

 gStyle->SetOptStat(0);

 double maxy_raw = std::max(
     std::max(raw_P_gen->GetMaximum(), raw_P_wgt->GetMaximum()),
     std::max(raw_H_gen->GetMaximum(), raw_H_wgt->GetMaximum())
 );
 raw_P_gen->SetMaximum(1.2 * maxy_raw);

 raw_P_gen->SetTitle("Q^{2} distribution for Helium-4 as a target (raw);Q^{2} [GeV^{2}];Counts");

 raw_P_gen->Draw("HIST");
 raw_P_wgt->Draw("HIST SAME");
 raw_H_gen->Draw("HIST SAME");
 raw_H_wgt->Draw("HIST SAME");

 TLegend *leg_raw = new TLegend(0.60, 0.65, 0.88, 0.88);
 leg_raw->AddEntry(raw_P_gen, "Proton generated", "l");
 leg_raw->AddEntry(raw_P_wgt, "Proton weighted", "l");
 leg_raw->AddEntry(raw_H_gen, "^{4}He generated", "l");
 leg_raw->AddEntry(raw_H_wgt, "^{4}He weighted", "l");
 leg_raw->Draw();

 c_raw->Update();
 
 
 TCanvas *c_norm = new TCanvas("c_norm", "Q2 comparison normalized", 1000, 800);

 norm_P_gen->SetLineColor(kBlack);
 norm_P_wgt->SetLineColor(kRed);
 norm_H_gen->SetLineColor(kBlue);
 norm_H_wgt->SetLineColor(kGreen+2);

 norm_P_gen->SetLineWidth(2);
 norm_P_wgt->SetLineWidth(2);
 norm_H_gen->SetLineWidth(2);
 norm_H_wgt->SetLineWidth(2);

 double maxy_norm = std::max(std::max(norm_P_gen->GetMaximum(), norm_P_wgt->GetMaximum()),std::max(norm_H_gen->GetMaximum(), norm_H_wgt->GetMaximum()));
 norm_P_gen->SetMaximum(1.2 * maxy_norm);

 norm_P_gen->SetTitle("Q^{2} distribution for Helium-4 as a target (normalized);Q^{2} [GeV^{2}];Normalized counts");

 norm_P_gen->Draw("HIST");
 norm_P_wgt->Draw("HIST SAME");
 norm_H_gen->Draw("HIST SAME");
 norm_H_wgt->Draw("HIST SAME");

 TLegend *leg_norm = new TLegend(0.60, 0.65, 0.88, 0.88);
 leg_norm->AddEntry(norm_P_gen, "Proton generated", "l");
 leg_norm->AddEntry(norm_P_wgt, "Proton weighted", "l");
 leg_norm->AddEntry(norm_H_gen, "^{4}He generated", "l");
 leg_norm->AddEntry(norm_H_wgt, "^{4}He weighted", "l");
 leg_norm->Draw();

 c_norm->Update();
 
 
 
 
 
 //distribution of W including normalized ones
 
 
 // ---------------------------
 // Clone W histograms for raw
 // ---------------------------
 TH1F *raw_wP_gen = (TH1F*)wP_generated->Clone("raw_wP_gen");
 TH1F *raw_wP_wgt = (TH1F*)wP_weighted->Clone("raw_wP_wgt");
 TH1F *raw_wH_gen = (TH1F*)wH_generated->Clone("raw_wH_gen");
 TH1F *raw_wH_wgt = (TH1F*)wH_weighted->Clone("raw_wH_wgt");

 // ----------------------------------
 // Clone W histograms for normalized
 // ----------------------------------
 TH1F *norm_wP_gen = (TH1F*)wP_generated->Clone("norm_wP_gen");
 TH1F *norm_wP_wgt = (TH1F*)wP_weighted->Clone("norm_wP_wgt");
 TH1F *norm_wH_gen = (TH1F*)wH_generated->Clone("norm_wH_gen");
 TH1F *norm_wH_wgt = (TH1F*)wH_weighted->Clone("norm_wH_wgt");

 // optional but useful
 raw_wP_gen->SetDirectory(0);
 raw_wP_wgt->SetDirectory(0);
 raw_wH_gen->SetDirectory(0);
 raw_wH_wgt->SetDirectory(0);

 norm_wP_gen->SetDirectory(0);
 norm_wP_wgt->SetDirectory(0);
 norm_wH_gen->SetDirectory(0);
 norm_wH_wgt->SetDirectory(0);

 // normalize only the normalized copies
 if (norm_wP_gen->Integral() > 0) norm_wP_gen->Scale(1.0 / norm_wP_gen->Integral());
 if (norm_wP_wgt->Integral() > 0) norm_wP_wgt->Scale(1.0 / norm_wP_wgt->Integral());
 if (norm_wH_gen->Integral() > 0) norm_wH_gen->Scale(1.0 / norm_wH_gen->Integral());
 if (norm_wH_wgt->Integral() > 0) norm_wH_wgt->Scale(1.0 / norm_wH_wgt->Integral());




 TCanvas *c_W = new TCanvas("c_W", "W distributions", 1800, 800);
 c_W->Divide(2,1);

 gStyle->SetOptStat(0);

 // =====================
 // Left pad: raw
 // =====================
 c_W->cd(1);
 // gPad->SetLogy();   // uncomment if needed
 // gPad->SetLogx();   // uncomment if needed

 raw_wP_gen->SetLineColor(kBlack);
 raw_wP_wgt->SetLineColor(kRed);
 raw_wH_gen->SetLineColor(kBlue);
 raw_wH_wgt->SetLineColor(kGreen+2);

 raw_wP_gen->SetLineWidth(2);
 raw_wP_wgt->SetLineWidth(2);
 raw_wH_gen->SetLineWidth(2);
 raw_wH_wgt->SetLineWidth(2);

 double maxy_W_raw = std::max(
     std::max(raw_wP_gen->GetMaximum(), raw_wP_wgt->GetMaximum()),
     std::max(raw_wH_gen->GetMaximum(), raw_wH_wgt->GetMaximum())
 );

 raw_wP_gen->SetMaximum(1.2 * maxy_W_raw);
 raw_wP_gen->SetMinimum(0.0);
 raw_wP_gen->SetTitle("W distribution (raw);W [GeV];Counts");

 raw_wP_gen->Draw("HIST");
 raw_wP_wgt->Draw("HIST SAME");
 raw_wH_gen->Draw("HIST SAME");
 raw_wH_wgt->Draw("HIST SAME");

 TLegend *leg_W_raw = new TLegend(0.58, 0.65, 0.88, 0.88);
 leg_W_raw->AddEntry(raw_wP_gen, "Proton generated", "l");
 leg_W_raw->AddEntry(raw_wP_wgt, "Proton weighted", "l");
 leg_W_raw->AddEntry(raw_wH_gen, "^{4}He generated", "l");
 leg_W_raw->AddEntry(raw_wH_wgt, "^{4}He weighted", "l");
 leg_W_raw->Draw();

 // ==========================
 // Right pad: normalized
 // ==========================
 c_W->cd(2);
 // gPad->SetLogy();   // uncomment if needed
 // gPad->SetLogx();   // uncomment if needed

 norm_wP_gen->SetLineColor(kBlack);
 norm_wP_wgt->SetLineColor(kRed);
 norm_wH_gen->SetLineColor(kBlue);
 norm_wH_wgt->SetLineColor(kGreen+2);

 norm_wP_gen->SetLineWidth(2);
 norm_wP_wgt->SetLineWidth(2);
 norm_wH_gen->SetLineWidth(2);
 norm_wH_wgt->SetLineWidth(2);

 double maxy_W_norm = std::max(
     std::max(norm_wP_gen->GetMaximum(), norm_wP_wgt->GetMaximum()),
     std::max(norm_wH_gen->GetMaximum(), norm_wH_wgt->GetMaximum())
 );

 norm_wP_gen->SetMaximum(1.2 * maxy_W_norm);
 norm_wP_gen->SetMinimum(0.0);
 norm_wP_gen->SetTitle("W distribution (normalized);W [GeV];Normalized counts");

 norm_wP_gen->Draw("HIST");
 norm_wP_wgt->Draw("HIST SAME");
 norm_wH_gen->Draw("HIST SAME");
 norm_wH_wgt->Draw("HIST SAME");

 TLegend *leg_W_norm = new TLegend(0.58, 0.65, 0.88, 0.88);
 leg_W_norm->AddEntry(norm_wP_gen, "Proton generated", "l");
 leg_W_norm->AddEntry(norm_wP_wgt, "Proton weighted", "l");
 leg_W_norm->AddEntry(norm_wH_gen, "^{4}He generated", "l");
 leg_W_norm->AddEntry(norm_wH_wgt, "^{4}He weighted", "l");
 leg_W_norm->Draw();

 c_W->Update();













 //generated and weighted histogram for XB for proton and helium-4 as target and normalized them
 
 
 // ---------------------------------
 // Clone x_B histograms for raw plot
 // ---------------------------------
 TH1F *raw_xp_gen = (TH1F*)hxp_generated->Clone("raw_xp_gen");
 TH1F *raw_xp_wgt = (TH1F*)hxp_weighted->Clone("raw_xp_wgt");
 TH1F *raw_xh_gen = (TH1F*)hxh_generated->Clone("raw_xh_gen");
 TH1F *raw_xh_wgt = (TH1F*)hxh_weighted->Clone("raw_xh_wgt");

 // ----------------------------------------
 // Clone x_B histograms for normalized plot
 // ----------------------------------------
 TH1F *norm_xp_gen = (TH1F*)hxp_generated->Clone("norm_xp_gen");
 TH1F *norm_xp_wgt = (TH1F*)hxp_weighted->Clone("norm_xp_wgt");
 TH1F *norm_xh_gen = (TH1F*)hxh_generated->Clone("norm_xh_gen");
 TH1F *norm_xh_wgt = (TH1F*)hxh_weighted->Clone("norm_xh_wgt");

 // normalize only the normalized copies
 if (norm_xp_gen->Integral() > 0) norm_xp_gen->Scale(1.0 / norm_xp_gen->Integral());
 if (norm_xp_wgt->Integral() > 0) norm_xp_wgt->Scale(1.0 / norm_xp_wgt->Integral());
 if (norm_xh_gen->Integral() > 0) norm_xh_gen->Scale(1.0 / norm_xh_gen->Integral());
 if (norm_xh_wgt->Integral() > 0) norm_xh_wgt->Scale(1.0 / norm_xh_wgt->Integral());



 TCanvas *c_xB = new TCanvas("c_xB", "x_B comparison", 1800, 800);
 c_xB->Divide(2,1);

 gStyle->SetOptStat(0);

 // =====================
 // Left pad: raw x_B
 // =====================
 c_xB->cd(1);

 raw_xp_gen->SetLineColor(kBlack);
 raw_xp_wgt->SetLineColor(kRed);
 raw_xh_gen->SetLineColor(kBlue);
 raw_xh_wgt->SetLineColor(kGreen+2);

 raw_xp_gen->SetLineWidth(2);
 raw_xp_wgt->SetLineWidth(2);
 raw_xh_gen->SetLineWidth(2);
 raw_xh_wgt->SetLineWidth(2);

 double maxy_x_raw = std::max(
     std::max(raw_xp_gen->GetMaximum(), raw_xp_wgt->GetMaximum()),
     std::max(raw_xh_gen->GetMaximum(), raw_xh_wgt->GetMaximum())
 );
 raw_xp_gen->SetMaximum(1.2 * maxy_x_raw);


 raw_xp_gen->SetTitle("x_{B} distribution (raw);x_{B};Counts");
 raw_xp_gen->Draw("HIST");
 raw_xp_wgt->Draw("HIST SAME");
 raw_xh_gen->Draw("HIST SAME");
 raw_xh_wgt->Draw("HIST SAME");

 TLegend *leg_x_raw = new TLegend(0.58, 0.65, 0.88, 0.88);
 leg_x_raw->AddEntry(raw_xp_gen, "Proton generated", "l");
 leg_x_raw->AddEntry(raw_xp_wgt, "Proton weighted", "l");
 leg_x_raw->AddEntry(raw_xh_gen, "^{4}He generated", "l");
 leg_x_raw->AddEntry(raw_xh_wgt, "^{4}He weighted", "l");
 leg_x_raw->Draw();

 // ==========================
 // Right pad: normalized x_B
 // ==========================
 c_xB->cd(2);

 norm_xp_gen->SetLineColor(kBlack);
 norm_xp_wgt->SetLineColor(kRed);
 norm_xh_gen->SetLineColor(kBlue);
 norm_xh_wgt->SetLineColor(kGreen+2);

 norm_xp_gen->SetLineWidth(2);
 norm_xp_wgt->SetLineWidth(2);
 norm_xh_gen->SetLineWidth(2);
 norm_xh_wgt->SetLineWidth(2);

 double maxy_x_norm = std::max(
     std::max(norm_xp_gen->GetMaximum(), norm_xp_wgt->GetMaximum()),
     std::max(norm_xh_gen->GetMaximum(), norm_xh_wgt->GetMaximum())
 );
 //norm_xp_gen->SetMaximum(1.2 * maxy_x_norm);
 
 norm_xp_gen->SetMinimum(0.0);
 norm_xp_gen->SetMaximum(0.06);

 norm_xp_gen->SetTitle("x_{B} distribution (normalized);x_{B};Normalized counts");
 norm_xp_gen->Draw("HIST");
 norm_xp_wgt->Draw("HIST SAME");
 norm_xh_gen->Draw("HIST SAME");
 norm_xh_wgt->Draw("HIST SAME");

 TLegend *leg_x_norm = new TLegend(0.58, 0.65, 0.88, 0.88);
 leg_x_norm->AddEntry(norm_xp_gen, "Proton generated", "l");
 leg_x_norm->AddEntry(norm_xp_wgt, "Proton weighted", "l");
 leg_x_norm->AddEntry(norm_xh_gen, "^{4}He generated", "l");
 leg_x_norm->AddEntry(norm_xh_wgt, "^{4}He weighted", "l");
 leg_x_norm->Draw();

 c_xB->Update();




 //for t_eff._generated and t_eff_weighted and t_h_generated and t_h_weighted
 
 
 // raw clones
 TH1F *raw_tEff_gen = (TH1F*)htEff_generated->Clone("raw_tEff_gen");
 TH1F *raw_tEff_wgt = (TH1F*)htEff_weighted->Clone("raw_tEff_wgt");
 TH1F *raw_tH_gen   = (TH1F*)htH_generated->Clone("raw_tH_gen");
 TH1F *raw_tH_wgt   = (TH1F*)htH_weighted->Clone("raw_tH_wgt");

 // normalized clones
 TH1F *norm_tEff_gen = (TH1F*)htEff_generated->Clone("norm_tEff_gen");
 TH1F *norm_tEff_wgt = (TH1F*)htEff_weighted->Clone("norm_tEff_wgt");
 TH1F *norm_tH_gen   = (TH1F*)htH_generated->Clone("norm_tH_gen");
 TH1F *norm_tH_wgt   = (TH1F*)htH_weighted->Clone("norm_tH_wgt");

 if (norm_tEff_gen->Integral() > 0) norm_tEff_gen->Scale(1.0 / norm_tEff_gen->Integral());
 if (norm_tEff_wgt->Integral() > 0) norm_tEff_wgt->Scale(1.0 / norm_tEff_wgt->Integral());
 if (norm_tH_gen->Integral() > 0)   norm_tH_gen->Scale(1.0 / norm_tH_gen->Integral());
 if (norm_tH_wgt->Integral() > 0)   norm_tH_wgt->Scale(1.0 / norm_tH_wgt->Integral());




 TCanvas *c_t = new TCanvas("c_t", "t comparison", 1800, 800);
 c_t->Divide(2,1);

 gStyle->SetOptStat(0);

 // =====================
 // Left pad: raw
 // =====================
 c_t->cd(1);

 raw_tEff_gen->SetLineColor(kBlack);
 raw_tEff_wgt->SetLineColor(kRed);
 raw_tH_gen->SetLineColor(kBlue);
 raw_tH_wgt->SetLineColor(kGreen+2);

 raw_tEff_gen->SetLineWidth(2);
 raw_tEff_wgt->SetLineWidth(2);
 raw_tH_gen->SetLineWidth(2);
 raw_tH_wgt->SetLineWidth(2);

 double maxy_t_raw = std::max(
     std::max(raw_tEff_gen->GetMaximum(), raw_tEff_wgt->GetMaximum()),
     std::max(raw_tH_gen->GetMaximum(), raw_tH_wgt->GetMaximum())
 );
 raw_tEff_gen->SetMaximum(1.2 * maxy_t_raw);

 raw_tEff_gen->SetTitle("t distribution (raw);t [GeV^{2}];Counts");
 raw_tEff_gen->Draw("HIST");
 raw_tEff_wgt->Draw("HIST SAME");
 raw_tH_gen->Draw("HIST SAME");
 raw_tH_wgt->Draw("HIST SAME");

 TLegend *leg_t_raw = new TLegend(0.58, 0.65, 0.88, 0.88);
 leg_t_raw->AddEntry(raw_tEff_gen, "t_{eff} generated", "l");
 leg_t_raw->AddEntry(raw_tEff_wgt, "t_{eff} weighted", "l");
 leg_t_raw->AddEntry(raw_tH_gen,   "t_{H} generated", "l");
 leg_t_raw->AddEntry(raw_tH_wgt,   "t_{H} weighted", "l");
 leg_t_raw->Draw();

 // =====================
 // Right pad: normalized
 // =====================
 c_t->cd(2);

 norm_tEff_gen->SetLineColor(kBlack);
 norm_tEff_wgt->SetLineColor(kRed);
 norm_tH_gen->SetLineColor(kBlue);
 norm_tH_wgt->SetLineColor(kGreen+2);

 norm_tEff_gen->SetLineWidth(2);
 norm_tEff_wgt->SetLineWidth(2);
 norm_tH_gen->SetLineWidth(2);
 norm_tH_wgt->SetLineWidth(2);

 double maxy_t_norm = std::max(
     std::max(norm_tEff_gen->GetMaximum(), norm_tEff_wgt->GetMaximum()),
     std::max(norm_tH_gen->GetMaximum(), norm_tH_wgt->GetMaximum())
 );
 norm_tEff_gen->SetMaximum(1.2 * maxy_t_norm);

 norm_tEff_gen->SetTitle("t distribution (normalized);t [GeV^{2}];Normalized counts");
 norm_tEff_gen->Draw("HIST");
 norm_tEff_wgt->Draw("HIST SAME");
 norm_tH_gen->Draw("HIST SAME");
 norm_tH_wgt->Draw("HIST SAME");

 TLegend *leg_t_norm = new TLegend(0.58, 0.65, 0.88, 0.88);
 leg_t_norm->AddEntry(norm_tEff_gen, "t_{eff} generated", "l");
 leg_t_norm->AddEntry(norm_tEff_wgt, "t_{eff} weighted", "l");
 leg_t_norm->AddEntry(norm_tH_gen,   "t_{H} generated", "l");
 leg_t_norm->AddEntry(norm_tH_wgt,   "t_{H} weighted", "l");
 leg_t_norm->Draw();

 c_t->Update();






 //difference of t:
 
 TH1F *raw_tdiff_gen  = (TH1F*)tdiff_generated->Clone("raw_tdiff_gen");
 TH1F *raw_tdiff_wgt  = (TH1F*)tdiff_weighted->Clone("raw_tdiff_wgt");

 TH1F *norm_tdiff_gen = (TH1F*)tdiff_generated->Clone("norm_tdiff_gen");
 TH1F *norm_tdiff_wgt = (TH1F*)tdiff_weighted->Clone("norm_tdiff_wgt");

 if (norm_tdiff_gen->Integral() > 0)
     norm_tdiff_gen->Scale(1.0 / norm_tdiff_gen->Integral());

 if (norm_tdiff_wgt->Integral() > 0)
     norm_tdiff_wgt->Scale(1.0 / norm_tdiff_wgt->Integral());
 
 TCanvas *c_tdiff = new TCanvas("c_tdiff", "tdiff comparison", 1600, 700);
 c_tdiff->Divide(2,1);

 gStyle->SetOptStat(0);

 // raw
 c_tdiff->cd(1);
 gPad->SetLogy();

 raw_tdiff_gen->SetLineColor(kBlue);
 raw_tdiff_wgt->SetLineColor(kGreen+2);

 raw_tdiff_gen->SetLineWidth(2);
 raw_tdiff_wgt->SetLineWidth(2);

 double maxy_tdiff_raw = std::max(raw_tdiff_gen->GetMaximum(), raw_tdiff_wgt->GetMaximum());
 raw_tdiff_gen->SetMaximum(1.2 * maxy_tdiff_raw);
 raw_tdiff_gen->SetMinimum(1.0);

 raw_tdiff_gen->SetTitle("t_{max}-t distribution for ^{4}He (raw);t_{max}-t [GeV^{2}];Counts");
 raw_tdiff_gen->Draw("HIST");
 raw_tdiff_wgt->Draw("HIST SAME");

 TLegend *leg_tdiff_raw = new TLegend(0.58, 0.72, 0.88, 0.88);
 leg_tdiff_raw->AddEntry(raw_tdiff_gen, "t_{max}-t generated", "l");
 leg_tdiff_raw->AddEntry(raw_tdiff_wgt, "t_{max}-t weighted", "l");
 leg_tdiff_raw->Draw();

 // normalized
 c_tdiff->cd(2);

 norm_tdiff_gen->SetLineColor(kBlue);
 norm_tdiff_wgt->SetLineColor(kGreen+2);

 norm_tdiff_gen->SetLineWidth(2);
 norm_tdiff_wgt->SetLineWidth(2);

 double maxy_tdiff_norm = std::max(norm_tdiff_gen->GetMaximum(), norm_tdiff_wgt->GetMaximum());
 norm_tdiff_gen->SetMaximum(1.2 * maxy_tdiff_norm);

 norm_tdiff_gen->SetTitle("t_{max}-t distribution for ^{4}He (normalized);t_{max}-t [GeV^{2}];Normalized counts");
 norm_tdiff_gen->Draw("HIST");
 norm_tdiff_wgt->Draw("HIST SAME");

 TLegend *leg_tdiff_norm = new TLegend(0.58, 0.72, 0.88, 0.88);
 leg_tdiff_norm->AddEntry(norm_tdiff_gen, "t_{min}-t generated", "l");
 leg_tdiff_norm->AddEntry(norm_tdiff_wgt, "t_{min}-t weighted", "l");
 leg_tdiff_norm->Draw();

 c_tdiff->Update();




 //test of weighted tmax-t = tdiff
 
 TCanvas *ctest = new TCanvas("ctest","tdiff weighted",800,600);
 tdiff_weighted->SetLineColor(kGreen+2);
 tdiff_weighted->SetLineWidth(2);
 tdiff_weighted->Draw("HIST");
 ctest->Update();



 //two-dimensional coverage of kinematic variable
 
 TCanvas *can2all = new TCanvas("can2all","Q2 phase-space plots",1800,600);
 can2all->Divide(3,1);

 gStyle->SetOptStat(0);
 

 // pad 1
 can2all->cd(1);
 gPad->SetLogz();
 gPad->SetRightMargin(0.15);
 Q2HW->SetTitle("W_{H} vs Q^{2}_{H};Q^{2}_{H} [GeV^{2}];W_{H} [GeV]");
 Q2HW->Draw("colz");

 // pad 2
 can2all->cd(2);
 gPad->SetLogz();
 gPad->SetRightMargin(0.15);
 Q2HXB->SetTitle("x_{B} vs Q^{2}_{H};x_{B};Q^{2}_{H} [GeV^{2}]");
 Q2HXB->Draw("colz");

 // pad 3
 can2all->cd(3);
 gPad->SetLogz();
 gPad->SetRightMargin(0.15);
 Q2HtH->SetTitle("-t vs Q^{2}_{H};Q^{2}_{H} [GeV^{2}];t_{min}-t [GeV^{2}]");
 Q2HtH->SetMinimum(1e-2);
 Q2HtH->Draw("colz");

 can2all->Update();
 
 
 
 //smearing in 2D histogram
 
 TCanvas *c_rec = new TCanvas("c_rec", "Reco weighted kinematics", 1400, 600);
 c_rec->Divide(2,1);

 c_rec->cd(1);
 gPad->SetLogz();
 Q2HW_rec->Draw("COLZ");

 c_rec->cd(2);
 gPad->SetLogz();
 Q2HXB_rec->Draw("COLZ");
 
 
 
 
 
 
 TProfile *raw_FCsq = (TProfile*)pFCsq_t_h->Clone("raw_FCsq");
 TProfile *norm_FCsq = (TProfile*)pFCsq_t_h->Clone("norm_FCsq");
 
 if (norm_FCsq->Integral() > 0)
	 norm_FCsq->Scale(1.0 / norm_FCsq->Integral());
 
 
 
 
 TCanvas *c_FC = new TCanvas("c_FC", "F_{C}^{2} vs -t", 1600, 700);
 c_FC->Divide(2,1);

 gStyle->SetOptStat(0);

 // =====================
 // Left pad: raw
 // =====================
 c_FC->cd(1);
 gPad->SetLogx();
 gPad->SetLogy();
 gPad->SetGrid();

 raw_FCsq->SetLineColor(kBlue);
 raw_FCsq->SetMarkerColor(kBlue);
 raw_FCsq->SetMarkerStyle(20);
 raw_FCsq->SetMarkerSize(0.9);
 raw_FCsq->SetLineWidth(2);

 raw_FCsq->SetMinimum(1e-8);
 raw_FCsq->SetTitle("F_{C,^{4}He}^{2} vs -t (raw);-t [GeV^{2}];F_{C,^{4}He}^{2}");
 raw_FCsq->Draw("PL");

 // ==========================
 // Right pad: normalized
 // ==========================
 c_FC->cd(2);
 gPad->SetLogx();
 gPad->SetLogy();
 gPad->SetGrid();

 norm_FCsq->SetLineColor(kRed);
 norm_FCsq->SetMarkerColor(kRed);
 norm_FCsq->SetMarkerStyle(21);
 norm_FCsq->SetMarkerSize(0.9);
 norm_FCsq->SetLineWidth(2);

 norm_FCsq->SetMinimum(1e-8);
 norm_FCsq->SetTitle("F_{C,^{4}He}^{2} vs -t (normalized);-t [GeV^{2}];Normalized F_{C,^{4}He}^{2}");
 norm_FCsq->Draw("PL");

 c_FC->Update();
 
 

 
 
 
 
 //plot of R for both protn and helium-4 into different pads into a single canvas
 
 TCanvas *c_R = new TCanvas("c_R", "R vs Q2 for proton and He4", 1600, 700);
 c_R->Divide(1,1);

 gStyle->SetOptStat(0);

 // =====================
 // Left pad: Proton
 // =====================
 c_R->cd(1);
 gPad->SetLogx();
 gPad->SetLogy();
 gPad->SetGrid();

 pRP->SetLineColor(kRed);
 pRP->SetMarkerColor(kRed);
 pRP->SetMarkerStyle(20);
 pRP->SetMarkerSize(1.0);
 pRP->SetLineWidth(2);

  //pRP->SetTitle("Proton: R = #sigma_{L}/#sigma_{T} vs Q^{2};Q^{2} [GeV^{2}];R");
 //pRP->SetTitle("Proton: R = #frac{#sigma_{L}}{#sigma_{T}} vs Q^{2};Q^{2} [GeV^{2}];R = #frac{#sigma_{L}}{#sigma_{T}}");
 pRP->SetMinimum(1e-1);
 pRP->GetXaxis()->SetTitle( "Q^{2} [GeV^{2}]");
 pRP->GetXaxis()->CenterTitle(true);
 pRP->GetYaxis()->SetTitle("R = #sigma_{L}/#sigma_{T}");
 pRP->GetYaxis()->CenterTitle(true);
 pRP->Draw("PL");

 //TLegend *legP = new TLegend(0.58, 0.78, 0.88, 0.88);
 //legP->AddEntry(pRP, "Proton", "lp");
 //legP-AddEntry(pRP,"#gamma*p -> #phi +p","lp")
 TLegend *legP = new TLegend(0.55, 0.76, 0.88, 0.88);
 legP->SetBorderSize(0);
 legP->SetFillStyle(0);
 legP->SetTextSize(0.035);

 legP->AddEntry(pRP,"#gamma*p #rightarrow #phi + p","lp");

 legP->Draw();
 legP->Draw();

 /*
 // =====================
 // Right pad: He-4
 // =====================
 c_R->cd(2);
 gPad->SetLogx();
 gPad->SetLogy();
 gPad->SetGrid();

 pRH->SetLineColor(kBlue);
 pRH->SetMarkerColor(kBlue);
 pRH->SetMarkerStyle(21);
 pRH->SetMarkerSize(1.0);
 pRH->SetLineWidth(2);

 pRH->SetTitle("^{4}He: R = #sigma_{L}/#sigma_{T} vs Q^{2};Q^{2} [GeV^{2}];R");
 pRH->SetMinimum(1e-3);
 pRH->Draw("PL");

 TLegend *legH = new TLegend(0.58, 0.78, 0.88, 0.88);
 legH->AddEntry(pRH, "^{4}He", "lp");
 legH->Draw();

 c_R->Update();
 

 TCanvas *c_comp = new TCanvas("c_comp", "Proton vs He-4 at fixed W", 1800, 600);
 c_comp->Divide(3,1);

 gStyle->SetOptStat(0);

 // -------- W ~ 2.5 --------
 c_comp->cd(1);
 gPad->SetLogx();
 gPad->SetLogy();
 gPad->SetGrid();

 pSigmaT_Q2_W45_P->SetLineColor(kRed);
 pSigmaT_Q2_W45_P->SetMarkerColor(kRed);
 pSigmaT_Q2_W45_P->SetMarkerStyle(20);
 pSigmaT_Q2_W45_P->SetLineWidth(2);

 pSigmaT_Q2_W45_H->SetLineColor(kBlue);
 pSigmaT_Q2_W45_H->SetMarkerColor(kBlue);
 pSigmaT_Q2_W45_H->SetMarkerStyle(21);
 pSigmaT_Q2_W45_H->SetLineWidth(2);

 pSigmaT_Q2_W45_P->SetTitle("#sigma_{T} vs Q^{2} at W #approx 2.5 GeV;Q^{2} [GeV^{2}];#sigma_{T} [nb]");
 pSigmaT_Q2_W45_P->SetMinimum(1e-2);
 pSigmaT_Q2_W45_P->Draw();
 pSigmaT_Q2_W45_H->Draw("SAME");

 TLegend *leg45 = new TLegend(0.55, 0.72, 0.88, 0.88);
 leg45->AddEntry(pSigmaT_Q2_W45_P, "Proton", "lp");
 leg45->AddEntry(pSigmaT_Q2_W45_H, "^{4}He", "lp");
 leg45->Draw();

 // -------- W ~ 5.0 --------
 c_comp->cd(2);
 gPad->SetLogx();
 gPad->SetLogy();
 gPad->SetGrid();

 pSigmaT_Q2_W50_P->SetLineColor(kRed);
 pSigmaT_Q2_W50_P->SetMarkerColor(kRed);
 pSigmaT_Q2_W50_P->SetMarkerStyle(20);
 pSigmaT_Q2_W50_P->SetLineWidth(2);

 pSigmaT_Q2_W50_H->SetLineColor(kBlue);
 pSigmaT_Q2_W50_H->SetMarkerColor(kBlue);
 pSigmaT_Q2_W50_H->SetMarkerStyle(21);
 pSigmaT_Q2_W50_H->SetLineWidth(2);

 pSigmaT_Q2_W50_P->SetTitle("#sigma_{T} vs Q^{2} at W #approx 5.0 GeV;Q^{2} [GeV^{2}];#sigma_{T} [nb]");
 pSigmaT_Q2_W50_P->SetMinimum(1e-2);
 pSigmaT_Q2_W50_P->Draw();
 pSigmaT_Q2_W50_H->Draw("SAME");

 TLegend *leg50 = new TLegend(0.55, 0.72, 0.88, 0.88);
 leg50->AddEntry(pSigmaT_Q2_W50_P, "Proton", "lp");
 leg50->AddEntry(pSigmaT_Q2_W50_H, "^{4}He", "lp");
 leg50->Draw();

 // -------- W ~ 7.5 --------
 c_comp->cd(3);
 gPad->SetLogx();
 gPad->SetLogy();
 gPad->SetGrid();

 pSigmaT_Q2_W75_P->SetLineColor(kRed);
 pSigmaT_Q2_W75_P->SetMarkerColor(kRed);
 pSigmaT_Q2_W75_P->SetMarkerStyle(20);
 pSigmaT_Q2_W75_P->SetLineWidth(2);

 pSigmaT_Q2_W75_H->SetLineColor(kBlue);
 pSigmaT_Q2_W75_H->SetMarkerColor(kBlue);
 pSigmaT_Q2_W75_H->SetMarkerStyle(21);
 pSigmaT_Q2_W75_H->SetLineWidth(2);

 pSigmaT_Q2_W75_P->SetTitle("#sigma_{T} vs Q^{2} at W #approx 7.5 GeV;Q^{2} [GeV^{2}];#sigma_{T} [nb]");
 pSigmaT_Q2_W75_P->SetMinimum(1e-2);
 pSigmaT_Q2_W75_P->Draw();
 pSigmaT_Q2_W75_H->Draw("SAME");

 TLegend *leg75 = new TLegend(0.55, 0.72, 0.88, 0.88);
 leg75->AddEntry(pSigmaT_Q2_W75_P, "Proton", "lp");
 leg75->AddEntry(pSigmaT_Q2_W75_H, "^{4}He", "lp");
 leg75->Draw();


 c_comp->Update();
 
 */
 
 
 //gStyle->SetOptStat(0);

 TCanvas *cSigmaT = new TCanvas("cSigmaT", "Proton transverse cross section",900, 70);
 gStyle->SetOptStat(0);
 cSigmaT->cd();

 gPad->SetLogx();
 gPad->SetLogy();
 gPad->SetGrid();

 //pSigmaT_Q2_P->SetTitle("Proton Target");

 pSigmaT_Q2_P->GetXaxis()->SetTitle( "Q^{2} [GeV^{2}]");
 pSigmaT_Q2_P->GetXaxis()->CenterTitle(true);

 pSigmaT_Q2_P->GetYaxis()->SetTitle("#sigma_{T} [nb]");
 pSigmaT_Q2_P->GetYaxis()->CenterTitle(true);

 pSigmaT_Q2_P->GetXaxis()->CenterTitle(true);
 pSigmaT_Q2_P->GetYaxis()->CenterTitle(true);

 pSigmaT_Q2_P->GetXaxis()->SetTitleOffset(1.2);
 pSigmaT_Q2_P->GetYaxis()->SetTitleOffset(1.3);

 pSigmaT_Q2_P->GetXaxis()->SetTitleSize(0.045);
 pSigmaT_Q2_P->GetYaxis()->SetTitleSize(0.045);

 pSigmaT_Q2_P->GetXaxis()->SetLabelSize(0.040);
 pSigmaT_Q2_P->GetYaxis()->SetLabelSize(0.040);

 pSigmaT_Q2_P->SetLineColor(kRed);
 pSigmaT_Q2_P->SetMarkerColor(kRed);
 pSigmaT_Q2_P->SetMarkerStyle(20);
 pSigmaT_Q2_P->SetMarkerSize(0.9);
 pSigmaT_Q2_P->SetLineWidth(2);
 pSigmaT_Q2_P->SetMinimum(1.0e-2);

 pSigmaT_Q2_P->Draw("E1");

 TLegend *legSigmaT = new TLegend(0.55, 0.78, 0.88, 0.88);

 legSigmaT->SetBorderSize(0);
 legSigmaT->SetFillStyle(0);
 legSigmaT->SetTextSize(0.035);

 legSigmaT->AddEntry(pSigmaT_Q2_P,"#gamma*p #rightarrow #phi + p", "lp");

 legSigmaT->Draw();

 cSigmaT->Modified();
 cSigmaT->Update();
 
 
 
 //gStyle->SetOptStat(0);

 TCanvas *cCTW = new TCanvas("cCTW", "Transverse cross section against W",900, 70);
 gStyle->SetOptStat(0);
 cCTW->cd();

 gPad->SetLogx();
 gPad->SetLogy();
 gPad->SetGrid();

 //pSigmaT_Q2_P->SetTitle("Proton Target");

 CTW->GetXaxis()->SetTitle( "W [GeV]");
 CTW->GetXaxis()->CenterTitle(true);

 CTW->GetYaxis()->SetTitle("#sigma_{T} [nb]");
 CTW->GetYaxis()->CenterTitle(true);

 CTW->GetXaxis()->CenterTitle(true);
 CTW->GetYaxis()->CenterTitle(true);

 CTW->GetXaxis()->SetTitleOffset(1.2);
 CTW->GetYaxis()->SetTitleOffset(1.3);

 CTW->GetXaxis()->SetTitleSize(0.045);
 CTW->GetYaxis()->SetTitleSize(0.045);

 CTW->GetXaxis()->SetLabelSize(0.040);
 CTW->GetYaxis()->SetLabelSize(0.040);

 CTW->SetLineColor(kRed);
 CTW->SetMarkerColor(kRed);
 CTW->SetMarkerStyle(20);
 CTW->SetMarkerSize(0.9);
 CTW->SetLineWidth(2);
 CTW->SetMinimum(1.0e-2);

 CTW->Draw("E1");

 TLegend *legCTW = new TLegend(0.55, 0.78, 0.88, 0.88);

 legCTW->SetBorderSize(0);
 legCTW->SetFillStyle(0);
 legCTW->SetTextSize(0.035);

 legCTW->AddEntry(pSigmaT_Q2_P,"#gamma*p #rightarrow #phi + p", "lp");

 legCTW->Draw();

 cCTW->Modified();
 cCTW->Update();
 
 
 
 
 
 TCanvas *c1_tdiff = new TCanvas("c1_tdiff","Recoil helium kinematics",1600, 1200);

 c1_tdiff->Divide(2,2);
 
 
 // Recoil He-4 momentum
 c1_tdiff->cd(3);
 gPad->SetLogy();

 pHe_generated->SetLineColor(kBlue);
 pHe_weighted->SetLineColor(kGreen+2);

 pHe_generated->SetLineWidth(2);
 pHe_weighted->SetLineWidth(2);

 double maxPHe = std::max(pHe_generated->GetMaximum(),pHe_weighted->GetMaximum());

 pHe_generated->SetMaximum(1.2 * maxPHe);
 pHe_generated->SetMinimum(1.0);

 pHe_generated->SetTitle("Recoil ^{4}He momentum;""p_{He} [GeV/c];Counts");

 pHe_generated->Draw("HIST");
 pHe_weighted->Draw("HIST SAME");

 TLegend *leg_pHe = new TLegend(0.58, 0.72, 0.88, 0.88);
 leg_pHe->AddEntry(pHe_generated, "Generated", "l");
 leg_pHe->AddEntry(pHe_weighted,"Weighted by F_{C}^{2}","l");
 leg_pHe->Draw();


 // Recoil He-4 polar angle
 c1_tdiff->cd(4);
 gPad->SetLogy();

 thetaHe_generated->SetLineColor(kBlue);
 thetaHe_weighted->SetLineColor(kGreen+2);

 thetaHe_generated->SetLineWidth(2);
 thetaHe_weighted->SetLineWidth(2);

 double maxThetaHe = std::max(thetaHe_generated->GetMaximum(),thetaHe_weighted->GetMaximum());

 thetaHe_generated->SetMaximum(1.2 * maxThetaHe);
 thetaHe_generated->SetMinimum(1.0);

 thetaHe_generated->SetTitle("Recoil ^{4}He polar angle;""#theta_{He} [degrees];Counts");

 thetaHe_generated->Draw("HIST");
 thetaHe_weighted->Draw("HIST SAME");

 TLegend *leg_thetaHe =new TLegend(0.58, 0.72, 0.88, 0.88);

 leg_thetaHe->AddEntry(thetaHe_generated,"Generated","l");

 leg_thetaHe->AddEntry(thetaHe_weighted,"Weighted by F_{C}^{2}","l");

 leg_thetaHe->Draw();

 c1_tdiff->Update();
 
 /*
 
 TCanvas *c_t_theta1 = new TCanvas("c_t_theta1","Momentum transfer vs recoil helium angle",1400, 600);

 c_t_theta1->Divide(2, 1);

 // t - t_min vs theta_He
 c_t_theta1->cd(1);
 gPad->SetRightMargin(0.15);
 gPad->SetGrid();

 tdiff_vs_thetaHe->SetStats(0);
 tdiff_vs_thetaHe->Draw("COLZ");

 // t_h vs theta_He
 c_t_theta1->cd(2);
 gPad->SetRightMargin(0.15);
 gPad->SetGrid();

 tH_vs_thetaHe->SetStats(0);
 tH_vs_thetaHe->Draw("COLZ");

 c_t_theta1->Update();
 c_t_theta1->Write();
 c_t_theta1->SaveAs("t_vs_thetaHe_unweighted.png");
 */
 
 TCanvas *c_recoil_corr = new TCanvas("c_recoil_corr","Generated recoil He-4 correlations",1800,600);

 c_recoil_corr->Divide(3, 1);

 // 1. Recoil momentum versus recoil theta
 c_recoil_corr->cd(1);
 gPad->SetRightMargin(0.15);
 gPad->SetGrid();
 h_p_vs_theta->Draw("COLZ");

 // 2. Recoil theta versus vertex z
 c_recoil_corr->cd(2);
 gPad->SetRightMargin(0.15);
 gPad->SetGrid();
 h_theta_vs_vz->Draw("COLZ");

 // 3. Recoil momentum versus vertex z
 c_recoil_corr->cd(3);
 gPad->SetRightMargin(0.15);
 gPad->SetGrid();
 h_p_vs_vz->Draw("COLZ");

 c_recoil_corr->Update();
 c_recoil_corr->SaveAs("generated_recoil_He4_correlations.png");
 
 
 //kaon's theta and momentum TCanvas plot
 
 TCanvas *cKaons = new TCanvas("cKaons","Generated kaon laboratory kinematics",1400, 1000);

 cKaons->Divide(2, 2);

 cKaons->cd(1);
 pKplus_generated->SetLineColor(kRed+1);
 pKplus_generated->SetLineWidth(2);
 pKplus_generated->Draw("HIST");

 cKaons->cd(2);
 pKminus_generated->SetLineColor(kBlue+1);
 pKminus_generated->SetLineWidth(2);
 pKminus_generated->Draw("HIST");

 cKaons->cd(3);
 thetaKplus_generated->SetLineColor(kRed+1);
 thetaKplus_generated->SetLineWidth(2);
 thetaKplus_generated->Draw("HIST");

 cKaons->cd(4);
 thetaKminus_generated->SetLineColor(kBlue+1);
 thetaKminus_generated->SetLineWidth(2);
 thetaKminus_generated->Draw("HIST");

 //cKaons->SaveAs("generated_kaon_momentum_angles.png");
 
 
 
 TCanvas *ceplusKaon2D = new TCanvas("ceplusKaon2D","Generated kaon momentum versus angle",1400, 600);

 ceplusKaon2D->Divide(3, 1);
 
 ceplusKaon2D->cd(1);
 gPad->SetRightMargin(0.15);
 pVsThetaElec->SetStats(0);
 pVsThetaElec->Draw("COLZ");

 ceplusKaon2D->cd(2);
 gPad->SetRightMargin(0.15);
 pVsThetaKplus->SetStats(0);
 pVsThetaKplus->Draw("COLZ");

 ceplusKaon2D->cd(3);
 gPad->SetRightMargin(0.15);
 pVsThetaKminus->SetStats(0);
 pVsThetaKminus->Draw("COLZ");
 
 

 //cKaon2D->SaveAs("generated_kaon_p_vs_theta.png");
 
 TGMainFrame *mainFrame = new TGMainFrame(gClient->GetRoot(), 1400, 1000);
 mainFrame->SetWindowName( "^{4}He momentum-transfer distributions before and after t' cut");

 TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas("embedded_t_summary", mainFrame, 1400, 1000);
 mainFrame->AddFrame(embeddedCanvas,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

 auto *canvas = embeddedCanvas->GetCanvas();
 canvas->Divide(2, 2);

 canvas->cd(1);
 gPad->SetGrid();
 hMinusTBefore->Draw("HIST");

 canvas->cd(2);
 gPad->SetGrid();
 hTprimeBefore->Draw("HIST");

 canvas->cd(3);
 gPad->SetGrid();
 hMinusTAfter->Draw("HIST");

 canvas->cd(4);
 gPad->SetGrid();
 hTprimeAfter->Draw("HIST");

 canvas->Modified();
 canvas->Update();

 mainFrame->MapSubwindows();
 mainFrame->Resize(mainFrame->GetDefaultSize());
 mainFrame->MapWindow();

 std::cout << "Events before t' cut: "
           << hTprimeBefore->GetEntries() << std::endl;
 std::cout << "Events after 0 < t' < 1.0 GeV^2: "
           << hTprimeAfter->GetEntries() << std::endl;
 
 
 TCanvas *cQ2Ht_h = new TCanvas("cQ2Ht_h","",1400, 600);
 gStyle->SetOptStat(0);
 cQ2Ht_h->Divide(1, 1);

 cQ2Ht_h->cd(1);
 gPad->SetRightMargin(0.15);
 //cQ2Ht_h->SetStats(0);
 Q2Hvst_h->Draw("COLZ");
 
 
 
 TCanvas *c_minus_tmin = new TCanvas("c_minus_tmin","Generated -tmin distribution",1000, 800);

 c_minus_tmin->SetLeftMargin(0.13);
 c_minus_tmin->SetRightMargin(0.05);
 c_minus_tmin->SetBottomMargin(0.12);
 c_minus_tmin->SetTopMargin(0.10);

 h_minus_tmin->SetStats(0);
 h_minus_tmin->SetLineColor(kBlue + 1);
 h_minus_tmin->SetLineWidth(3);
 h_minus_tmin->SetFillColorAlpha(kAzure + 1, 0.30);

 h_minus_tmin->GetXaxis()->CenterTitle();
 h_minus_tmin->GetYaxis()->CenterTitle();

 h_minus_tmin->GetXaxis()->SetTitleSize(0.045);
 h_minus_tmin->GetYaxis()->SetTitleSize(0.045);

 h_minus_tmin->GetXaxis()->SetLabelSize(0.040);
 h_minus_tmin->GetYaxis()->SetLabelSize(0.040);

 h_minus_tmin->GetXaxis()->SetTitleOffset(1.10);
 h_minus_tmin->GetYaxis()->SetTitleOffset(1.35);

 h_minus_tmin->Draw("HIST");

 c_minus_tmin->Modified();
 c_minus_tmin->Update();

 c_minus_tmin->SaveAs("minus_tmin_distribution_before_cuts.png");
 
 
 
 
 
 TCanvas *c_tmin_vs_th = new TCanvas( "c_tmin_vs_th","-tmin vs -th before cuts",1000, 850);

 c_tmin_vs_th->SetLeftMargin(0.12);
 c_tmin_vs_th->SetRightMargin(0.15);
 c_tmin_vs_th->SetBottomMargin(0.12);
 c_tmin_vs_th->SetTopMargin(0.10);

 h_tmin_vs_th->SetStats(0);

 h_tmin_vs_th->GetXaxis()->CenterTitle();
 h_tmin_vs_th->GetYaxis()->CenterTitle();

 h_tmin_vs_th->GetXaxis()->SetTitleSize(0.045);
 h_tmin_vs_th->GetYaxis()->SetTitleSize(0.045);

 h_tmin_vs_th->GetXaxis()->SetLabelSize(0.040);
 h_tmin_vs_th->GetYaxis()->SetLabelSize(0.040);

 h_tmin_vs_th->GetXaxis()->SetTitleOffset(1.15);
 h_tmin_vs_th->GetYaxis()->SetTitleOffset(1.20);

 h_tmin_vs_th->Draw("COLZ");

 gPad->Update();

 double axisMin = 0.0;
 double axisMax = 5.0;

 // t' = 0 boundary: -t_h = -t_min
 TLine *lineTdiff0 = new TLine(axisMin, axisMin,axisMax, axisMax);

 lineTdiff0->SetLineColor(kRed + 1);
 lineTdiff0->SetLineWidth(3);
 lineTdiff0->SetLineStyle(2);
 lineTdiff0->Draw("SAME");

 // t' = 1 boundary: -t_h = -t_min + 1
 TLine *lineTdiff1 = new TLine(axisMin, axisMin + 1.0,axisMax - 1.0, axisMax);

 lineTdiff1->SetLineColor(kMagenta + 1);
 lineTdiff1->SetLineWidth(3);
 lineTdiff1->SetLineStyle(2);
 lineTdiff1->Draw("SAME");

 c_tmin_vs_th->Modified();
 c_tmin_vs_th->Update();

 c_tmin_vs_th->SaveAs( "minus_tmin_vs_minus_th_before_cuts.png");
 
 
 
 
 
 
 
 
 
 
 TCanvas *c_tmin_theta = new TCanvas("c_tmin_theta","-tmin vs recoil helium angle",1000, 800);

 c_tmin_theta->SetLeftMargin(0.12);
 c_tmin_theta->SetRightMargin(0.15);
 c_tmin_theta->SetBottomMargin(0.12);

 h_tmin_vs_recoil_theta->SetStats(0);

 h_tmin_vs_recoil_theta->GetXaxis()->CenterTitle();
 h_tmin_vs_recoil_theta->GetYaxis()->CenterTitle();

 h_tmin_vs_recoil_theta->GetXaxis()->SetTitleSize(0.045);
 h_tmin_vs_recoil_theta->GetYaxis()->SetTitleSize(0.045);
 h_tmin_vs_recoil_theta->GetYaxis()->SetRangeUser(0.0, 90.0);

 h_tmin_vs_recoil_theta->GetXaxis()->SetLabelSize(0.040);
 h_tmin_vs_recoil_theta->GetYaxis()->SetLabelSize(0.040);

 // Draw the histogram first
 h_tmin_vs_recoil_theta->Draw("COLZ");

 gPad->Update();

 // Obtain the displayed x-axis range
 double xMin = h_tmin_vs_recoil_theta->GetXaxis()->GetXmin();

 double xMax = h_tmin_vs_recoil_theta->GetXaxis()->GetXmax();

 // Lower ALERT angular boundary
 TLine *line45 = new TLine(xMin, 45.0,xMax, 45.0);

 line45->SetLineColor(kRed + 1);
 line45->SetLineStyle(2);
 line45->SetLineWidth(3);
 //line45->Draw("SAME");

 // Upper ALERT angular boundary
 TLine *line90 = new TLine(xMin, 90.0,xMax, 90.0);

 line90->SetLineColor(kRed + 1);
 line90->SetLineStyle(2);
 line90->SetLineWidth(3);
 //line90->Draw("SAME");


 TBox *acceptedRegion = new TBox(xMin, 45.0, xMax, 90.0);
 acceptedRegion->SetFillColorAlpha(kRed, 0.10);
 acceptedRegion->SetLineColor(kRed + 1);
 acceptedRegion->SetLineStyle(2);
 acceptedRegion->Draw("SAME");

 h_tmin_vs_recoil_theta->Draw("COLZ SAME");
 line45->Draw("SAME");
 line90->Draw("SAME");




 c_tmin_theta->Modified();
 c_tmin_theta->Update();

 c_tmin_theta->SaveAs( "minus_tmin_vs_recoil_He4_theta_before_cuts.png");
 
 
 
 TCanvas *c_t_theta = new TCanvas("c_t_theta","-t vs recoil helium angle",1000, 800);

 c_t_theta->SetLeftMargin(0.12);
 c_t_theta->SetRightMargin(0.15);  // space for color scale
 c_t_theta->SetBottomMargin(0.12);
 c_t_theta->SetTopMargin(0.10);

 h_t_vs_recoil_theta->SetStats(0);

 h_t_vs_recoil_theta->GetXaxis()->CenterTitle();
 h_t_vs_recoil_theta->GetYaxis()->CenterTitle();

 h_t_vs_recoil_theta->GetXaxis()->SetTitleSize(0.045);
 h_t_vs_recoil_theta->GetYaxis()->SetTitleSize(0.045);

 h_t_vs_recoil_theta->GetXaxis()->SetLabelSize(0.040);
 h_t_vs_recoil_theta->GetYaxis()->SetLabelSize(0.040);

 h_t_vs_recoil_theta->GetXaxis()->SetTitleOffset(1.15);
 h_t_vs_recoil_theta->GetYaxis()->SetTitleOffset(1.25);

 // Draw the histogram first
 h_t_vs_recoil_theta->Draw("COLZ");

 gPad->Update();

 double xMinT = h_t_vs_recoil_theta->GetXaxis()->GetXmin();
 double xMaxT = h_t_vs_recoil_theta->GetXaxis()->GetXmax();

 // Draw the assumed ALERT angular-selection boundaries
 TLine *line45_t = new TLine(xMinT, 45.0,xMaxT, 45.0);

 line45_t->SetLineColor(kRed + 1);
 line45_t->SetLineStyle(2);
 line45_t->SetLineWidth(3);
 line45_t->Draw("SAME");

 TLine *line90_t = new TLine(xMinT, 90.0, xMaxT, 90.0);

 line90_t->SetLineColor(kRed + 1);
 line90_t->SetLineStyle(2);
 line90_t->SetLineWidth(3);
 line90_t->Draw("SAME");

 c_t_theta->Modified();
 c_t_theta->Update();

 c_t_theta->SaveAs("minus_t_vs_recoil_He4_theta_before_cuts.png");
 
 
 
 TCanvas *c_tdiff_theta_after = new TCanvas("c_tdiff_theta_after", "tprime vs recoil angle after cuts",1000, 800);

 c_tdiff_theta_after->SetLeftMargin(0.12);
 c_tdiff_theta_after->SetRightMargin(0.15);
 c_tdiff_theta_after->SetBottomMargin(0.12);
 c_tdiff_theta_after->SetTopMargin(0.10);

 h_tdiff_vs_theta_after->SetStats(0);

 h_tdiff_vs_theta_after->GetXaxis()->CenterTitle();
 h_tdiff_vs_theta_after->GetYaxis()->CenterTitle();

 h_tdiff_vs_theta_after->GetXaxis()->SetTitleSize(0.045);
 h_tdiff_vs_theta_after->GetYaxis()->SetTitleSize(0.045);

 h_tdiff_vs_theta_after->GetXaxis()->SetLabelSize(0.040);
 h_tdiff_vs_theta_after->GetYaxis()->SetLabelSize(0.040);

 h_tdiff_vs_theta_after->GetXaxis()->SetTitleOffset(1.15);
 h_tdiff_vs_theta_after->GetYaxis()->SetTitleOffset(1.20);

 h_tdiff_vs_theta_after->Draw("COLZ");

 c_tdiff_theta_after->Modified();
 c_tdiff_theta_after->Update();

 c_tdiff_theta_after->SaveAs(
     "tprime_vs_recoil_He4_theta_after_cuts.png"
 );
 
 
 
 
 
 TCanvas *c_tmin_vs_th_after = new TCanvas( "c_tmin_vs_th_after", "-tmin,h vs -th after tprime cut",1000, 850);

 c_tmin_vs_th_after->SetLeftMargin(0.12);
 c_tmin_vs_th_after->SetRightMargin(0.15);
 c_tmin_vs_th_after->SetBottomMargin(0.12);
 c_tmin_vs_th_after->SetTopMargin(0.10);

 h_tmin_vs_th_after->SetStats(0);

 h_tmin_vs_th_after->GetXaxis()->CenterTitle();
 h_tmin_vs_th_after->GetYaxis()->CenterTitle();

 h_tmin_vs_th_after->GetXaxis()->SetTitleSize(0.045);
 h_tmin_vs_th_after->GetYaxis()->SetTitleSize(0.045);

 h_tmin_vs_th_after->GetXaxis()->SetLabelSize(0.040);
 h_tmin_vs_th_after->GetYaxis()->SetLabelSize(0.040);

 h_tmin_vs_th_after->GetXaxis()->SetTitleOffset(1.15);
 h_tmin_vs_th_after->GetYaxis()->SetTitleOffset(1.20);

 h_tmin_vs_th_after->Draw("COLZ");

 double axisMin1 = 0.0;
 double axisMax1 = 5.0;

 // t' = 0: -t_h = -t_min,h
 TLine *afterLine0 = new TLine(axisMin1, axisMin1,axisMax1, axisMax1);

 afterLine0->SetLineColor(kRed + 1);
 afterLine0->SetLineStyle(2);
 afterLine0->SetLineWidth(3);
 afterLine0->Draw("SAME");

 // t' = 1 GeV^2: -t_h = -t_min,h + 1
 TLine *afterLine1 = new TLine(axisMin1, axisMin1 + 1.0,axisMax1 - 1.0, axisMax1);

 afterLine1->SetLineColor(kMagenta + 1);
 afterLine1->SetLineStyle(2);
 afterLine1->SetLineWidth(3);
 afterLine1->Draw("SAME");

 c_tmin_vs_th_after->Modified();
 c_tmin_vs_th_after->Update();

 c_tmin_vs_th_after->SaveAs("minus_tmin_h_vs_minus_th_after_tprime_cut.png");
 
 
 
 //MC validation starts
 
 
 
 TCanvas *c_MKK = new TCanvas("c_MKK","Generated K+K- invariant mass",1000, 800);

 c_MKK->SetLeftMargin(0.13);
 c_MKK->SetRightMargin(0.05);
 c_MKK->SetBottomMargin(0.12);
 c_MKK->SetTopMargin(0.10);

 h_MKK->SetStats(0);
 h_MKK->SetLineColor(kBlue + 1);
 h_MKK->SetLineWidth(3);
 h_MKK->SetFillColorAlpha(kAzure + 1, 0.30);

 h_MKK->GetXaxis()->CenterTitle();
 h_MKK->GetYaxis()->CenterTitle();

 h_MKK->GetXaxis()->SetTitleSize(0.045);
 h_MKK->GetYaxis()->SetTitleSize(0.045);

 h_MKK->GetXaxis()->SetLabelSize(0.040);
 h_MKK->GetYaxis()->SetLabelSize(0.040);

 h_MKK->GetXaxis()->SetTitleOffset(1.15);
 h_MKK->GetYaxis()->SetTitleOffset(1.35);

 // Draw histogram first
 h_MKK->Draw("HIST");

 // Nominal phi-meson mass
 double nominalPhiMass = 1.019461;  // GeV/c^2
 double yMaximumMKK = h_MKK->GetMaximum();

 TLine *phiMassLine = new TLine(nominalPhiMass, 0.0,nominalPhiMass, 1.05*yMaximumMKK);

 phiMassLine->SetLineColor(kRed + 1);
 phiMassLine->SetLineStyle(2);
 phiMassLine->SetLineWidth(3);
 phiMassLine->Draw("SAME");

 // Optional label
 TLatex *phiMassLabel = new TLatex(nominalPhiMass + 0.002,0.85*yMaximumMKK,"m_{#phi} = 1.0195 GeV/c^{2}");

 phiMassLabel->SetTextColor(kRed + 1);
 phiMassLabel->SetTextSize(0.035);
 phiMassLabel->Draw("SAME");

 c_MKK->Modified();
 c_MKK->Update();

 c_MKK->SaveAs("generator_validation_M_Kplus_Kminus.png");
 
 
 TCanvas *c_deltaE = new TCanvas("c_deltaE","Energy conservation",1000, 800);

 c_deltaE->SetLeftMargin(0.13);
 c_deltaE->SetRightMargin(0.05);
 c_deltaE->SetBottomMargin(0.12);
 c_deltaE->SetTopMargin(0.10);

 h_deltaE->SetStats(0);
 h_deltaE->SetLineColor(kBlue + 1);
 h_deltaE->SetLineWidth(3);
 h_deltaE->SetFillColorAlpha(kAzure + 1, 0.30);

 h_deltaE->GetXaxis()->CenterTitle();
 h_deltaE->GetYaxis()->CenterTitle();

 h_deltaE->GetXaxis()->SetTitleSize(0.045);
 h_deltaE->GetYaxis()->SetTitleSize(0.045);

 h_deltaE->GetXaxis()->SetLabelSize(0.040);
 h_deltaE->GetYaxis()->SetLabelSize(0.040);

 h_deltaE->GetXaxis()->SetTitleOffset(1.15);
 h_deltaE->GetYaxis()->SetTitleOffset(1.35);

 h_deltaE->Draw("HIST");

 // Draw the expected Delta E = 0 position
 double yMaximumDeltaE = h_deltaE->GetMaximum();

 TLine *deltaEZeroLine = new TLine(0.0, 0.0,0.0, 1.05*yMaximumDeltaE);

 deltaEZeroLine->SetLineColor(kRed + 1);
 deltaEZeroLine->SetLineStyle(2);
 deltaEZeroLine->SetLineWidth(3);
 deltaEZeroLine->Draw("SAME");

 c_deltaE->Modified();
 c_deltaE->Update();

 c_deltaE->SaveAs("generator_validation_deltaE.png");
 
 
 TCanvas *c_deltaP = new TCanvas("c_deltaP", "Momentum conservation", 1000, 800);

 c_deltaP->SetLeftMargin(0.13);
 c_deltaP->SetRightMargin(0.05);
 c_deltaP->SetBottomMargin(0.12);
 c_deltaP->SetTopMargin(0.10);

 h_deltaP->SetStats(0);
 h_deltaP->SetLineColor(kBlue + 1);
 h_deltaP->SetLineWidth(3);
 h_deltaP->SetFillColorAlpha(kAzure + 1, 0.30);

 h_deltaP->GetXaxis()->CenterTitle();
 h_deltaP->GetYaxis()->CenterTitle();

 h_deltaP->GetXaxis()->SetTitleSize(0.045);
 h_deltaP->GetYaxis()->SetTitleSize(0.045);

 h_deltaP->GetXaxis()->SetLabelSize(0.040);
 h_deltaP->GetYaxis()->SetLabelSize(0.040);

 h_deltaP->GetXaxis()->SetTitleOffset(1.15);
 h_deltaP->GetYaxis()->SetTitleOffset(1.35);

 h_deltaP->Draw("HIST");

 c_deltaP->Modified();
 c_deltaP->Update();

 c_deltaP->SaveAs("generator_validation_deltaP_magnitude.png");
 
 
 
 TCanvas *c_delta_t_recoil_phi = new TCanvas("c_delta_t_recoil_phi", "Comparison of t definitions",1000, 800);

 c_delta_t_recoil_phi->SetLeftMargin(0.13);
 c_delta_t_recoil_phi->SetRightMargin(0.05);
 c_delta_t_recoil_phi->SetBottomMargin(0.12);
 c_delta_t_recoil_phi->SetTopMargin(0.10);

 h_delta_t_recoil_phi->SetStats(0);
 h_delta_t_recoil_phi->SetLineColor(kBlue + 1);
 h_delta_t_recoil_phi->SetLineWidth(3);
 h_delta_t_recoil_phi->SetFillColorAlpha(kAzure + 1, 0.30);

 h_delta_t_recoil_phi->GetXaxis()->CenterTitle();
 h_delta_t_recoil_phi->GetYaxis()->CenterTitle();

 h_delta_t_recoil_phi->GetXaxis()->SetTitleSize(0.045);
 h_delta_t_recoil_phi->GetYaxis()->SetTitleSize(0.045);

 h_delta_t_recoil_phi->GetXaxis()->SetLabelSize(0.040);
 h_delta_t_recoil_phi->GetYaxis()->SetLabelSize(0.040);

 h_delta_t_recoil_phi->GetXaxis()->SetTitleOffset(1.15);
 h_delta_t_recoil_phi->GetYaxis()->SetTitleOffset(1.35);

 h_delta_t_recoil_phi->Draw("HIST");

 // Expected position: t_recoil - t_phi = 0
 double yMaximumDeltaT = h_delta_t_recoil_phi->GetMaximum();

 TLine *deltaTZeroLine = new TLine(0.0, 0.0,0.0, 1.05*yMaximumDeltaT);

 deltaTZeroLine->SetLineColor(kRed + 1);
 deltaTZeroLine->SetLineStyle(2);
 deltaTZeroLine->SetLineWidth(3);
 deltaTZeroLine->Draw("SAME");

 c_delta_t_recoil_phi->Modified();
 c_delta_t_recoil_phi->Update();

 c_delta_t_recoil_phi->SaveAs("generator_validation_t_recoil_minus_t_phi.png");
 
 
 
 
 
 
 
 auto normalizeModelHistogram = [](TH1D *hist)
 {
     double normalization =
      hist->Integral("width");

     if(std::isfinite(normalization) &&
        normalization > 0.0) {

         hist->Scale(1.0/normalization);
     }
 };

 normalizeModelHistogram(hMinusTCompactModel);

 normalizeModelHistogram(hMinusTChargeModel);

 normalizeModelHistogram(hMinusTBroadModel);

 normalizeModelHistogram(hMinusTSmoothModel);
 
 
 
 
 // Compact
 hMinusTCompactModel->SetLineColor(kRed + 1);
 hMinusTCompactModel->SetLineWidth(4);
 hMinusTCompactModel->SetLineStyle(1);
 hMinusTCompactModel->SetMarkerColor(kRed + 1);

 // Charge-like
 hMinusTChargeModel->SetLineColor(kBlue + 1);
 hMinusTChargeModel->SetLineWidth(4);
 hMinusTChargeModel->SetLineStyle(1);
 hMinusTChargeModel->SetMarkerColor(kBlue + 1);

 // Broad
 hMinusTBroadModel->SetLineColor(kGreen + 2);
 hMinusTBroadModel->SetLineWidth(4);
 hMinusTBroadModel->SetLineStyle(1);
 hMinusTBroadModel->SetMarkerColor(kGreen + 2);

 // Smooth exponential
 hMinusTSmoothModel->SetLineColor(kBlack);
 hMinusTSmoothModel->SetLineWidth(4);
 hMinusTSmoothModel->SetLineStyle(2);
 hMinusTSmoothModel->SetMarkerColor(kBlack);
 
 
 // Find the largest normalized value among all four models
 double modelMaximum = std::max({
     hMinusTCompactModel->GetMaximum(),
     hMinusTChargeModel->GetMaximum(),
     hMinusTBroadModel->GetMaximum(),
     hMinusTSmoothModel->GetMaximum()
 });
 
 

 TCanvas *cNormalizedModels = new TCanvas("cNormalizedModels","Normalized helium model comparison",1100, 850);

 cNormalizedModels->SetLeftMargin(0.15);
 cNormalizedModels->SetRightMargin(0.05);
 cNormalizedModels->SetBottomMargin(0.13);
 cNormalizedModels->SetTopMargin(0.10);

 // Essential for seeing the minima
 cNormalizedModels->SetLogy();

 hMinusTCompactModel->SetStats(0);

 //hMinusTCompactModel->SetTitle("Generator-Level ^{4}He Model Comparison;" "-t [GeV^{2}];" "1/#sigma d#sigma/d(-t) [GeV^{-2}]");
 
 //hMinusTCompactModel->SetTitle("Generator-Level ^{4}He Model Comparison;""-t [GeV^{2}];""Normalized d#sigma/d(-t) [GeV^{-2}]");
 hMinusTCompactModel->SetTitle("Generator-Level ^{4}He Model Comparison;" "-t [GeV^{2}];" "Normalized d#sigma/d(-t) [GeV^{-2}]");

 //hMinusTCompactModel->SetMinimum(1.0e-8);
 //hMinusTCompactModel->SetMaximum(1.0e2);
 hMinusTCompactModel->SetMinimum(1.0e-4);
 
 hMinusTCompactModel->SetMaximum(10.0*modelMaximum);
 
 
 hMinusTCompactModel->GetXaxis()->CenterTitle();
 hMinusTCompactModel->GetYaxis()->CenterTitle();

 hMinusTCompactModel->GetXaxis()->SetTitleSize(0.045);
 hMinusTCompactModel->GetYaxis()->SetTitleSize(0.045);

 hMinusTCompactModel->GetXaxis()->SetLabelSize(0.040);
 hMinusTCompactModel->GetYaxis()->SetLabelSize(0.040);

 hMinusTCompactModel->GetXaxis()->SetTitleOffset(1.15);
 hMinusTCompactModel->GetYaxis()->SetTitleOffset(1.55);
 
 
 hMinusTCompactModel->GetXaxis()->SetRangeUser( 0.0, 0.8);
 

 hMinusTCompactModel->Draw("HIST");
 hMinusTChargeModel->Draw("HIST SAME");
 hMinusTBroadModel->Draw("HIST SAME");
 hMinusTSmoothModel->Draw("HIST SAME");

 // Legend
 TLegend *modelLegend = new TLegend( 0.48, 0.61,0.91, 0.88);

 modelLegend->SetHeader("#lambda = R_{g}/R_{C}","C");

 modelLegend->SetBorderSize(0);
 modelLegend->SetFillStyle(0);
 modelLegend->SetTextSize(0.032);

 modelLegend->AddEntry(hMinusTCompactModel,"Compact: #lambda = 0.93","l");

 modelLegend->AddEntry(hMinusTChargeModel,"Charge-like: #lambda = 1.00","l");

 modelLegend->AddEntry(hMinusTBroadModel,"Broad: #lambda = 1.10","l");

 modelLegend->AddEntry(hMinusTSmoothModel,"Smooth exponential: no minimum","l");

 modelLegend->Draw();

 // Final canvas update
 cNormalizedModels->Modified();
 cNormalizedModels->Update();

 // Save
 cNormalizedModels->SaveAs( "normalized_He4_dsigma_dt_model_comparison.png");
 
 
 
 
 
 TCanvas *cNormalizedChargeModels = new TCanvas("cNormalizedModels","Normalized helium model comparison",1100, 850);

 cNormalizedChargeModels->SetLeftMargin(0.15);
 cNormalizedChargeModels->SetRightMargin(0.05);
 cNormalizedChargeModels->SetBottomMargin(0.13);
 cNormalizedChargeModels->SetTopMargin(0.10);

 // Essential for seeing the minima
 cNormalizedChargeModels->SetLogy();

 hMinusTChargeModel->SetStats(0);

 //hMinusTCompactModel->SetTitle("Generator-Level ^{4}He Model Comparison;" "-t [GeV^{2}];" "1/#sigma d#sigma/d(-t) [GeV^{-2}]");
 
 //hMinusTCompactModel->SetTitle("Generator-Level ^{4}He Model Comparison;""-t [GeV^{2}];""Normalized d#sigma/d(-t) [GeV^{-2}]");
 //hMinusTChargeModel->SetTitle("Generator-Level ^{4}He Model Comparison;" "-t [GeV^{2}];" "Normalized d#sigma/d(-t) [GeV^{-2}]");
 hMinusTChargeModel->SetTitle("Generator-Level ^{4}He Model Comparison;" "-t [GeV^{2}];""d^{3}#sigma/(dQ^{2} dx_{B} dt) [nb/GeV^{4}]");

 //hMinusTCompactModel->SetMinimum(1.0e-8);
 //hMinusTCompactModel->SetMaximum(1.0e2);
 hMinusTChargeModel->SetMinimum(1.0e-4);
 
 hMinusTChargeModel->SetMaximum(10.0*modelMaximum);
 
 
 hMinusTChargeModel->GetXaxis()->CenterTitle();
 hMinusTChargeModel->GetYaxis()->CenterTitle();

 hMinusTChargeModel->GetXaxis()->SetTitleSize(0.045);
 hMinusTChargeModel->GetYaxis()->SetTitleSize(0.045);

 hMinusTChargeModel->GetXaxis()->SetLabelSize(0.040);
 hMinusTChargeModel->GetYaxis()->SetLabelSize(0.040);

 hMinusTChargeModel->GetXaxis()->SetTitleOffset(1.15);
 hMinusTChargeModel->GetYaxis()->SetTitleOffset(1.55);
 
 
 hMinusTChargeModel->GetXaxis()->SetRangeUser( 0.0, 0.8);
 

 //hMinusTCompactModel->Draw("HIST");
 hMinusTChargeModel->Draw("HIST");
 //hMinusTBroadModel->Draw("HIST SAME");
 //hMinusTSmoothModel->Draw("HIST SAME")
 
 
 


}