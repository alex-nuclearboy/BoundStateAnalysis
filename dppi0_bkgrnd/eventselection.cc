/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2019-06
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2020-05
***********************************************/

//Selection criteria for the analysis of the simulation results of the pd -> (3He-eta)_bound -> pdpi0 -> pd2g reaction

#include "Wasa.hh"
#include "CDataManager.hh"
#include "CHistoManager.hh"
#include "CParameterManager.hh"
#include "CLog.hh"
#include "CConst.hh"
#include "EmsEvent.hh"

#include "WHitBank.hh"
#include "WHitScint.hh"
#include "WVertex.hh"
#include "WVertexBank.hh"
#include "WCluster.hh"
#include "WClusterBank.hh"
#include "WClusterChamb.hh"
#include <WClusterFinder.hh>
#include "WTrack.hh"
#include "WTrackBank.hh"
#include "CDTracksSimple.hh"
#include "FDFTHTracks.hh"

#include "TString.h"
#include "TMath.h"
#include <TVector3.h>
#include <TLorentzVector.h>
#include <Riostream.h>
#include <iostream>
#include <fstream>

#include "eventselection.hh"

ClassImp(eventselection);

eventselection::eventselection() {}

eventselection::eventselection(const char * name):CAnalysisModule(name) {

/////////////////////////////////////////GRAPHICAL CUT//////////////////////////////////////////

    Double_t x[11];
    Double_t y[11];

    x[0] = 0.009269; y[0] = 0.00819116;
    x[1] = 0.009269; y[1] = 0.00591918;
    x[2] = 0.078264; y[2] = 0.00390835;
    x[3] = 0.184810; y[3] = 0.00244919;
    x[4] = 0.279872; y[4] = 0.00198307;
    x[5] = 0.423476; y[5] = 0.00198307;
    x[6] = 0.423476; y[6] = 0.00350519;
    x[7] = 0.349413; y[7] = 0.00372213;
    x[8] = 0.260732; y[8] = 0.00417181;
    x[9] = 0.176516; y[9] = 0.00519560;
    x[10] = 0.009269; y[10] = 0.00819116;

    //main cut
    cutg0 = new TCutG("CUTG0",11);
    cutg0->SetVarX("");
    cutg0->SetVarY("");
    cutg0->SetTitle("Graph");
    cutg0->SetFillColor(0);
    cutg0->SetLineColor(2);
    cutg0->SetLineWidth(2);
    cutg0->SetPoint(0,x[0],y[0]);
    cutg0->SetPoint(1,x[1],y[1]);
    cutg0->SetPoint(2,x[2],y[2]);
    cutg0->SetPoint(3,x[3],y[3]);
    cutg0->SetPoint(4,x[4],y[4]);
    cutg0->SetPoint(5,x[5],y[5]);
    cutg0->SetPoint(6,x[6],y[6]);
    cutg0->SetPoint(7,x[7],y[7]);
    cutg0->SetPoint(8,x[8],y[8]);
    cutg0->SetPoint(9,x[9],y[9]);
    cutg0->SetPoint(10,x[10],y[10]);

    //for systematics
    Double_t dX = 0.0017;
    Double_t dY = 0.000105;

    cutg1 = new TCutG("CUTG1",11);
    cutg1->SetVarX("");
    cutg1->SetVarY("");
    cutg1->SetTitle("Graph");
    cutg1->SetFillColor(0);
    cutg1->SetLineColor(5);
    cutg1->SetLineWidth(2);
    cutg1->SetPoint(0,x[0]-dX,y[0]+dY);
    cutg1->SetPoint(1,x[1]-dX,y[1]-dY);
    cutg1->SetPoint(2,x[2],y[2]-dY);
    cutg1->SetPoint(3,x[3],y[3]-dY);
    cutg1->SetPoint(4,x[4],y[4]-dY);
    cutg1->SetPoint(5,x[5]+dX,y[5]-dY);
    cutg1->SetPoint(6,x[6]+dX,y[6]+dY);
    cutg1->SetPoint(7,x[7],y[7]+dY);
    cutg1->SetPoint(8,x[8],y[8]+dY);
    cutg1->SetPoint(9,x[9],y[9]+dY);
    cutg1->SetPoint(10,x[10]-dX,y[10]+dY);

    cutg2 = new TCutG("CUTG2",11);
    cutg2->SetVarX("");
    cutg2->SetVarY("");
    cutg2->SetTitle("Graph");
    cutg2->SetFillColor(1);
    cutg2->SetLineColor(3);
    cutg2->SetLineWidth(1);
    cutg2->SetPoint(0,x[0]+dX,y[0]-dY);
    cutg2->SetPoint(1,x[1]+dX,y[1]+dY);
    cutg2->SetPoint(2,x[2],y[2]+dY);
    cutg2->SetPoint(3,x[3],y[3]+dY);
    cutg2->SetPoint(4,x[4],y[4]+dY);
    cutg2->SetPoint(5,x[5]-dX,y[5]+dY);
    cutg2->SetPoint(6,x[6]-dX,y[6]-dY);
    cutg2->SetPoint(7,x[7],y[7]-dY);
    cutg2->SetPoint(8,x[8],y[8]-dY);
    cutg2->SetPoint(9,x[9],y[9]-dY);
    cutg2->SetPoint(10,x[10]+dX,y[10]-dY);
/*
/////////////////////////////////TIME IN CYKLE TO BEAM MOMENTUM/////////////////////////////////

    BMfile = fopen("Time.2.PBeam.calibration.dat", "r");    //file with beam momentum
    gBeamMomentum = new TGraph();
    Int_t time_cycle;
    Float_t beam_mom;
    Int_t i = -1;
    while (!feof(BMfile)) {
        i++;
        fscanf(BMfile, "%d %f\n", &time_cycle, &beam_mom);
        beam_mom = beam_mom/1000.;
        gBeamMomentum->SetPoint(i, time_cycle, beam_mom);
    }
*/
////////////////////////////////////////////////////////////////////////////////////////////////

    fDetectorTable = dynamic_cast<CCardWDET*>(gParameterManager->GetParameterObject("CCardWDET","default"));
    //FD table
    kFTH1_old = fDetectorTable->GetDet(CConst::kFTH)->Get1stPlane();
    kFRH1_old = fDetectorTable->GetDet(CConst::kFRH)->Get1stPlane();
    //kFVH1_old = fDetectorTable->GetDet(CConst::kFVH)->Get1stPlane();
    kFWC1_old = fDetectorTable->GetDet(CConst::kFWC)->Get1stPlane();
    //printf("Plane Numbers FTH,FRH,FVH,FWC: %i,%i,%i,%i \n",kFTH1_old,kFRH1_old,kFVH1_old,kFWC1_old);
    printf("Plane Numbers FTH,FRH,FVH,FWC: %i,%i,%i \n",kFTH1_old,kFRH1_old,kFWC1_old);

    //CD table
    //Get PlaneNumbers of first and last Planes of PS and SE
    kPSfirst_old = fDetectorTable->GetDet(CConst::kPSB)->Get1stPlane();     //141
    kPSlast_old = fDetectorTable->GetDet(CConst::kPSF)->Get1stPlane();      //143
    kSEfirst_old = fDetectorTable->GetDet(CConst::kSEB)->Get1stPlane();     //151
    kSElast_old = fDetectorTable->GetDet(CConst::kSEF)->Get1stPlane() + fDetectorTable->GetDet(CConst::kSEF)->GetWasaPlanes()-1;  //174
    gScreen<<kPSfirst_old<<"\t"<<kPSlast_old<<"\t"<<kSEfirst_old<<"\t"<<kSElast_old<<CLog::endl;

    fCDTrackFinder = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));    //CD
    if(fCDTrackFinder!=0) fCDTrackBank = fCDTrackFinder->GetTrackBank();

    fFDTrackFinder = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));          //FD
    if(fFDTrackFinder!=0) fFDTrackBank = fFDTrackFinder->GetTrackBank();

    WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
    fMCTrackBank  = MCTrf->GetTrackBank();
    fMCVertexBank = MCTrf->GetVertexBank();

    fEventHeader = dynamic_cast<REventWmcHeader*>(gDataManager->GetDataObject("REventWmcHeader","EventHeader"));    //WMC Event header
    fHeader = dynamic_cast<REventHeader*>(gDataManager->GetDataObject("REventHeader","Header"));                    //DATA Event Header

    //change Edep to Ekin (WasaParameters)
    //fFDEdep2Ekin = dynamic_cast<FDEdep2Ekin*>(gParameterManager->GetParameterObject("FDEdep2Ekin","default"));    //"default" is for protons

    SetupSpectra(name);

}

////////////////////////////////////////////////////////////////////////////////////////////////

eventselection::~eventselection() {}

void eventselection::ProcessEvent() {    //01//

    if (fProcessed) return;
    fProcessed = kTRUE;

/////////////////////////////////////////ANALYSIS START/////////////////////////////////////////
/*
    //get event weight. For real data and Pluto, weight should be always one,
    //but for MC data, histograms have to be filled with according event weight
    Double_t ww=1;
    if (gWasa->IsAnalysisMode(Wasa::kMCRaw)||
            gWasa->IsAnalysisMode(Wasa::kMCReco)||
            gWasa->IsAnalysisMode(Wasa::kMC))
        ww=fEventHeader->GetWeight();
    cout<<"ww= "<<ww<<endl;
    //Int_t RunNumber = SorterOption::GetIntValue("RunNumber");
*/
    Double_t beamMom;
    Double_t Q;

    //beam momentum offset for DATA analysis
    //const Double_t pbeam_offset = 0.004;  //(from O. Rundel)
    //const Double_t pbeam_offset = 0.;

    Double_t E_calib = 1.;      //energy correction factor for WMC
    //Double_t E_calib = 1.526;   //energy correction factor for DATA

    ////PARTICLE MASSES////
    const Double_t m_target = 1.875613; //deuteron target mass  [GeV]
    const Double_t m_beam = 0.938272;   //proton beam mass      [GeV]
    const Double_t m_3He = 2.808950;    //3He mass              [GeV]
    const Double_t m_n = 0.939565;      //neutron mass          [GeV]
    const Double_t m_p = 0.938272;      //proton mass           [GeV]
    const Double_t m_d = 1.875613;      //deuteron mass         [GeV]
    const Double_t m_pi0 = 0.13497;     //neutral pion mass     [GeV]
    const Double_t m_pi = 0.13957;      //charged pion mass     [GeV]
    const Double_t m_eta = 0.547853;    //eta mass              [GeV]

    const Double_t s_thr = m_3He + m_eta;   //invariant mass on threshold [GeV]

///////////////////////////////GENERATED (TRUE) EVENTS FROM PLUTO///////////////////////////////

    if (gWasa->IsAnalysisMode(Wasa::kMCRaw)||
            gWasa->IsAnalysisMode(Wasa::kMCReco)||
            gWasa->IsAnalysisMode(Wasa::kMC)) {     //A01//

        TVector3 vec_d_MC;
        TVector3 vec_p_MC;
        TVector3 vec_g1_MC;
        TVector3 vec_g2_MC;
        TVector3 vec_pi0_MC;

        TLorentzVector P_d_MC;
        TLorentzVector P_p_MC;
        TLorentzVector P_g1_MC;
        TLorentzVector P_g2_MC;
        TLorentzVector P_pi0_MC;

        Double_t Ekin_d_lab_MC, Ekin_p_lab_MC;
        Double_t E_d_lab_MC, E_p_lab_MC;
        Double_t p_d_lab_MC, p_p_lab_MC;
        Double_t Theta_d_lab_MC, Theta_p_lab_MC;
        Double_t Phi_d_lab_MC, Phi_p_lab_MC;

        Double_t Ekin_g1_lab_MC, Ekin_g2_lab_MC;
        Double_t p_g1_lab_MC, p_g2_lab_MC;
        Double_t Theta_g1_lab_MC, Theta_g2_lab_MC;
        Double_t Phi_g1_lab_MC, Phi_g2_lab_MC;

        Double_t Ekin_pi0_lab_MC;
        Double_t E_pi0_lab_MC;
        Double_t p_pi0_lab_MC;
        Double_t Theta_pi0_lab_MC;
        Double_t Phi_pi0_lab_MC;

        Int_t PType;
        WParticle *part = 0;

        WVertexIter vIt(fMCVertexBank);
        Int_t NrVertex = 0;

        while (WVertex* vert = dynamic_cast<WVertex*>(vIt.Next())) {    //A02//

            NrVertex++;

            for (Int_t particleindex = 0; particleindex < vert->NumberOfParticles(); particleindex++) {     //A03//

                part = vert->GetParticle(particleindex);
                PType = part->GetType();

                //cout<<"NrVertex: "<<NrVertex<<endl;
                //cout<<"particleindex: "<<particleindex<<endl;
                //cout<<"PType: "<<PType<<endl;

                ////LAB////

                //deuteron//
                if ((NrVertex == 2) && (PType == 45)) { //d//

                    Ekin_d_lab_MC = part->GetEkin();
                    Theta_d_lab_MC = part->GetTheta();
                    Phi_d_lab_MC = part->GetPhi();

                    //cout<<"Ekin_d_lab_MC = "<<Ekin_d_lab_MC<<endl;
                    //cout<<"Theta_d_lab_MC = "<<Theta_d_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"Phi_d_lab_MC = "<<Phi_d_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"m_d = "<<part->GetMass()<<endl;

                    p_d_lab_MC = TMath::Sqrt(Ekin_d_lab_MC*(Ekin_d_lab_MC + 2*m_d));    //E^2 = m^2 + p^2, E = Ekin + m => p^2 = Ekin^2 + 2*m*Ekin
                    //E_d_lab_MC = TMath::Sqrt(p_d_lab_MC*p_d_lab_MC + m_d*m_d);    //total energy
                    E_d_lab_MC = TMath::Sqrt(TMath::Power(p_d_lab_MC,2) + TMath::Power(m_d,2));   //total energy

                    vec_d_MC.SetMagThetaPhi(p_d_lab_MC,Theta_d_lab_MC,Phi_d_lab_MC);
                    P_d_MC.SetVectM(vec_d_MC,m_d);

                    //histograms
                    hEkin_d_lab_MC->Fill(Ekin_d_lab_MC);
                    hp_d_lab_MC->Fill(p_d_lab_MC);
                    hE_d_lab_MC->Fill(E_d_lab_MC);
                    hTheta_d_lab_MC->Fill((Theta_d_lab_MC*TMath::RadToDeg()));
                    hPhi_d_lab_MC->Fill(Phi_d_lab_MC*TMath::RadToDeg());
                    hEkin_vs_Theta_d_lab_MC->Fill(Ekin_d_lab_MC,(Theta_d_lab_MC*TMath::RadToDeg()));

                }   //d//

                //proton//
                if ((NrVertex == 2) && (PType == 14)) { //p//

                    Ekin_p_lab_MC = part->GetEkin();
                    Theta_p_lab_MC = part->GetTheta();
                    Phi_p_lab_MC = part->GetPhi();

                    //cout<<"Ekin_p_lab_MC = "<<Ekin_p_lab_MC<<endl;
                    //cout<<"Theta_p_lab_MC = "<<Theta_p_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"Phi_p_lab_MC = "<<Phi_p_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"m_p = "<<part->GetMass()<<endl;

                    p_p_lab_MC = TMath::Sqrt(Ekin_p_lab_MC*(Ekin_p_lab_MC + 2*m_p));    //E^2 = m^2 + p^2, E = Ekin + m => p^2 = Ekin^2 + 2*m*Ekin
                    //E_p_lab_MC = TMath::Sqrt(p_p_lab_MC*p_p_lab_MC + m_p*m_p);    //total energy
                    E_p_lab_MC = TMath::Sqrt(TMath::Power(p_p_lab_MC,2) + TMath::Power(m_p,2));   //total energy

                    vec_p_MC.SetMagThetaPhi(p_p_lab_MC,Theta_p_lab_MC,Phi_p_lab_MC);
                    P_p_MC.SetVectM(vec_p_MC,m_p);

                    //histograms
                    hEkin_p_lab_MC->Fill(Ekin_p_lab_MC);
                    hp_p_lab_MC->Fill(p_p_lab_MC);
                    hE_p_lab_MC->Fill(E_p_lab_MC);
                    hTheta_p_lab_MC->Fill((Theta_p_lab_MC*TMath::RadToDeg()));
                    hPhi_p_lab_MC->Fill(Phi_p_lab_MC*TMath::RadToDeg());
                    hEkin_vs_Theta_p_lab_MC->Fill(Ekin_p_lab_MC,(Theta_p_lab_MC*TMath::RadToDeg()));

                }   //p//

                //gamma #1//
                if ((NrVertex == 3) && (PType == 1) && (particleindex == 1)) {  //g1//

                    Ekin_g1_lab_MC = part->GetEkin();
                    Theta_g1_lab_MC = part->GetTheta();
                    Phi_g1_lab_MC = part->GetPhi();

                    //cout<<"Ekin_g1_lab_MC = "<<Ekin_g1_lab_MC<<endl;
                    //cout<<"Theta_g1_lab_MC = "<<Theta_g1_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"Phi_g1_lab_MC = "<<Phi_g1_lab_MC*TMath::RadToDeg()<<endl;

                    p_g1_lab_MC = Ekin_g1_lab_MC;

                    vec_g1_MC.SetMagThetaPhi(p_g1_lab_MC,Theta_g1_lab_MC,Phi_g1_lab_MC);
                    P_g1_MC.SetVectM(vec_g1_MC,0.);

                    //histograms
                    //hEkin_g1_lab_MC->Fill(Ekin_g1_lab_MC);
                    hp_g1_lab_MC->Fill(p_g1_lab_MC);
                    hTheta_g1_lab_MC->Fill(Theta_g1_lab_MC*TMath::RadToDeg());
                    hPhi_g1_lab_MC->Fill(Phi_g1_lab_MC*TMath::RadToDeg());
                    hEkin_vs_Theta_g1_lab_MC->Fill(Ekin_g1_lab_MC,Theta_g1_lab_MC*TMath::RadToDeg());

                }   //g1//

                //gamma #2//
                if ((NrVertex == 3) && (PType == 1) && (particleindex == 2)) {  //g2//

                    Ekin_g2_lab_MC = part->GetEkin();
                    Theta_g2_lab_MC = part->GetTheta();
                    Phi_g2_lab_MC = part->GetPhi();

                    //cout<<"Ekin_g2_lab_MC = "<<Ekin_g2_lab_MC<<endl;
                    //cout<<"Theta_g2_lab_MC = "<<Theta_g2_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"Phi_g2_lab_MC = "<<Phi_g2_lab_MC*TMath::RadToDeg()<<endl;

                    p_g2_lab_MC = Ekin_g2_lab_MC;

                    vec_g2_MC.SetMagThetaPhi(p_g2_lab_MC,Theta_g2_lab_MC,Phi_g2_lab_MC);
                    P_g2_MC.SetVectM(vec_g2_MC,0.);

                    //histograms
                    //hEkin_g2_lab_MC->Fill(Ekin_g1_lab_MC);
                    hp_g2_lab_MC->Fill(p_g2_lab_MC);
                    hTheta_g2_lab_MC->Fill(Theta_g2_lab_MC*TMath::RadToDeg());
                    hPhi_g2_lab_MC->Fill(Phi_g2_lab_MC*TMath::RadToDeg());
                    hEkin_vs_Theta_g2_lab_MC->Fill(Ekin_g2_lab_MC,Theta_g2_lab_MC*TMath::RadToDeg());

                }   //g2//

            }   //A03//

        }   //A02//

        hEkin_vs_Theta_gammas_lab_MC->Fill(Ekin_g1_lab_MC,Theta_g1_lab_MC*TMath::RadToDeg());
        hEkin_vs_Theta_gammas_lab_MC->Fill(Ekin_g2_lab_MC,Theta_g2_lab_MC*TMath::RadToDeg());

        ////pion////
        vec_pi0_MC = (P_g1_MC + P_g2_MC).Vect();
        P_pi0_MC.SetVectM(vec_pi0_MC,m_pi0);

        p_pi0_lab_MC = vec_pi0_MC.Mag();
        //E_pi0_lab_MC = TMath::Sqrt(p_pi0_lab_MC*p_pi0_lab_MC + m_pi0*m_pi0);
        E_pi0_lab_MC = (P_g1_MC + P_g2_MC).E();
        Ekin_pi0_lab_MC = E_pi0_lab_MC - m_pi0;
        Theta_pi0_lab_MC = (P_g1_MC + P_g2_MC).Theta();
        Phi_pi0_lab_MC = (P_g1_MC + P_g2_MC).Phi();

        //histograms
        hEkin_pi0_lab_MC->Fill(Ekin_pi0_lab_MC);
        hp_pi0_lab_MC->Fill(p_pi0_lab_MC);
        hE_pi0_lab_MC->Fill(E_pi0_lab_MC);
        hTheta_pi0_lab_MC->Fill(Theta_pi0_lab_MC*TMath::RadToDeg());
        hPhi_pi0_lab_MC->Fill(Phi_pi0_lab_MC*TMath::RadToDeg());
        hEkin_vs_Theta_pi0_lab_MC->Fill(Ekin_pi0_lab_MC,Theta_pi0_lab_MC*TMath::RadToDeg());

        ////Angles////
        Double_t OpeningAngle_g1_g2_lab_MC = (vec_g1_MC.Angle(vec_g2_MC))*TMath::RadToDeg();   //angle between gammas [deg]
        Double_t OpeningAngle_pi0_p_lab_MC = (vec_pi0_MC.Angle(vec_p_MC))*TMath::RadToDeg();   //angle between pion & proton [deg]

        hOpeningAngle_g1_g2_lab_MC->Fill(OpeningAngle_g1_g2_lab_MC);
        hOpeningAngle_pi0_p_lab_MC->Fill(OpeningAngle_pi0_p_lab_MC);

        hTheta_g1_vs_Theta_g2_lab_MC->Fill(Theta_g1_lab_MC*TMath::RadToDeg(),Theta_g2_lab_MC*TMath::RadToDeg());
        hTheta_p_vs_Theta_d_lab_MC->Fill((Theta_d_lab_MC*TMath::RadToDeg()),Theta_p_lab_MC*TMath::RadToDeg());

        if( (Theta_d_lab_MC >= 0.052) && (Theta_d_lab_MC <= 0.314) && (Theta_p_lab_MC >= 0.349) && (Theta_p_lab_MC <= 2.950) ) {
            hTheta_p_vs_Theta_d_lab_cut_MC->Fill((Theta_d_lab_MC*TMath::RadToDeg()),Theta_p_lab_MC*TMath::RadToDeg());
        }

        //FOUR-VECTORS//

        TVector3 vec_beam_MC = vec_d_MC + vec_p_MC + vec_g1_MC + vec_g2_MC;

        //beam//
        beamMom = vec_beam_MC.Mag();

        hp_beam_MC->Fill(beamMom);

        TLorentzVector P_b_MC;      //4-vector of the beam
        P_b_MC.SetVectM(vec_beam_MC,m_beam);

        //target//
        TVector3 vec_target_MC;
        vec_target_MC.SetMagThetaPhi(0.,0.,0.);

        TLorentzVector P_t_MC;      //4-vector of the target
        P_t_MC.SetVectM(vec_target_MC,m_target);

        //total//
        TLorentzVector P_tot_MC = P_b_MC + P_t_MC;  //total 4-vector

        //Double_t s = TMath::Sqrt(m_beam*m_beam + m_target*m_target + 2*m_target*TMath::Sqrt(m_beam*m_beam + beamMom*beamMom));  //corresponds to Q from Q=-70MeV to Q=30MeV
        Double_t sqMass = TMath::Power(m_beam,2) + TMath::Power(m_target,2);
        Double_t beamEnergy = TMath::Sqrt( TMath::Power(beamMom,2) + TMath::Power(m_beam,2) );
        Double_t s = TMath::Sqrt(sqMass + 2*m_target*beamEnergy);   //corresponds to Q from Q=-70MeV to Q=30MeV

        //cout<<"beamMom = "<<beamMom<<endl;
        //cout<<"s = "<<s<<endl;
        //cout<<"s_thr = "<<s_thr<<endl;

        Q = 1000*(s - s_thr);   //excess energy

        //cout<<"Q= "<<Q<<endl;

        hGenerated_Q->Fill(Q);

        //Geometrical acceptance of WASA-at-COSY detector
        if( (Theta_d_lab_MC >= 0.052) && (Theta_d_lab_MC <= 0.314) ) {
            hEkin_vs_Theta_d_lab_acc_MC->Fill(Ekin_d_lab_MC,Theta_d_lab_MC*TMath::RadToDeg());
        }

        if( (Theta_p_lab_MC >= 0.052) && (Theta_p_lab_MC <= 0.314) ) {
            hEkin_vs_Theta_p_lab_acc_MC->Fill(Ekin_p_lab_MC,Theta_p_lab_MC*TMath::RadToDeg());
        }

        if( (Theta_p_lab_MC >= 0.349) && (Theta_p_lab_MC <= 2.950) ) {
            hEkin_vs_Theta_p_lab_acc_MC->Fill(Ekin_p_lab_MC,(Theta_p_lab_MC*TMath::RadToDeg()));
        }

        if ( (Theta_g1_lab_MC >= 0.349) && (Theta_g1_lab_MC <= 2.950) ) {
            hEkin_vs_Theta_g1_lab_acc_MC->Fill(Ekin_g1_lab_MC,Theta_g1_lab_MC*TMath::RadToDeg());
        }

        if ( (Theta_g2_lab_MC >= 0.349) && (Theta_g2_lab_MC <= 2.950) ) {
            hEkin_vs_Theta_g2_lab_acc_MC->Fill(Ekin_g2_lab_MC,Theta_g2_lab_MC*TMath::RadToDeg());
        }

        if ( (Theta_g1_lab_MC >= 0.349) && (Theta_g1_lab_MC <= 2.950) && (Theta_g2_lab_MC >= 0.349) && (Theta_g2_lab_MC <= 2.950) ) {
            hEkin_vs_Theta_gammas_lab_acc_MC->Fill(Ekin_g1_lab_MC,Theta_g1_lab_MC*TMath::RadToDeg());
            hEkin_vs_Theta_gammas_lab_acc_MC->Fill(Ekin_g2_lab_MC,Theta_g2_lab_MC*TMath::RadToDeg());
        }

        if( (Theta_d_lab_MC >= 0.052) && (Theta_d_lab_MC <= 0.314) ) {
            if ( (Theta_g1_lab_MC >= 0.349) && (Theta_g1_lab_MC <= 2.950) && (Theta_g2_lab_MC >= 0.349) && (Theta_g2_lab_MC <= 2.950) ) {
                if ( (Theta_p_lab_MC >= 0.052) && (Theta_p_lab_MC <= 0.314) ) {
                    hAccepted_Q->Fill(Q);
                }
                else {
                    if ( (Theta_p_lab_MC >= 0.349) && (Theta_p_lab_MC <= 2.950) ) {
                        hAccepted_Q->Fill(Q);
                    }
                }
            }
        }

////////////////////////////////////////////CM FRAME////////////////////////////////////////////

        ////boost to CM////
        TVector3 b_MC;
        b_MC = P_tot_MC.BoostVector();  //boost to LAB

        ////beam and target////
        P_b_MC.Boost(-b_MC);    //to CM
        P_t_MC.Boost(-b_MC);    //to CM

        ////deuteron////
        P_d_MC.Boost(-b_MC);    //to CM

        TVector3 vec_d_cm_MC = P_d_MC.Vect();

        Double_t p_d_cm_MC = vec_d_cm_MC.Mag();
        //Double_t E_d_cm_MC = TMath::Sqrt(p_d_cm_MC*p_d_cm_MC + m_d*m_d);
        Double_t E_d_cm_MC = P_d_MC.E();
        Double_t Ekin_d_cm_MC = E_d_cm_MC - m_d;;
        Double_t Theta_d_cm_MC = P_d_MC.Theta();
        Double_t Phi_d_cm_MC = P_d_MC.Phi();

        //histograms
        hp_d_cm_MC->Fill(p_d_cm_MC);
        hE_d_cm_MC->Fill(E_d_cm_MC);
        hEkin_d_cm_MC->Fill(Ekin_d_cm_MC);
        hTheta_d_cm_MC->Fill(Theta_d_cm_MC*TMath::RadToDeg());
        hPhi_d_cm_MC->Fill(Phi_d_cm_MC*TMath::RadToDeg());
        hEkin_vs_Theta_d_cm_MC->Fill(Ekin_d_cm_MC,Theta_d_cm_MC*TMath::RadToDeg());

        ////proton////
        P_p_MC.Boost(-b_MC);    //to CM

        TVector3 vec_p_cm_MC = P_p_MC.Vect();
        Double_t p_p_cm_MC = vec_p_cm_MC.Mag();
        //Double_t E_p_cm_MC = TMath::Sqrt(p_p_cm_MC*p_p_cm_MC + m_p*m_p);
        Double_t E_p_cm_MC = P_p_MC.E();
        Double_t Ekin_p_cm_MC = E_p_cm_MC - m_p;
        Double_t Theta_p_cm_MC = P_p_MC.Theta();
        Double_t Phi_p_cm_MC = P_p_MC.Phi();

        //histograms
        hp_p_cm_MC->Fill(p_p_cm_MC);
        hE_p_cm_MC->Fill(E_p_cm_MC);
        hEkin_p_cm_MC->Fill(Ekin_p_cm_MC);
        hTheta_p_cm_MC->Fill(Theta_p_cm_MC*TMath::RadToDeg());
        hPhi_p_cm_MC->Fill(Phi_p_cm_MC*TMath::RadToDeg());
        hEkin_vs_Theta_p_cm_MC->Fill(Ekin_p_cm_MC,Theta_p_cm_MC*TMath::RadToDeg());

        //
        hTheta_p_vs_Theta_d_cm_MC->Fill(Theta_d_cm_MC*TMath::RadToDeg(),Theta_p_cm_MC*TMath::RadToDeg());

        ////gamma 1////
        P_g1_MC.Boost(-b_MC);   //to CM
        TVector3 vec_g1_cm_MC = P_g1_MC.Vect();
        Double_t p_g1_cm_MC = vec_g1_cm_MC.Mag();
        Double_t E_g1_cm_MC = P_g1_MC.E();
        Double_t Theta_g1_cm_MC = P_g1_MC.Theta();
        Double_t Phi_g1_cm_MC = P_g1_MC.Phi();

        //histograms
        hp_g1_cm_MC->Fill(p_g1_cm_MC);
        //hE_g1_cm_MC->Fill(E_g1_cm_MC);
        hTheta_g1_cm_MC->Fill(Theta_g1_cm_MC*TMath::RadToDeg());
        hPhi_g1_cm_MC->Fill(Phi_g1_cm_MC*TMath::RadToDeg());
        hEkin_vs_Theta_g1_cm_MC->Fill(p_g1_cm_MC,Theta_g1_cm_MC*TMath::RadToDeg());

        ////gamma 2////
        P_g2_MC.Boost(-b_MC);   //to CM
        TVector3 vec_g2_cm_MC = P_g2_MC.Vect();
        Double_t p_g2_cm_MC = vec_g2_cm_MC.Mag();
        Double_t E_g2_cm_MC = P_g2_MC.E();
        Double_t Theta_g2_cm_MC = P_g2_MC.Theta();
        Double_t Phi_g2_cm_MC = P_g2_MC.Phi();

        Double_t OpeningAngle_g1_g2_cm_MC = (vec_g1_cm_MC.Angle(vec_g2_cm_MC))*TMath::RadToDeg();   //angle between gammas in CM [deg]

        //histograms
        hp_g2_cm_MC->Fill(p_g2_cm_MC);
        //hE_g2_cm_MC->Fill(E_g2_cm_MC);
        hTheta_g2_cm_MC->Fill(Theta_g2_cm_MC*TMath::RadToDeg());
        hPhi_g2_cm_MC->Fill(Phi_g2_cm_MC*TMath::RadToDeg());
        hEkin_vs_Theta_g2_cm_MC->Fill(p_g2_cm_MC,Theta_g2_cm_MC*TMath::RadToDeg());

        hOpeningAngle_g1_g2_cm_MC->Fill(OpeningAngle_g1_g2_cm_MC);
        hTheta_g1_vs_Theta_g2_cm_MC->Fill(Theta_g1_cm_MC*TMath::RadToDeg(),Theta_g2_cm_MC*TMath::RadToDeg());

        ////pion////
        TVector3 vec_pi0_cm_MC = (P_g1_MC + P_g2_MC).Vect();
        Double_t p_pi0_cm_MC = vec_pi0_cm_MC.Mag();
        //Double_t E_pi0_cm_MC = TMath::Sqrt(p_pi0_cm_MC*p_pi0_cm_MC+m_pi0*m_pi0);
        Double_t E_pi0_cm_MC = (P_g1_MC + P_g2_MC).E();
        Double_t Ekin_pi0_cm_MC = E_pi0_cm_MC - m_pi0;
        Double_t Theta_pi0_cm_MC = (P_g1_MC + P_g2_MC).Theta();
        Double_t Phi_pi0_cm_MC = (P_g1_MC + P_g2_MC).Phi();

        //histograms
        hp_pi0_cm_MC->Fill(p_pi0_cm_MC);
        hE_pi0_cm_MC->Fill(E_pi0_cm_MC);
        hEkin_pi0_cm_MC->Fill(Ekin_pi0_cm_MC);
        hTheta_pi0_cm_MC->Fill(Theta_pi0_cm_MC*TMath::RadToDeg());
        hPhi_pi0_cm_MC->Fill(Phi_pi0_cm_MC*TMath::RadToDeg());
        hEkin_vs_Theta_pi0_cm_MC->Fill(Ekin_pi0_cm_MC,Theta_pi0_cm_MC*TMath::RadToDeg());

        //Invariant and Missing Masses in CM//
        Double_t InvariantMass_MC = (P_g1_MC + P_g2_MC).M();
        Double_t MissingMass_d_MC = (P_b_MC + P_t_MC - P_p_MC - P_g1_MC - P_g2_MC).M();
        Double_t MissingMass_p_MC = (P_b_MC + P_t_MC - P_d_MC - P_g1_MC - P_g2_MC).M();
        Double_t MissingEnergy_d_MC = (P_b_MC + P_t_MC - P_p_MC - P_g1_MC - P_g2_MC).E();
        Double_t MissingEnergy_p_MC = (P_b_MC + P_t_MC - P_d_MC - P_g1_MC - P_g2_MC).E();
        Double_t TotalMass_MC = (P_d_MC + P_p_MC + P_g1_MC + P_g2_MC).M();

        Double_t OpeningAngle_pi0_p_cm_MC = (vec_pi0_cm_MC.Angle(vec_p_cm_MC))*TMath::RadToDeg();

        //histograms
        hIM_pion_MC->Fill(InvariantMass_MC);
        hMM_deuteron_MC->Fill(MissingMass_d_MC);
        hMM_proton_MC->Fill(MissingMass_p_MC);
        hME_deuteron_MC->Fill(MissingEnergy_d_MC);
        hME_proton_MC->Fill(MissingEnergy_p_MC);
        hMMvsME_nucleon_MC->Fill(MissingMass_d_MC,MissingEnergy_d_MC);

        hOpeningAngle_pi0_p_cm_MC->Fill(OpeningAngle_pi0_p_cm_MC);

    }   //A01//

//////////////////////////////////////RECONSTRUCTED EVENTS//////////////////////////////////////

    ///////TRIGGER///////

    hStatistics[0]->Fill(0);

    //if (!(fHeader->TriggerNumSet(10)))  return; //trigger #10

    hStatistics[0]->Fill(1);

    //
    Int_t NumNeutTrackCD = fCDTrackBank->GetEntries(11);    //Neutral Tracks in CD have Type 11
    Int_t NumCharTrackCD = fCDTrackBank->GetEntries(12);    //Charged Tracks in CD have Type 12
    Int_t NumCharTrackFD = fFDTrackBank->GetEntries(2);     //Charged Tracks in FD have type 2

    //////LEVEL 0//////

    hNeutralTracksCD[0][0]->Fill(NumNeutTrackCD);
    hChargedTracksCD[0][0]->Fill(NumCharTrackCD);
    hChargedTracksFD[0][0]->Fill(NumCharTrackFD);
/*
    //BEAM MOMENTUM//
    Double_t t_incycle = 1000*fHeader->GetTimeInCycle();        //time in cycle (t_axis - t_start)
    beamMom = (gBeamMomentum->Eval(t_incycle)) + pbeam_offset;  //beam momentum [GeV/c]

    Double_t sqMass = TMath::Power(m_beam, 2) + TMath::Power(m_target, 2);
    Double_t beamEnergy = TMath::Sqrt(TMath::Power(beamMom, 2) + TMath::Power(m_beam, 2));
    Double_t s = TMath::Sqrt(sqMass + 2*m_target*beamEnergy);  //corresponds to Q from Q=-70MeV to Q=30MeV
    //Double_t s = TMath::Sqrt(m_beam*m_beam + m_target*m_target +2*m_target*TMath::Sqrt(m_beam*m_beam + beamMom*beamMom)); //corresponds to Q from Q=-70MeV to Q=30MeV

    Q = 1000*(s - s_thr);   //excess energy [MeV]
*/
    hp_beam[0][0]->Fill(beamMom);
    hQ[0][0][0][0]->Fill(Q);

    //FOUR-VECTORS//
    //beam and target//
    TVector3 vec_beam;
    vec_beam.SetMagThetaPhi(beamMom,0.,0.);
    TLorentzVector P_b;                 //4-vector of the beam
    P_b.SetVectM(vec_beam,m_beam);

    TVector3 vec_target;
    vec_target.SetMagThetaPhi(0.,0.,0.);
    TLorentzVector P_t;                 //4-vector of the target
    P_t.SetVectM(vec_target,m_target);

    //total 4-vector//
    TLorentzVector P_tot = P_b + P_t;

////////////////////////////////////////////////////////////////////////////////////////////////

    //Forward Detector//
    WTrackIter TrackIterFD(fFDTrackBank);   //define an Iterator for th FD TrackBank
    TrackIterFD.SetType(2);                 //prepare Iterator to look at charged tracks only

    while (WTrack *FDTrack = dynamic_cast<WTrack*> (TrackIterFD.Next())) {

        if(FDTrack->Theta() != 0.125) {

            hTime_FDC[0][0]->Fill(FDTrack->Time());
            hTheta_FDC[0][0]->Fill(FDTrack->Theta()*TMath::RadToDeg());
            hPhi_FDC[0][0]->Fill(FDTrack->Phi()*TMath::RadToDeg());

            //energy deposited in different FD layers
            hEdepFWC1vsFRH1[0][0]->Fill(FDTrack->Edep(kFWC1),FDTrack->Edep(kFRH1));
            hEdepFWC2vsFRH1[0][0]->Fill(FDTrack->Edep(kFWC2),FDTrack->Edep(kFRH1));
            hEdepFTH1vsFRH1[0][0]->Fill(FDTrack->Edep(kFTH1),FDTrack->Edep(kFRH1));
            hEdepFRH1vsFRH2[0][0]->Fill(FDTrack->Edep(kFRH1),FDTrack->Edep(kFRH2));
            hEdepFRH2vsFRH3[0][0]->Fill(FDTrack->Edep(kFRH2),FDTrack->Edep(kFRH3));
            hEdepFWC1vsFRH1FRH2FRH3[0][0]->Fill(FDTrack->Edep(kFWC1),(FDTrack->Edep(kFRH1)+FDTrack->Edep(kFRH2)+FDTrack->Edep(kFRH3)));

        }

    }

    //Central Detector//
    WTrackIter TrackIterCD(fCDTrackBank);   //define an Iterator for th CD TrackBank
    TrackIterCD.SetType(12);                //prepare Iterator to look at charged tracks only

    while (WTrack *CDTrack = dynamic_cast<WTrack*> (TrackIterCD.Next())) {

        hTime_CDC[0][0]->Fill(CDTrack->Time());
        hTheta_CDC[0][0]->Fill(CDTrack->Theta()*TMath::RadToDeg());
        hPhi_CDC[0][0]->Fill(CDTrack->Phi()*TMath::RadToDeg());
        hMom_CDC[0][0]->Fill(CDTrack->Momentum());

        hEdepPSBvsSEC[0][0]->Fill((CDTrack->Edep(151,174)),(CDTrack->GetSpecELossPS()*E_calib));
        hEdepPSBvsSigMom[0][0]->Fill((CDTrack->Momentum()*CDTrack->Charge()),(CDTrack->GetSpecELossPS()*E_calib));
        hEdepSECvsSigMom[0][0]->Fill((CDTrack->Momentum()*CDTrack->Charge()),(CDTrack->Edep(151,174)));

    }

    TrackIterCD.Reset();        //reset Iterator and
    TrackIterCD.SetType(11);    //prepare for looking at neutral tracks

    while (WTrack *CDTrack = dynamic_cast<WTrack*> (TrackIterCD.Next())) {

        hTime_CDN[0][0]->Fill(CDTrack->Time());
        hTheta_CDN[0][0]->Fill(CDTrack->Theta()*TMath::RadToDeg());
        hPhi_CDN[0][0]->Fill(CDTrack->Phi()*TMath::RadToDeg());
        hMom_CDN[0][0]->Fill(CDTrack->Momentum());
    }

////////////////////////////////////////////////////////////////////////////////////////////////

    //////LEVEL 1//////
    //exactly 1 charged particle in FD & exactly 1 charged particle in CD & at least 2 netral clusters in CD//

    if((NumCharTrackCD == 1) && (NumCharTrackFD == 1) && (NumNeutTrackCD >=2)) {    //B01//

        hStatistics[0]->Fill(2);
        hStatistics[1]->Fill(0);
        hStatistics[2]->Fill(0);

        hNeutralTracksCD[1][0]->Fill(NumNeutTrackCD);
        hChargedTracksCD[1][0]->Fill(NumCharTrackCD);
        hChargedTracksFD[1][0]->Fill(NumCharTrackFD);

        //Forward Detector//
        Double_t EdepFWC1, EdepFWC2;
        Double_t EdepFTH1;
        Double_t EdepFRH1, EdepFRH2, EdepFRH3;
        Double_t TimeFD;
        Double_t ThetaFD_lab;
        Double_t PhiFD_lab;

        WTrackIter FDTrackIter(fFDTrackBank);   //define an iterator for FD TrackBank
        FDTrackIter.SetType(2);                 //prepare Iterator to look at charged tracks only

        while (WTrack *TrackFD = dynamic_cast<WTrack*> (FDTrackIter.Next())) {

            TimeFD = TrackFD->Time();
            ThetaFD_lab = TrackFD->Theta();
            PhiFD_lab = TrackFD->Phi();

            EdepFWC1 = TrackFD->Edep(kFWC1);
            EdepFWC2 = TrackFD->Edep(kFWC2);
            EdepFTH1 = TrackFD->Edep(kFTH1);
            EdepFRH1 = TrackFD->Edep(kFRH1);
            EdepFRH2 = TrackFD->Edep(kFRH2);
            EdepFRH3 = TrackFD->Edep(kFRH3);

        }

        //Central Detector//
        Double_t EdepSEC;
        Double_t EdepPSB;
        Double_t SgnMom;

        Double_t TimeCD;
        Double_t ThetaCD_lab;
        Double_t PhiCD_lab;
        Double_t MomCD_lab;

        WTrackIter CDTrackIter(fCDTrackBank);   //define an Iterator for CD TrackBank
        CDTrackIter.SetType(12);                //prepare Iterator to look at charged tracks only

        while (WTrack *TrackCD = dynamic_cast<WTrack*> (CDTrackIter.Next())) {

            TimeCD = TrackCD->Time();
            ThetaCD_lab = TrackCD->Theta();
            PhiCD_lab = TrackCD->Phi();
            MomCD_lab = TrackCD->Momentum();

            EdepSEC = TrackCD->Edep(151,174);
            EdepPSB = TrackCD->GetSpecELossPS()*E_calib;
            SgnMom = TrackCD->Momentum()*(TrackCD->Charge());

        }

        //Loop over tracks in CD to look the best gamma pair

        TLorentzVector Gamma1;
        TLorentzVector Gamma2;
        WTrack *TrackGamma1;
        WTrack *TrackGamma2;

        Double_t best_delta = 99999;

        for (Int_t l = 0; l < fCDTrackBank->GetEntries(); l++) {

            WTrack *track1 = fCDTrackBank->GetTrack(l);

            if (track1->Type() == kCDN) {

                for (Int_t k = l + 1; k < fCDTrackBank->GetEntries(); k++) {

                    WTrack *track2 = fCDTrackBank->GetTrack(k);

                    if (track2->Type() == kCDN) {

                        Double_t Mom1 = track1->Momentum();
                        Double_t Mom2 = track2->Momentum();
                        Double_t Theta1 = track1->Theta();
                        Double_t Theta2 = track2->Theta();
                        Double_t Phi1 = track1->Phi();
                        Double_t Phi2 = track2->Phi();

                        TVector3 vec_1;
                        vec_1.SetMagThetaPhi(Mom1,Theta1,Phi1);
                        TLorentzVector P_1;
                        P_1.SetVectM(vec_1,0.);

                        TVector3 vec_2;
                        vec_2.SetMagThetaPhi(Mom2,Theta2,Phi2);
                        TLorentzVector P_2;
                        P_2.SetVectM(vec_2,0.);

                        Double_t inv_m = (P_1 + P_2).M();
                        Double_t delta = TMath::Abs(inv_m - m_pi0);

                        if (delta < best_delta) {

                            best_delta = delta;

                            Gamma1 = (Mom1 >= Mom2)?P_1:P_2;    //higher energy gamma has number one
                            Gamma2 = (Mom1 >= Mom2)?P_2:P_1;    //lower energy gamma has number two
                            TrackGamma1 = (Mom1 >= Mom2)?track1:track2;
                            TrackGamma2 = (Mom1 >= Mom2)?track2:track1;

                        }
                    }
                }
            }
        }

        //gamma quanta//
        Double_t p_g_lab[2], Theta_g_lab[2], Phi_g_lab[2];

        p_g_lab[0] = TrackGamma1->Momentum();
        Theta_g_lab[0] = TrackGamma1->Theta();
        Phi_g_lab[0] = TrackGamma1->Phi();

        TVector3 vec_g1_lab;
        vec_g1_lab.SetMagThetaPhi(p_g_lab[0],Theta_g_lab[0],Phi_g_lab[0]);
        TLorentzVector P_g1;
        P_g1.SetVectM(vec_g1_lab,0.);

        p_g_lab[1] = TrackGamma2->Momentum();
        Theta_g_lab[1] = TrackGamma2->Theta();
        Phi_g_lab[1] = TrackGamma2->Phi();

        TVector3 vec_g2_lab;
        vec_g2_lab.SetMagThetaPhi(p_g_lab[1],Theta_g_lab[1],Phi_g_lab[1]);
        TLorentzVector P_g2;
        P_g2.SetVectM(vec_g2_lab,0.);

        Double_t OpeningAngle_g1_g2_lab = (vec_g1_lab.Angle(vec_g2_lab))*TMath::RadToDeg(); //angle between gamma quanta in LAB

        //pion//
        TVector3 vec_pi0_lab;
        vec_pi0_lab = (P_g1 + P_g2).Vect();
        TLorentzVector P_pi0;
        P_pi0.SetVectM(vec_pi0_lab,m_pi0);

        Double_t p_pi0_lab = vec_pi0_lab.Mag();
        Double_t E_pi0_lab = TMath::Sqrt(TMath::Power(p_pi0_lab,2) + TMath::Power(m_pi0,2));
        //Double_t E_pi0_lab = (P_g1 + P_g2).E();
        Double_t Ekin_pi0_lab = E_pi0_lab - m_pi0;
        Double_t Theta_pi0_lab = (P_g1 + P_g2).Theta();
        Double_t Phi_pi0_lab = (P_g1 + P_g2).Phi();

///////////////////////////////////ADDITIONAL GAMMA QUANTA CUT//////////////////////////////////

        for (Int_t j = 0; j < fCDTrackBank->GetEntries(); j++) {
            WTrack *track = fCDTrackBank->GetTrack(j);

            if (track->Type() == kCDN) {

                Double_t p_g = track->Momentum();

                if(TMath::Abs(p_g-p_g_lab[0])>1e-9 && TMath::Abs(p_g-p_g_lab[1])>1e-9) {

                    hEnergy_additional_gammas->Fill(p_g);

                    if(p_g>0.03) return;

                    hEnergy_additional_gammas_cut->Fill(p_g);

                }
            }
        }

////////////////////////////////////////////////////////////////////////////////////////////////

        if(ThetaFD_lab == 0.125) return;

        hStatistics[1]->Fill(1);
        hStatistics[2]->Fill(1);

        hQ[1][0][0][0]->Fill(Q);

        hTime_FDC[1][0]->Fill(TimeFD);
        hTheta_FDC[1][0]->Fill(ThetaFD_lab*TMath::RadToDeg());
        hPhi_FDC[1][0]->Fill(PhiFD_lab*TMath::RadToDeg());

        hEdepFWC1vsFRH1[1][0]->Fill(EdepFWC1,EdepFRH1);
        hEdepFWC2vsFRH1[1][0]->Fill(EdepFWC2,EdepFRH1);
        hEdepFTH1vsFRH1[1][0]->Fill(EdepFTH1,EdepFRH1);
        hEdepFRH1vsFRH2[1][0]->Fill(EdepFRH1,EdepFRH2);
        hEdepFRH2vsFRH3[1][0]->Fill(EdepFRH2,EdepFRH3);
        hEdepFWC1vsFRH1FRH2FRH3[1][0]->Fill(EdepFWC1,(EdepFRH1 + EdepFRH2 + EdepFRH3));

        hTime_CDC[1][0]->Fill(TimeCD);
        hMom_CDC[1][0]->Fill(MomCD_lab);
        hTheta_CDC[1][0]->Fill(ThetaCD_lab*TMath::RadToDeg());
        hPhi_CDC[1][0]->Fill(PhiCD_lab*TMath::RadToDeg());

        hEdepPSBvsSEC[1][0]->Fill(EdepSEC,EdepPSB);
        hEdepPSBvsSigMom[1][0]->Fill(SgnMom,EdepPSB);
        hEdepSECvsSigMom[1][0]->Fill(SgnMom,EdepSEC);

        hp_g1_lab[1][0]->Fill(p_g_lab[0]);
        hTheta_g1_lab[1][0]->Fill(Theta_g_lab[0]*TMath::RadToDeg());
        hPhi_g1_lab[1][0]->Fill(Phi_g_lab[0]*TMath::RadToDeg());
        hEkin_vs_Theta_g1_lab[1][0]->Fill(p_g_lab[0],Theta_g_lab[0]*TMath::RadToDeg());

        hp_g2_lab[1][0]->Fill(p_g_lab[1]);
        hTheta_g2_lab[1][0]->Fill(Theta_g_lab[1]*TMath::RadToDeg());
        hPhi_g2_lab[1][0]->Fill(Phi_g_lab[1]*TMath::RadToDeg());
        hEkin_vs_Theta_g2_lab[1][0]->Fill(p_g_lab[1],Theta_g_lab[1]*TMath::RadToDeg());

        hOpeningAngle_g1_g2_lab[1][0]->Fill(OpeningAngle_g1_g2_lab);
        hTheta_g1_vs_Theta_g2_lab[1][0]->Fill(Theta_g_lab[0]*TMath::RadToDeg(),Theta_g_lab[1]*TMath::RadToDeg());

        hp_pi0_lab[1][0]->Fill(p_pi0_lab);
        hE_pi0_lab[1][0]->Fill(E_pi0_lab);
        hEkin_pi0_lab[1][0]->Fill(Ekin_pi0_lab);
        hTheta_pi0_lab[1][0]->Fill(Theta_pi0_lab*TMath::RadToDeg());
        hPhi_pi0_lab[1][0]->Fill(Phi_pi0_lab*TMath::RadToDeg());
        hEkin_vs_Theta_pi0_lab[1][0]->Fill(Ekin_pi0_lab,Theta_pi0_lab*TMath::RadToDeg());

////////////////////////////////////////////////////////////////////////////////////////////////

        //proton//
        Double_t p_p_lab = MomCD_lab;
        Double_t Theta_p_lab = ThetaCD_lab;
        Double_t Phi_p_lab = PhiCD_lab;
        Double_t E_p_lab = TMath::Sqrt(TMath::Power(p_p_lab,2) + TMath::Power(m_p,2));
        Double_t Ekin_p_lab = E_p_lab - m_p;

        TVector3 vec_p_lab;
        vec_p_lab.SetMagThetaPhi(p_p_lab,Theta_p_lab,Phi_p_lab);
        TLorentzVector P_p;
        P_p.SetVectM(vec_p_lab,m_p);

        Double_t OpeningAngle_pi0_p_lab = (vec_pi0_lab.Angle(vec_p_lab))*TMath::RadToDeg(); //angle between pion & proton in LAB

        //deuteron (reconstructed)//
        TVector3 vec_d;
        vec_d = (P_b + P_t - P_p - P_g1 - P_g2).Vect();
        TLorentzVector P_d;
        P_d.SetVectM(vec_d,m_d);

        Double_t p_d_lab = vec_d.Mag();
        Double_t E_d_lab = TMath::Sqrt(TMath::Power(p_d_lab,2) + TMath::Power(m_d,2));
        Double_t Ekin_d_lab = E_d_lab - m_d;
        Double_t Theta_d_lab = P_d.Theta();
        Double_t Phi_d_lab = P_d.Phi();

        ////CM FRAME////

        TVector3 b;
        b = P_tot.BoostVector();    //boost to CM

        //beam and target//
        P_t.Boost(-b);  //to CM
        P_b.Boost(-b);  //to CM

        //gamma quanta//
        P_g1.Boost(-b); //to CM
        TVector3 vec_g1_cm = P_g1.Vect();

        Double_t p_g_cm[2], Theta_g_cm[2], Phi_g_cm[2];

        p_g_cm[0] = vec_g1_cm.Mag();
        Theta_g_cm[0] = P_g1.Theta();
        Phi_g_cm[0] = P_g1.Phi();

        P_g2.Boost(-b); //to CM
        TVector3 vec_g2_cm = P_g2.Vect();

        p_g_cm[1] = vec_g2_cm.Mag();
        Theta_g_cm[1] = P_g2.Theta();
        Phi_g_cm[1] = P_g2.Phi();

        Double_t OpeningAngle_g1_g2_cm = (vec_g1_cm.Angle(vec_g2_cm))*TMath::RadToDeg();    //angle between gamma quanta in CM

        //pion//
        P_pi0.Boost(-b);    //to CM
        TVector3 vec_pi0_cm = P_pi0.Vect();

        Double_t p_pi0_cm = vec_pi0_cm.Mag();
        Double_t E_pi0_cm = TMath::Sqrt(TMath::Power(p_pi0_cm,2) + TMath::Power(m_pi0,2));
        Double_t Ekin_pi0_cm = E_pi0_cm - m_pi0;
        Double_t Theta_pi0_cm = P_pi0.Theta();
        Double_t Phi_pi0_cm = P_pi0.Phi();

        //proton//
        P_p.Boost(-b);  //to CM
        TVector3 vec_p_cm = P_p.Vect();

        Double_t p_p_cm = vec_p_cm.Mag();
        Double_t E_p_cm = TMath::Sqrt(TMath::Power(p_p_cm,2) + TMath::Power(m_p,2));
        Double_t Ekin_p_cm = E_p_cm - m_p;
        Double_t Theta_p_cm = P_p.Theta();
        Double_t Phi_p_cm = P_p.Phi();

        Double_t OpeningAngle_pi0_p_cm = (vec_pi0_cm.Angle(vec_p_cm))*TMath::RadToDeg();    //angle between pion & proton in CM

        //deuteron//
        P_d.Boost(-b);  //to CM
        TVector3 vec_d_cm = P_d.Vect();

        Double_t p_d_cm = vec_d_cm.Mag();
        Double_t E_d_cm = TMath::Sqrt(TMath::Power(p_d_cm,2) + TMath::Power(m_d,2));
        Double_t Ekin_d_cm = E_d_cm - m_d;
        Double_t Theta_d_cm = P_d.Theta();
        Double_t Phi_d_cm = P_d.Phi();

        //Invariant and Missing Masses in CM//
        Double_t InvariantMass = (P_g1 + P_g2).M();
        Double_t MissingMass = (P_b + P_t - P_p - P_g1 - P_g2).M();
        Double_t MissingEnergy = (P_b + P_t - P_p - P_g1 - P_g2).E();
        Double_t TotalMass = (P_d + P_p + P_g1 + P_g2).M();

///////////////////////////////////////////MAIN CUTS////////////////////////////////////////////

        Bool_t lev[4][3];

        //////LEVEL 2//////
        //graphical cut on EdepPSBvsEdepSEC spectrum//

        lev[2][0] = (cutg0->IsInside(EdepSEC,EdepPSB));
        lev[2][1] = (cutg1->IsInside(EdepSEC,EdepPSB));
        lev[2][2] = (cutg2->IsInside(EdepSEC,EdepPSB));

        //////LEVEL 3//////
        //positively charged in CD && graphical cut on EdepPSBvsEdepSEC spectrum//

        lev[3][0] = ( (SgnMom > 0.) && cutg0->IsInside(EdepSEC,EdepPSB) );
        lev[3][1] = ( (SgnMom > 0.) && cutg1->IsInside(EdepSEC,EdepPSB) );
        lev[3][2] = ( (SgnMom > 0.) && cutg2->IsInside(EdepSEC,EdepPSB) );

        for (Int_t l = 2; l < 4; l++) {   //C01//

            if (lev[l][0]) {   //C02//

                ////FILLING HISTOGRAMS////

                if (l == 2) {hStatistics[1]->Fill(2);}
                if (l == 3) {hStatistics[2]->Fill(2);}

                hQ[l][0][0][0]->Fill(Q);

                hEdepFWC1vsFRH1[l][0]->Fill(EdepFWC1,EdepFRH1);
                hEdepFWC2vsFRH1[l][0]->Fill(EdepFWC2,EdepFRH1);
                hEdepFTH1vsFRH1[l][0]->Fill(EdepFTH1,EdepFRH1);
                hEdepFRH1vsFRH2[l][0]->Fill(EdepFRH1,EdepFRH2);
                hEdepFRH2vsFRH3[l][0]->Fill(EdepFRH2,EdepFRH3);
                hEdepFWC1vsFRH1FRH2FRH3[l][0]->Fill(EdepFWC1,(EdepFRH1 + EdepFRH2 + EdepFRH3));

                hEdepPSBvsSEC[l][0]->Fill(EdepSEC,EdepPSB);
                hEdepPSBvsSigMom[l][0]->Fill(SgnMom,EdepPSB);
                hEdepSECvsSigMom[l][0]->Fill(SgnMom,EdepSEC);

                //gamma quanta//
                hp_g1_lab[l][0]->Fill(p_g_lab[0]);
                hTheta_g1_lab[l][0]->Fill(Theta_g_lab[0]*TMath::RadToDeg());
                hPhi_g1_lab[l][0]->Fill(Phi_g_lab[0]*TMath::RadToDeg());
                hEkin_vs_Theta_g1_lab[l][0]->Fill(p_g_lab[0],Theta_g_lab[0]*TMath::RadToDeg());

                hp_g1_cm[l][0]->Fill(p_g_cm[0]);
                hTheta_g1_cm[l][0]->Fill(Theta_g_cm[0]*TMath::RadToDeg());
                hPhi_g1_cm[l][0]->Fill(Phi_g_cm[0]*TMath::RadToDeg());
                hEkin_vs_Theta_g1_cm[l][0]->Fill(p_g_cm[0],Theta_g_cm[0]*TMath::RadToDeg());

                hp_g2_lab[l][0]->Fill(p_g_lab[1]);
                hTheta_g2_lab[l][0]->Fill(Theta_g_lab[1]*TMath::RadToDeg());
                hPhi_g2_lab[l][0]->Fill(Phi_g_lab[1]*TMath::RadToDeg());
                hEkin_vs_Theta_g2_lab[l][0]->Fill(p_g_lab[1],Theta_g_lab[1]*TMath::RadToDeg());

                hp_g2_cm[l][0]->Fill(p_g_cm[1]);
                hTheta_g2_cm[l][0]->Fill(Theta_g_cm[1]*TMath::RadToDeg());
                hPhi_g2_cm[l][0]->Fill(Phi_g_cm[1]*TMath::RadToDeg());
                hEkin_vs_Theta_g2_cm[l][0]->Fill(p_g_cm[1],Theta_g_cm[1]*TMath::RadToDeg());

                hOpeningAngle_g1_g2_lab[l][0]->Fill(OpeningAngle_g1_g2_lab);
                hOpeningAngle_g1_g2_cm[l][0]->Fill(OpeningAngle_g1_g2_cm);

                hTheta_g1_vs_Theta_g2_lab[l][0]->Fill(Theta_g_lab[0]*TMath::RadToDeg(),Theta_g_lab[1]*TMath::RadToDeg());
                hTheta_g1_vs_Theta_g2_cm[l][0]->Fill(Theta_g_cm[0]*TMath::RadToDeg(),Theta_g_cm[1]*TMath::RadToDeg());

                //pion//
                hp_pi0_lab[l][0]->Fill(p_pi0_lab);
                hE_pi0_lab[l][0]->Fill(E_pi0_lab);
                hEkin_pi0_lab[l][0]->Fill(Ekin_pi0_lab);
                hTheta_pi0_lab[l][0]->Fill(Theta_pi0_lab*TMath::RadToDeg());
                hPhi_pi0_lab[l][0]->Fill(Phi_pi0_lab*TMath::RadToDeg());
                hEkin_vs_Theta_pi0_lab[l][0]->Fill(Ekin_pi0_lab,Theta_pi0_lab*TMath::RadToDeg());

                hp_pi0_cm[l][0]->Fill(p_pi0_cm);
                hE_pi0_cm[l][0]->Fill(E_pi0_cm);
                hEkin_pi0_cm[l][0]->Fill(Ekin_pi0_cm);
                hTheta_pi0_cm[l][0]->Fill(Theta_pi0_cm*TMath::RadToDeg());
                hPhi_pi0_cm[l][0]->Fill(Phi_pi0_cm*TMath::RadToDeg());
                hEkin_vs_Theta_pi0_cm[l][0]->Fill(Ekin_pi0_cm,Theta_pi0_cm*TMath::RadToDeg());

                //proton//
                hp_p_lab[l][0]->Fill(p_p_lab);
                hE_p_lab[l][0]->Fill(E_p_lab);
                hEkin_p_lab[l][0]->Fill(Ekin_p_lab);
                hTheta_p_lab[l][0]->Fill(Theta_p_lab*TMath::RadToDeg());
                hPhi_p_lab[l][0]->Fill(Phi_p_lab*TMath::RadToDeg());
                hEkin_vs_Theta_p_lab[l][0]->Fill(Ekin_p_lab,Theta_p_lab*TMath::RadToDeg());

                hp_p_cm[l][0]->Fill(p_p_cm);
                hE_p_cm[l][0]->Fill(E_p_cm);
                hEkin_p_cm[l][0]->Fill(Ekin_p_cm);
                hTheta_p_cm[l][0]->Fill(Theta_p_cm*TMath::RadToDeg());
                hPhi_p_cm[l][0]->Fill(Phi_p_cm*TMath::RadToDeg());
                hEkin_vs_Theta_p_cm[l][0]->Fill(Ekin_p_cm,Theta_p_cm*TMath::RadToDeg());

                hOpeningAngle_pi0_p_lab[l][0]->Fill(OpeningAngle_pi0_p_lab);
                hOpeningAngle_pi0_p_cm[l][0]->Fill(OpeningAngle_pi0_p_cm);

                //deuteron//
                hp_d_lab[l][0]->Fill(p_d_lab);
                hE_d_lab[l][0]->Fill(E_d_lab);
                hEkin_d_lab[l][0]->Fill(Ekin_d_lab);
                hTheta_d_lab[l][0]->Fill(Theta_d_lab*TMath::RadToDeg());
                hPhi_d_lab[l][0]->Fill(Phi_d_lab*TMath::RadToDeg());
                hEkin_vs_Theta_d_lab[l][0]->Fill(Ekin_d_lab,Theta_d_lab*TMath::RadToDeg());

                hp_d_cm[l][0]->Fill(p_d_cm);
                hE_d_cm[l][0]->Fill(E_d_cm);
                hEkin_d_cm[l][0]->Fill(Ekin_d_cm);
                hTheta_d_cm[l][0]->Fill(Theta_d_cm*TMath::RadToDeg());
                hPhi_d_cm[l][0]->Fill(Phi_d_cm*TMath::RadToDeg());
                hEkin_vs_Theta_d_cm[l][0]->Fill(Ekin_d_cm,Theta_d_cm*TMath::RadToDeg());

                if( (Theta_d_lab >= 0.052) && (Theta_d_lab <= 0.314) && (Theta_p_lab >= 0.349) && (Theta_p_lab <= 2.950) ) {
                    hTheta_p_vs_Theta_d_lab[l][0]->Fill(Theta_d_lab*TMath::RadToDeg(),Theta_p_lab*TMath::RadToDeg());
                }

                hTheta_p_vs_Theta_d_cm[l][0]->Fill(Theta_d_cm*TMath::RadToDeg(),Theta_p_cm*TMath::RadToDeg());

                //Invariant and Missing Masses in CM//
                hIM_pion[l][0]->Fill(InvariantMass);
                hMM_nucleon[l][0]->Fill(MissingMass);
                hME_nucleon[l][0]->Fill(MissingEnergy);
                hMMvsME_nucleon[l][0]->Fill(MissingMass,MissingEnergy);
                hIMvsMM[l][0]->Fill(InvariantMass,MissingMass);

////////////////////////////////////////////////////////////////////////////////////////////////

                //CONDITIONS FOR CUTS
                Bool_t cond[5][3];

                //Invariant Mass
                Double_t beginCutIM[3],endCutIM[3],deltaCutIM[1];
                beginCutIM[0] = 0.1;
                endCutIM[0] = 0.17;
                deltaCutIM[0] = 0.003;

                beginCutIM[1] = beginCutIM[0] - deltaCutIM[0];
                endCutIM[1] = endCutIM[0] + deltaCutIM[0];
                beginCutIM[2] = beginCutIM[0] + deltaCutIM[0];
                endCutIM[2] = endCutIM[0] - deltaCutIM[0];

                cond[1][0] = ( (InvariantMass >= beginCutIM[0]) && (InvariantMass <= endCutIM[0]) );
                cond[1][1] = ( (InvariantMass >= beginCutIM[1]) && (InvariantMass <= endCutIM[1]) );
                cond[1][2] = ( (InvariantMass >= beginCutIM[2]) && (InvariantMass <= endCutIM[2]) );

                //Open Angle between pi0 & proton
                Double_t beginCutOA[3],deltaCutOA[1];
                beginCutOA[0] = 155.;
                deltaCutOA[0] = 1.;

                beginCutOA[1] = beginCutOA[0] - deltaCutOA[0];
                beginCutOA[2] = beginCutOA[0] + deltaCutOA[0];

                cond[2][0] = ( (OpeningAngle_pi0_p_cm >= beginCutOA[0]) );
                cond[2][1] = ( (OpeningAngle_pi0_p_cm >= beginCutOA[1]) );
                cond[2][2] = ( (OpeningAngle_pi0_p_cm >= beginCutOA[2]) );

                //Missing Mass
                Double_t beginCutMM[3],endCutMM[3],deltaCutMM[1];
                beginCutMM[0] = 1.7;
                endCutMM[0] = 2.05;
                deltaCutMM[0] = 0.01;

                beginCutMM[1] = beginCutMM[0] - deltaCutMM[0];
                endCutMM[1] = endCutMM[0] + deltaCutMM[0];
                beginCutMM[2] = beginCutMM[0] + deltaCutMM[0];
                endCutMM[2] = endCutMM[0] - deltaCutMM[0];

                cond[3][0] = ( (MissingMass >= beginCutMM[0]) && (MissingMass <= endCutMM[0]) );
                cond[3][1] = ( (MissingMass >= beginCutMM[1]) && (MissingMass <= endCutMM[1]) );
                cond[3][2] = ( (MissingMass >= beginCutMM[2]) && (MissingMass <= endCutMM[2]) );

                //Deuteron Momentum
                Double_t beginCutDM[3],endCutDM[3],deltaCutDM[1];
                beginCutDM[0] = 0.6;
                endCutDM[0] = 1.1;
                deltaCutDM[0] = 0.02;

                beginCutDM[1] = beginCutDM[0] - deltaCutDM[0];
                endCutDM[1] = endCutDM[0] + deltaCutDM[0];
                beginCutDM[2] = beginCutDM[0] + deltaCutDM[0];
                endCutDM[2] = endCutDM[0] - deltaCutDM[0];

                cond[4][0] = ( (p_d_lab >= beginCutDM[0]) && (p_d_lab <= endCutDM[0]) );
                cond[4][1] = ( (p_d_lab >= beginCutDM[1]) && (p_d_lab <= endCutDM[1]) );
                cond[4][2] = ( (p_d_lab >= beginCutDM[2]) && (p_d_lab <= endCutDM[2]) );

                //Deuteron Kinetic Energy in CM
                Double_t beginCutDE[3],deltaCutDE[1];
                beginCutDE[0] = 0.03;
                deltaCutDE[0] = 0.002;

                beginCutDE[1] = beginCutDE[0] - deltaCutDE[0];
                beginCutDE[2] = beginCutDE[0] + deltaCutDE[0];

                cond[5][0] = ( (Ekin_d_cm <= beginCutDE[0]) );
                cond[5][1] = ( (Ekin_d_cm <= beginCutDE[1]) );
                cond[5][2] = ( (Ekin_d_cm <= beginCutDE[2]) );

                /////////////////////////////////
                Bool_t cut[7][10][10];
                cut[1][0][0] = ( cond[1][0] );                                             //cut on Invariant Mass
                cut[2][0][0] = ( cond[1][0] && cond[2][0] );                               //cut1 + cut on pi0-p opening angle in CM
                cut[3][0][0] = ( cond[1][0] && cond[2][0] && cond[3][0] );                 //cut2 + cut on Missing Mass
                cut[4][0][0] = ( cond[1][0] && cond[2][0] && cond[3][0] && cond[4][0] );   //cut3 + cut on deuteron momentum
                cut[5][0][0] = ( cond[1][0] && cond[2][0] && cond[3][0] && cond[5][0] );   //cut3 + cut on deuteron kinetic energy
                cut[6][0][0] = ( cut[4][0][0] && cond[5][0] );                             //cut4 + cut on deuteron kinetic energy

                for (Int_t k = 1; k < 7; k++) { //C03//

                    if (cut[k][0][0]) {               //C04//

                        if ( (l == 2) && (k < 5) ) {hStatistics[1]->Fill(2+k);}
                        if ( (l == 3) && (k < 5) ) {hStatistics[2]->Fill(2+k);}

                        hQ[l][k][0][0]->Fill(Q);

                        hEdepFWC1vsFRH1[l][k]->Fill(EdepFWC1,EdepFRH1);
                        hEdepFWC2vsFRH1[l][k]->Fill(EdepFWC2,EdepFRH1);
                        hEdepFTH1vsFRH1[l][k]->Fill(EdepFTH1,EdepFRH1);
                        hEdepFRH1vsFRH2[l][k]->Fill(EdepFRH1,EdepFRH2);
                        hEdepFRH2vsFRH3[l][k]->Fill(EdepFRH2,EdepFRH3);
                        hEdepFWC1vsFRH1FRH2FRH3[l][k]->Fill(EdepFWC1,(EdepFRH1 + EdepFRH2 + EdepFRH3));

                        hEdepPSBvsSEC[l][k]->Fill(EdepSEC,EdepPSB);
                        hEdepPSBvsSigMom[l][k]->Fill(SgnMom,EdepPSB);
                        hEdepSECvsSigMom[l][k]->Fill(SgnMom,EdepSEC);

                        //gamma quanta//
                        hp_g1_lab[l][k]->Fill(p_g_lab[0]);
                        hTheta_g1_lab[l][k]->Fill(Theta_g_lab[0]*TMath::RadToDeg());
                        hPhi_g1_lab[l][k]->Fill(Phi_g_lab[0]*TMath::RadToDeg());
                        hEkin_vs_Theta_g1_lab[l][k]->Fill(p_g_lab[0],Theta_g_lab[0]*TMath::RadToDeg());

                        hp_g1_cm[l][k]->Fill(p_g_cm[0]);
                        hTheta_g1_cm[l][k]->Fill(Theta_g_cm[0]*TMath::RadToDeg());
                        hPhi_g1_cm[l][k]->Fill(Phi_g_cm[0]*TMath::RadToDeg());
                        hEkin_vs_Theta_g1_cm[l][k]->Fill(p_g_cm[0],Theta_g_cm[0]*TMath::RadToDeg());

                        hp_g2_lab[l][k]->Fill(p_g_lab[1]);
                        hTheta_g2_lab[l][k]->Fill(Theta_g_lab[1]*TMath::RadToDeg());
                        hPhi_g2_lab[l][k]->Fill(Phi_g_lab[1]*TMath::RadToDeg());
                        hEkin_vs_Theta_g2_lab[l][k]->Fill(p_g_lab[1],Theta_g_lab[1]*TMath::RadToDeg());

                        hp_g2_cm[l][k]->Fill(p_g_cm[1]);
                        hTheta_g2_cm[l][k]->Fill(Theta_g_cm[1]*TMath::RadToDeg());
                        hPhi_g2_cm[l][k]->Fill(Phi_g_cm[1]*TMath::RadToDeg());
                        hEkin_vs_Theta_g2_cm[l][k]->Fill(p_g_cm[1],Theta_g_cm[1]*TMath::RadToDeg());

                        hOpeningAngle_g1_g2_lab[l][k]->Fill(OpeningAngle_g1_g2_lab);
                        hOpeningAngle_g1_g2_cm[l][k]->Fill(OpeningAngle_g1_g2_cm);

                        hTheta_g1_vs_Theta_g2_lab[l][k]->Fill(Theta_g_lab[0]*TMath::RadToDeg(),Theta_g_lab[1]*TMath::RadToDeg());
                        hTheta_g1_vs_Theta_g2_cm[l][k]->Fill(Theta_g_cm[0]*TMath::RadToDeg(),Theta_g_cm[1]*TMath::RadToDeg());

                        //pion//
                        hp_pi0_lab[l][k]->Fill(p_pi0_lab);
                        hE_pi0_lab[l][k]->Fill(E_pi0_lab);
                        hEkin_pi0_lab[l][k]->Fill(Ekin_pi0_lab);
                        hTheta_pi0_lab[l][k]->Fill(Theta_pi0_lab*TMath::RadToDeg());
                        hPhi_pi0_lab[l][k]->Fill(Phi_pi0_lab*TMath::RadToDeg());
                        hEkin_vs_Theta_pi0_lab[l][k]->Fill(Ekin_pi0_lab,Theta_pi0_lab*TMath::RadToDeg());

                        hp_pi0_cm[l][k]->Fill(p_pi0_cm);
                        hE_pi0_cm[l][k]->Fill(E_pi0_cm);
                        hEkin_pi0_cm[l][k]->Fill(Ekin_pi0_cm);
                        hTheta_pi0_cm[l][k]->Fill(Theta_pi0_cm*TMath::RadToDeg());
                        hPhi_pi0_cm[l][k]->Fill(Phi_pi0_cm*TMath::RadToDeg());
                        hEkin_vs_Theta_pi0_cm[l][k]->Fill(Ekin_pi0_cm,Theta_pi0_cm*TMath::RadToDeg());

                        //proton//
                        hp_p_lab[l][k]->Fill(p_p_lab);
                        hE_p_lab[l][k]->Fill(E_p_lab);
                        hEkin_p_lab[l][k]->Fill(Ekin_p_lab);
                        hTheta_p_lab[l][k]->Fill(Theta_p_lab*TMath::RadToDeg());
                        hPhi_p_lab[l][k]->Fill(Phi_p_lab*TMath::RadToDeg());
                        hEkin_vs_Theta_p_lab[l][k]->Fill(Ekin_p_lab,Theta_p_lab*TMath::RadToDeg());

                        hp_p_cm[l][k]->Fill(p_p_cm);
                        hE_p_cm[l][k]->Fill(E_p_cm);
                        hEkin_p_cm[l][k]->Fill(Ekin_p_cm);
                        hTheta_p_cm[l][k]->Fill(Theta_p_cm*TMath::RadToDeg());
                        hPhi_p_cm[l][k]->Fill(Phi_p_cm*TMath::RadToDeg());
                        hEkin_vs_Theta_p_cm[l][k]->Fill(Ekin_p_cm,Theta_p_cm*TMath::RadToDeg());

                        hOpeningAngle_pi0_p_lab[l][k]->Fill(OpeningAngle_pi0_p_lab);
                        hOpeningAngle_pi0_p_cm[l][k]->Fill(OpeningAngle_pi0_p_cm);

                        //deuteron//
                        hp_d_lab[l][k]->Fill(p_d_lab);
                        hE_d_lab[l][k]->Fill(E_d_lab);
                        hEkin_d_lab[l][k]->Fill(Ekin_d_lab);
                        hTheta_d_lab[l][k]->Fill(Theta_d_lab*TMath::RadToDeg());
                        hPhi_d_lab[l][k]->Fill(Phi_d_lab*TMath::RadToDeg());
                        hEkin_vs_Theta_d_lab[l][k]->Fill(Ekin_d_lab,Theta_d_lab*TMath::RadToDeg());

                        hp_d_cm[l][k]->Fill(p_d_cm);
                        hE_d_cm[l][k]->Fill(E_d_cm);
                        hEkin_d_cm[l][k]->Fill(Ekin_d_cm);
                        hTheta_d_cm[l][k]->Fill(Theta_d_cm*TMath::RadToDeg());
                        hPhi_d_cm[l][k]->Fill(Phi_d_cm*TMath::RadToDeg());
                        hEkin_vs_Theta_d_cm[l][k]->Fill(Ekin_d_cm,Theta_d_cm*TMath::RadToDeg());

                        if( (Theta_d_lab >= 0.052) && (Theta_d_lab <= 0.314) && (Theta_p_lab >= 0.349) && (Theta_p_lab <= 2.950) ) {
                            hTheta_p_vs_Theta_d_lab[l][k]->Fill(Theta_d_lab*TMath::RadToDeg(),Theta_p_lab*TMath::RadToDeg());
                        }

                        hTheta_p_vs_Theta_d_cm[l][k]->Fill(Theta_d_cm*TMath::RadToDeg(),Theta_p_cm*TMath::RadToDeg());

                        //Invariant and Missing Masses in CM//
                        hIM_pion[l][k]->Fill(InvariantMass);
                        hMM_nucleon[l][k]->Fill(MissingMass);
                        hME_nucleon[l][k]->Fill(MissingEnergy);
                        hMMvsME_nucleon[l][k]->Fill(MissingMass,MissingEnergy);
                        hIMvsMM[l][k]->Fill(InvariantMass,MissingMass);

                    }   //C04//

                }       //C03//

                //SYSTEMATICS//
                cut[4][1][1] = ( cond[1][1] && cond[2][0] && cond[3][0] && cond[4][0] );
                cut[4][1][2] = ( cond[1][2] && cond[2][0] && cond[3][0] && cond[4][0] );
                cut[4][2][1] = ( cond[1][0] && cond[2][1] && cond[3][0] && cond[4][0] );
                cut[4][2][2] = ( cond[1][0] && cond[2][2] && cond[3][0] && cond[4][0] );
                cut[4][3][1] = ( cond[1][0] && cond[2][0] && cond[3][1] && cond[4][0] );
                cut[4][3][2] = ( cond[1][0] && cond[2][0] && cond[3][2] && cond[4][0] );
                cut[4][4][1] = ( cond[1][0] && cond[2][0] && cond[3][0] && cond[4][1] );
                cut[4][4][2] = ( cond[1][0] && cond[2][0] && cond[3][0] && cond[4][2] );

                cut[5][1][1] = ( cond[1][1] && cond[2][0] && cond[3][0] && cond[5][0] );
                cut[5][1][2] = ( cond[1][2] && cond[2][0] && cond[3][0] && cond[5][0] );
                cut[5][2][1] = ( cond[1][0] && cond[2][1] && cond[3][0] && cond[5][0] );
                cut[5][2][2] = ( cond[1][0] && cond[2][2] && cond[3][0] && cond[5][0] );
                cut[5][3][1] = ( cond[1][0] && cond[2][0] && cond[3][1] && cond[5][0] );
                cut[5][3][2] = ( cond[1][0] && cond[2][0] && cond[3][2] && cond[5][0] );
                cut[5][4][1] = ( cond[1][0] && cond[2][0] && cond[3][0] && cond[5][1] );
                cut[5][4][2] = ( cond[1][0] && cond[2][0] && cond[3][0] && cond[5][2] );

                for (Int_t m = 4; m < 6; m++) {
                    for (Int_t n = 1; n < 5; n++) {
                        for (Int_t p = 1; p < 3; p++) {
                            if (cut[m][n][p]) {
                                hQ[l][m][n][p]->Fill(Q);
                            }
                        }
                    }
                }

            }   //C02//

        }   //C01//

        for (Int_t q = 2; q < 4; q++) { //D01//
            for (Int_t r = 1; r < 3; r++) { //D02//
                if(lev[q][r]) {   //D03//

                    if ( (InvariantMass >= 0.1) && (InvariantMass <= 0.17) ) {
                        if (OpeningAngle_pi0_p_cm >= 155.) {
                            if ( (MissingMass >= 1.7) && (MissingMass <= 2.05) ) {
                                if ( (p_d_lab >= 0.6) && (p_d_lab <= 1.1) ) {
                                    hQ[q+1][r][4][0]->Fill(Q);
                                }
                                if ( Ekin_d_cm <= 0.03 ) {
                                    hQ[q+1][r][5][0]->Fill(Q);
                                }
                            }
                        }
                    }
                }   //D03//
            }   //D02//
        }   //D01//

    }   //B01//

    return;

}   //01//

////////////////////////////////////////////////////////////////////////////////////////////////

void eventselection::SetupSpectra(const char * lpath) {   //02//

    TString hpathMC = "WMC";
    TString h_st = "Statistics";

    TString hpath[10][10];
    TString hpathP[10][10];
    TString hpathD[10][10];
    TString hpathG[10][10];
    TString hpathPi[10][10];

    for (Int_t i = 0; i < 10; i++) {
        for (Int_t j = 0; j < 10; j++) {
            hpath[i][j] = Form("DATA_lev%d_cut%d",i,j);
            hpathP[i][j] = hpath[i][j] + Form("/proton");
            hpathD[i][j] = hpath[i][j] + Form("/deuteron");
            hpathG[i][j] = hpath[i][j] + Form("/gammas");
            hpathPi[i][j] = hpath[i][j] + Form("/pion");
        }
    }

////////////////////////////////////TRUE EVENTS (MONTE CARLO////////////////////////////////////

    hp_beam_MC = new TH1F("hp_beam_MC","",200,1.41,1.65);
    hp_beam_MC->GetXaxis()->SetTitle("p_{beam} [GeV/c]");
    hp_beam_MC->GetYaxis()->SetTitle("N");
    gHistoManager->Add(hp_beam_MC,hpathMC);

    hGenerated_Q = new TH1F("hGenerated_Q","",40,-70.,30.);
    hGenerated_Q->GetXaxis()->SetTitle("Q [MeV]");
    hGenerated_Q->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hGenerated_Q,hpathMC);

    hAccepted_Q = new TH1F("hAccepted_Q","",40,-70.,30.);
    hAccepted_Q->GetXaxis()->SetTitle("Q [MeV]");
    hAccepted_Q->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hAccepted_Q,hpathMC);

    //deuteron//
    hEkin_d_lab_MC = new TH1F("hEkin_d_lab_MC","",500,0.,0.4);
    hEkin_d_lab_MC->GetXaxis()->SetTitle("E^{kin}_{d} [GeV]");
    hEkin_d_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_d_lab_MC,hpathMC);

    hEkin_d_cm_MC = new TH1F("hEkin_d_cm_MC","",500,0.,0.02);
    hEkin_d_cm_MC->GetXaxis()->SetTitle("E^{kin}_{d} [GeV]");
    hEkin_d_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_d_cm_MC,hpathMC);

    hE_d_lab_MC = new TH1F("hE_d_lab_MC","",500,1.8,2.8);
    hE_d_lab_MC->GetXaxis()->SetTitle("E^{lab}_{d} [GeV]");
    hE_d_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_d_lab_MC,hpathMC);

    hE_d_cm_MC = new TH1F("hE_d_cm_MC","",500,1.86,1.91);
    hE_d_cm_MC->GetXaxis()->SetTitle("E^{cm}_{d} [GeV]");
    hE_d_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_d_cm_MC,hpathMC);

    hp_d_lab_MC = new TH1F("hp_d_lab_MC","",500,0.,2.);
    hp_d_lab_MC->GetXaxis()->SetTitle("p^{lab}_{d} [GeV/c]");
    hp_d_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_d_lab_MC,hpathMC);

    hp_d_cm_MC = new TH1F("hp_d_cm_MC","",500,0.,0.4);
    hp_d_cm_MC->GetXaxis()->SetTitle("p^{cm}_{d} [GeV/c]");
    hp_d_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_d_cm_MC,hpathMC);

    hTheta_d_lab_MC = new TH1F("hTheta_d_lab_MC","",500,0.,20.);
    hTheta_d_lab_MC->GetXaxis()->SetTitle("#theta^{lab}_{d} [deg]");
    hTheta_d_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_d_lab_MC,hpathMC);

    hTheta_d_cm_MC = new TH1F("hTheta_d_cm_MC","",360,0.,180.);
    hTheta_d_cm_MC->GetXaxis()->SetTitle("#theta^{cm}_{d} [deg]");
    hTheta_d_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_d_cm_MC,hpathMC);

    hPhi_d_lab_MC = new TH1F("hPhi_d_lab_MC","",360,0.,360.);
    hPhi_d_lab_MC->GetXaxis()->SetTitle("#phi^{lab}_{d} [deg]");
    hPhi_d_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_d_lab_MC,hpathMC);

    hPhi_d_cm_MC = new TH1F("hPhi_d_cm_MC","",360,-180.,180.);
    hPhi_d_cm_MC->GetXaxis()->SetTitle("#phi^{cm}_{d} [deg]");
    hPhi_d_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_d_cm_MC,hpathMC);

    hEkin_vs_Theta_d_lab_MC = new TH2F("hEkin_vs_Theta_d_lab_MC","",500,0.,0.4,500,0.,20.);
    hEkin_vs_Theta_d_lab_MC->GetXaxis()->SetTitle("E^{kin}_{d} [GeV]");
    hEkin_vs_Theta_d_lab_MC->GetYaxis()->SetTitle("#theta^{lab}_{d} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_d_lab_MC,hpathMC);

    hEkin_vs_Theta_d_lab_acc_MC = new TH2F("hEkin_vs_Theta_d_lab_acc_MC","",500,0.,0.4,500,0.,20.);
    hEkin_vs_Theta_d_lab_acc_MC->GetXaxis()->SetTitle("E^{kin}_{d} [GeV]");
    hEkin_vs_Theta_d_lab_acc_MC->GetYaxis()->SetTitle("#theta^{lab}_{d} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_d_lab_acc_MC,hpathMC);

    hEkin_vs_Theta_d_cm_MC = new TH2F("hEkin_vs_Theta_d_cm_MC","",500,0.,0.05,500,0.,180.);
    hEkin_vs_Theta_d_cm_MC->GetXaxis()->SetTitle("E^{kin,cm}_{d} [GeV]");
    hEkin_vs_Theta_d_cm_MC->GetYaxis()->SetTitle("#theta^{cm}_{d} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_d_cm_MC,hpathMC);

    //proton//
    hEkin_p_lab_MC = new TH1F("hEkin_p_lab_MC","",500,0.,0.6);
    hEkin_p_lab_MC->GetXaxis()->SetTitle("E^{kin}_{p} [GeV]");
    hEkin_p_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_p_lab_MC,hpathMC);

    hEkin_p_cm_MC = new TH1F("hEkin_p_cm_MC","",500,0.,0.25);
    hEkin_p_cm_MC->GetXaxis()->SetTitle("E^{kin}_{p} [GeV]");
    hEkin_p_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_p_cm_MC,hpathMC);

    hE_p_lab_MC = new TH1F("hE_p_lab_MC","",500,0.9,1.5);
    hE_p_lab_MC->GetXaxis()->SetTitle("E^{lab}_{p} [GeV]");
    hE_p_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_p_lab_MC,hpathMC);

    hE_p_cm_MC = new TH1F("hE_p_cm_MC","",500,0.9,1.25);
    hE_p_cm_MC->GetXaxis()->SetTitle("E^{cm}_{p} [GeV]");
    hE_p_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_p_cm_MC,hpathMC);

    hp_p_lab_MC = new TH1F("hp_p_lab_MC","",750,0.,1.5);
    hp_p_lab_MC->GetXaxis()->SetTitle("p^{lab}_{p} [GeV/c]");
    hp_p_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_p_lab_MC,hpathMC);

    hp_p_cm_MC = new TH1F("hp_p_cm_MC","",500,0.,1.);
    hp_p_cm_MC->GetXaxis()->SetTitle("p^{cm}_{p} [GeV/c]");
    hp_p_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_p_cm_MC,hpathMC);

    hTheta_p_lab_MC = new TH1F("hTheta_p_lab_MC","",360,0.,180.);
    hTheta_p_lab_MC->GetXaxis()->SetTitle("#theta^{lab}_{p} [deg]");
    hTheta_p_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_p_lab_MC,hpathMC);

    hTheta_p_cm_MC = new TH1F("hTheta_p_cm_MC","",360,0.,180.);
    hTheta_p_cm_MC->GetXaxis()->SetTitle("#theta^{cm}_{p} [deg]");
    hTheta_p_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_p_cm_MC,hpathMC);

    hPhi_p_lab_MC = new TH1F("hPhi_p_lab_MC","",360,0.,360.);
    hPhi_p_lab_MC->GetXaxis()->SetTitle("#phi^{lab}_{p} [deg]");
    hPhi_p_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_p_lab_MC,hpathMC);

    hPhi_p_cm_MC = new TH1F("hPhi_p_cm_MC","",360,-180.,180.);
    hPhi_p_cm_MC->GetXaxis()->SetTitle("#phi^{cm}_{p} [deg]");
    hPhi_p_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_p_cm_MC,hpathMC);

    hEkin_vs_Theta_p_lab_MC = new TH2F("hEkin_vs_Theta_p_lab_MC","",500,0.,0.6,500,0.,180.);
    hEkin_vs_Theta_p_lab_MC->GetXaxis()->SetTitle("E^{kin}_{p} [GeV]");
    hEkin_vs_Theta_p_lab_MC->GetYaxis()->SetTitle("#theta^{lab}_{p} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_p_lab_MC,hpathMC);

    hEkin_vs_Theta_p_lab_acc_MC = new TH2F("hEkin_vs_Theta_p_lab_acc_MC","",500,0.,0.6,500,0.,180.);
    hEkin_vs_Theta_p_lab_acc_MC->GetXaxis()->SetTitle("E^{kin}_{p} [GeV]");
    hEkin_vs_Theta_p_lab_acc_MC->GetYaxis()->SetTitle("#theta^{lab}_{p} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_p_lab_acc_MC,hpathMC);

    hEkin_vs_Theta_p_cm_MC = new TH2F("hEkin_vs_Theta_p_cm_MC","",500,0.,0.25,500,0.,180.);
    hEkin_vs_Theta_p_cm_MC->GetXaxis()->SetTitle("E^{kin,cm}_{p} [GeV]");
    hEkin_vs_Theta_p_cm_MC->GetYaxis()->SetTitle("#theta^{cm}_{p} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_p_cm_MC,hpathMC);

    ////
    hTheta_p_vs_Theta_d_lab_MC = new TH2F("hTheta_p_vs_Theta_d_lab_MC","",500,0.,20.,500,0.,180.);
    hTheta_p_vs_Theta_d_lab_MC->GetXaxis()->SetTitle("#theta^{lab}_{d} [deg]");
    hTheta_p_vs_Theta_d_lab_MC->GetYaxis()->SetTitle("#theta^{lab}_{p} [deg]");
    gHistoManager->Add(hTheta_p_vs_Theta_d_lab_MC,hpathMC);

    hTheta_p_vs_Theta_d_lab_cut_MC = new TH2F("hTheta_p_vs_Theta_d_lab_cut_MC","",500,0.,20.,500,0.,180.);
    hTheta_p_vs_Theta_d_lab_cut_MC->GetXaxis()->SetTitle("#theta^{lab}_{d} [deg]");
    hTheta_p_vs_Theta_d_lab_cut_MC->GetYaxis()->SetTitle("#theta^{lab}_{p} [deg]");
    gHistoManager->Add(hTheta_p_vs_Theta_d_lab_cut_MC,hpathMC);

    hTheta_p_vs_Theta_d_cm_MC = new TH2F("hTheta_p_vs_Theta_d_cm_MC","",540,0.,180.,540,0.,180.);
    hTheta_p_vs_Theta_d_cm_MC->GetXaxis()->SetTitle("#theta^{lab}_{d} [deg]");
    hTheta_p_vs_Theta_d_cm_MC->GetYaxis()->SetTitle("#theta^{lab}_{p} [deg]");
    gHistoManager->Add(hTheta_p_vs_Theta_d_cm_MC,hpathMC);

    //gamma #1//
    //hEkin_g1_lab_MC = new TH1F("hEkin_g1_lab_MC","",500,0.,0.8);
    //hEkin_g1_lab_MC->GetXaxis()->SetTitle("E^{kin}_{#gamma_{1}} [GeV/c]");
    //hEkin_g1_lab_MC->GetYaxis()->SetTitle("counts");
    //gHistoManager->Add(hEkin_g1_lab_MC,hpathMC);

    hp_g1_lab_MC = new TH1F("hp_g1_lab_MC","",500,0.,0.8);
    hp_g1_lab_MC->GetXaxis()->SetTitle("p^{lab}_{#gamma_{1}} [GeV/c]");
    hp_g1_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_g1_lab_MC,hpathMC);

    hp_g1_cm_MC = new TH1F("hp_g1_cm_MC","",500,0.,0.6);
    hp_g1_cm_MC->GetXaxis()->SetTitle("p^{cm}_{#gamma_{1}} [GeV/c]");
    hp_g1_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_g1_cm_MC,hpathMC);

    //hE_g1_cm_MC = new TH1F("hE_g1_cm_MC","",500,0.,0.6);
    //hE_g1_cm_MC->GetXaxis()->SetTitle("E^{cm}_{#gamma_{1}} [GeV]");
    //hE_g1_cm_MC->GetYaxis()->SetTitle("counts");
    //gHistoManager->Add(hE_g1_cm_MC,hpathMC);

    hTheta_g1_lab_MC = new TH1F("hTheta_g1_lab_MC","",360,0.,180.);
    hTheta_g1_lab_MC->GetXaxis()->SetTitle("#theta^{lab}_{#gamma_{1}} [deg]");
    hTheta_g1_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_g1_lab_MC,hpathMC);

    hTheta_g1_cm_MC = new TH1F("hTheta_g1_cm_MC","",360,0.,180.);
    hTheta_g1_cm_MC->GetXaxis()->SetTitle("#theta^{cm}_{#gamma_{1}} [deg]");
    hTheta_g1_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_g1_cm_MC,hpathMC);

    hPhi_g1_lab_MC = new TH1F("hPhi_g1_lab_MC","",360,0.,360.);
    hPhi_g1_lab_MC->GetXaxis()->SetTitle("#phi^{lab}_{#gamma_{1}} [deg]");
    hPhi_g1_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_g1_lab_MC,hpathMC);

    hPhi_g1_cm_MC = new TH1F("hPhi_g1_cm_MC","",360,-180.,180.);
    hPhi_g1_cm_MC->GetXaxis()->SetTitle("#phi^{cm}_{#gamma_{1}} [deg]");
    hPhi_g1_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_g1_cm_MC,hpathMC);

    hEkin_vs_Theta_g1_lab_MC = new TH2F("hEkin_vs_Theta_g1_lab_MC","",500,0.,0.85,500,0.,180.);
    hEkin_vs_Theta_g1_lab_MC->GetXaxis()->SetTitle("E^{kin}_{#gamma_{1}} [GeV]");
    hEkin_vs_Theta_g1_lab_MC->GetYaxis()->SetTitle("#theta^{lab}_{#gamma_{1}} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_g1_lab_MC,hpathMC);

    hEkin_vs_Theta_g1_lab_acc_MC = new TH2F("hEkin_vs_Theta_g1_lab_acc_MC","",500,0.,0.85,500,0.,180.);
    hEkin_vs_Theta_g1_lab_acc_MC->GetXaxis()->SetTitle("E^{kin}_{#gamma_{1}} [GeV]");
    hEkin_vs_Theta_g1_lab_acc_MC->GetYaxis()->SetTitle("#theta^{lab}_{#gamma_{1}} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_g1_lab_acc_MC,hpathMC);

    hEkin_vs_Theta_g1_cm_MC = new TH2F("hEkin_vs_Theta_g1_cm_MC","",500,0.,0.6,500,0.,180.);
    hEkin_vs_Theta_g1_cm_MC->GetXaxis()->SetTitle("E^{kin,cm}_{#gamma_{1}} [GeV]");
    hEkin_vs_Theta_g1_cm_MC->GetYaxis()->SetTitle("#theta^{cm}_{#gamma_{1}} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_g1_cm_MC,hpathMC);

    //gamma #2//
    //hEkin_g2_lab_MC = new TH1F("hEkin_g2_lab_MC","",500,0.,0.8);
    //hEkin_g2_lab_MC->GetXaxis()->SetTitle("E^{kin}_{#gamma_{2}} [GeV/c]");
    //hEkin_g2_lab_MC->GetYaxis()->SetTitle("counts");
    //gHistoManager->Add(hEkin_g2_lab_MC,hpathMC);

    hp_g2_lab_MC = new TH1F("hp_g2_lab_MC","",500,0.,0.8);
    hp_g2_lab_MC->GetXaxis()->SetTitle("p^{lab}_{#gamma_{2}} [GeV/c]");
    hp_g2_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_g2_lab_MC,hpathMC);

    hp_g2_cm_MC = new TH1F("hp_g2_cm_MC","",500,0.,0.6);
    hp_g2_cm_MC->GetXaxis()->SetTitle("p^{cm}_{#gamma_{2}} [GeV/c]");
    hp_g2_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_g2_cm_MC,hpathMC);

    //hE_g2_cm_MC = new TH1F("hE_g2_cm_MC","",500,0.,0.6);
    //hE_g2_cm_MC->GetXaxis()->SetTitle("E^{cm}_{#gamma_{2}} [GeV]");
    //hE_g2_cm_MC->GetYaxis()->SetTitle("counts");
    //gHistoManager->Add(hE_g2_cm_MC,hpathMC);

    hTheta_g2_lab_MC = new TH1F("hTheta_g2_lab_MC","",360,0.,180.);
    hTheta_g2_lab_MC->GetXaxis()->SetTitle("#theta^{lab}_{#gamma_{2}} [deg]");
    hTheta_g2_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_g2_lab_MC,hpathMC);

    hTheta_g2_cm_MC = new TH1F("hTheta_g2_cm_MC","",360,0.,180.);
    hTheta_g2_cm_MC->GetXaxis()->SetTitle("#theta^{cm}_{#gamma_{2}} [deg]");
    hTheta_g2_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_g2_cm_MC,hpathMC);

    hPhi_g2_lab_MC = new TH1F("hPhi_g2_lab_MC","",360,0.,360.);
    hPhi_g2_lab_MC->GetXaxis()->SetTitle("#phi^{lab}_{#gamma_{2}} [deg]");
    hPhi_g2_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_g2_lab_MC,hpathMC);

    hPhi_g2_cm_MC = new TH1F("hPhi_g2_cm_MC","",360,-180.,180.);
    hPhi_g2_cm_MC->GetXaxis()->SetTitle("#phi^{cm}_{#gamma_{1}} [deg]");
    hPhi_g2_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_g2_cm_MC,hpathMC);

    hEkin_vs_Theta_g2_lab_MC = new TH2F("hEkin_vs_Theta_g2_lab_MC","",500,0.,0.85,500,0.,180.);
    hEkin_vs_Theta_g2_lab_MC->GetXaxis()->SetTitle("E^{kin}_{#gamma_{2}} [GeV]");
    hEkin_vs_Theta_g2_lab_MC->GetYaxis()->SetTitle("#theta^{lab}_{#gamma_{2}} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_g2_lab_MC,hpathMC);

    hEkin_vs_Theta_g2_lab_acc_MC = new TH2F("hEkin_vs_Theta_g2_lab_acc_MC","",500,0.,0.85,500,0.,180.);
    hEkin_vs_Theta_g2_lab_acc_MC->GetXaxis()->SetTitle("E^{kin}_{#gamma_{2}} [GeV]");
    hEkin_vs_Theta_g2_lab_acc_MC->GetYaxis()->SetTitle("#theta^{lab}_{#gamma_{2}} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_g2_lab_acc_MC,hpathMC);

    hEkin_vs_Theta_g2_cm_MC = new TH2F("hEkin_vs_Theta_g2_cm_MC","",500,0.,0.6,500,0.,180.);
    hEkin_vs_Theta_g2_cm_MC->GetXaxis()->SetTitle("E^{kin,cm}_{#gamma_{2}} [GeV]");
    hEkin_vs_Theta_g2_cm_MC->GetYaxis()->SetTitle("#theta^{cm}_{#gamma_{2}} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_g2_cm_MC,hpathMC);

    ////
    hEkin_vs_Theta_gammas_lab_MC = new TH2F("hEkin_vs_Theta_gammas_lab_MC","",500,0.,0.85,500,0.,180.);
    hEkin_vs_Theta_gammas_lab_MC->GetXaxis()->SetTitle("E^{kin}_{#gamma} [GeV]");
    hEkin_vs_Theta_gammas_lab_MC->GetYaxis()->SetTitle("#theta^{lab}_{#gamma} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_gammas_lab_MC,hpathMC);

    hEkin_vs_Theta_gammas_lab_acc_MC = new TH2F("hEkin_vs_Theta_gammas_lab_acc_MC","",500,0.,0.85,500,0.,180.);
    hEkin_vs_Theta_gammas_lab_acc_MC->GetXaxis()->SetTitle("E^{kin}_{#gamma} [GeV]");
    hEkin_vs_Theta_gammas_lab_acc_MC->GetYaxis()->SetTitle("#theta^{lab}_{#gamma} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_gammas_lab_acc_MC,hpathMC);

    ////
    hOpeningAngle_g1_g2_lab_MC=new TH1F("hOpeningAngle_g1_g2_lab_MC","",360,0.,180.);
    hOpeningAngle_g1_g2_lab_MC->GetXaxis()->SetTitle("#theta_{#gamma_{1},#gamma_{2}} [deg]");
    hOpeningAngle_g1_g2_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hOpeningAngle_g1_g2_lab_MC,hpathMC);

    hOpeningAngle_g1_g2_cm_MC=new TH1F("hOpeningAngle_g1_g2_cm_MC","",360,0.,180.);
    hOpeningAngle_g1_g2_cm_MC->GetXaxis()->SetTitle("#theta_{#gamma_{1},#gamma_{2}} [deg]");
    hOpeningAngle_g1_g2_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hOpeningAngle_g1_g2_cm_MC,hpathMC);

    hTheta_g1_vs_Theta_g2_lab_MC = new TH2F("hTheta_g1_vs_Theta_g2_lab_MC","",540,0.,180.,540,0.,180.);
    hTheta_g1_vs_Theta_g2_lab_MC->GetXaxis()->SetTitle("#theta^{lab}_{#gamma_{1}} [deg]");
    hTheta_g1_vs_Theta_g2_lab_MC->GetYaxis()->SetTitle("#theta^{lab}_{#gamma_{2}} [deg]");
    gHistoManager->Add(hTheta_g1_vs_Theta_g2_lab_MC,hpathMC);

    hTheta_g1_vs_Theta_g2_cm_MC = new TH2F("hTheta_g1_vs_Theta_g2_cm_MC","",540,0.,180.,540,0.,180.);
    hTheta_g1_vs_Theta_g2_cm_MC->GetXaxis()->SetTitle("#theta^{cm}_{#gamma_{1}} [deg]");
    hTheta_g1_vs_Theta_g2_cm_MC->GetYaxis()->SetTitle("#theta^{cm}_{#gamma_{2}} [deg]");
    gHistoManager->Add(hTheta_g1_vs_Theta_g2_cm_MC,hpathMC);

    //pion//
    hEkin_pi0_lab_MC = new TH1F("hEkin_pi0_lab_MC","",500,0.,0.75);
    hEkin_pi0_lab_MC->GetXaxis()->SetTitle("E^{kin}_{#pi^{0}} [GeV]");
    hEkin_pi0_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_pi0_lab_MC,hpathMC);

    hEkin_pi0_cm_MC = new TH1F("hEkin_pi0_cm_MC","",500,0.,0.5);
    hEkin_pi0_cm_MC->GetXaxis()->SetTitle("E^{kin}_{#pi^{0}} [GeV]");
    hEkin_pi0_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_pi0_cm_MC,hpathMC);

    hE_pi0_lab_MC = new TH1F("hE_pi0_lab_MC","",500,0.,1.);
    hE_pi0_lab_MC->GetXaxis()->SetTitle("E^{lab}_{#pi^{0}} [GeV]");
    hE_pi0_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_pi0_lab_MC,hpathMC);

    hE_pi0_cm_MC = new TH1F("hE_pi0_cm_MC","",500,0.,0.6);
    hE_pi0_cm_MC->GetXaxis()->SetTitle("E^{cm}_{#pi^{0}} [GeV]");
    hE_pi0_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_pi0_cm_MC,hpathMC);

    hp_pi0_lab_MC = new TH1F("hp_pi0_lab_MC","",500,0.,0.8);
    hp_pi0_lab_MC->GetXaxis()->SetTitle("p^{lab}_{#pi^{0}} [GeV/c]");
    hp_pi0_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_pi0_lab_MC,hpathMC);

    hp_pi0_cm_MC = new TH1F("hp_pi0_cm_MC","",500,0.,0.6);
    hp_pi0_cm_MC->GetXaxis()->SetTitle("p^{cm}_{#pi^{0}} [GeV/c]");
    hp_pi0_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_pi0_cm_MC,hpathMC);

    hTheta_pi0_lab_MC = new TH1F("hTheta_pi0_lab_MC","",360,0.,180.);
    hTheta_pi0_lab_MC->GetXaxis()->SetTitle("#theta^{lab}_{#pi^{0}} [deg]");
    hTheta_pi0_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_pi0_lab_MC,hpathMC);

    hTheta_pi0_cm_MC = new TH1F("hTheta_pi0_cm_MC","",360,0.,180.);
    hTheta_pi0_cm_MC->GetXaxis()->SetTitle("#theta^{cm}_{#pi^{0}} [deg]");
    hTheta_pi0_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_pi0_cm_MC,hpathMC);

    hPhi_pi0_lab_MC = new TH1F("hPhi_pi0_lab_MC","",360,-180.,180.);
    hPhi_pi0_lab_MC->GetXaxis()->SetTitle("#phi^{lab}_{#pi^{0}} [deg]");
    hPhi_pi0_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_pi0_lab_MC,hpathMC);

    hPhi_pi0_cm_MC = new TH1F("hPhi_pi0_cm_MC","",360,-180.,180.);
    hPhi_pi0_cm_MC->GetXaxis()->SetTitle("#phi^{cm}_{#pi^{0}} [deg]");
    hPhi_pi0_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_pi0_cm_MC,hpathMC);

    hEkin_vs_Theta_pi0_lab_MC = new TH2F("hEkin_vs_Theta_pi0_lab_MC","",500,0.,0.8,500,0.,180.);
    hEkin_vs_Theta_pi0_lab_MC->GetXaxis()->SetTitle("E^{kin}_{#pi^{0}} [GeV]");
    hEkin_vs_Theta_pi0_lab_MC->GetYaxis()->SetTitle("#theta^{lab}_{#pi^{0}} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_pi0_lab_MC,hpathMC);

    hEkin_vs_Theta_pi0_cm_MC = new TH2F("hEkin_vs_Theta_pi0_cm_MC","",500,0.,0.5,500,0.,180.);
    hEkin_vs_Theta_pi0_cm_MC->GetXaxis()->SetTitle("E^{kin,cm}_{#pi^{0}} [GeV]");
    hEkin_vs_Theta_pi0_cm_MC->GetYaxis()->SetTitle("#theta^{cm}_{#pi^{0}} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_pi0_cm_MC,hpathMC);

    ////
    hOpeningAngle_pi0_p_lab_MC = new TH1F("hOpeningAngle_pi0_p_lab_MC","",400,0.,200.);
    hOpeningAngle_pi0_p_lab_MC->GetXaxis()->SetTitle("#theta_{p,#pi^{0}} [deg]");
    hOpeningAngle_pi0_p_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hOpeningAngle_pi0_p_lab_MC,hpathMC);

    hOpeningAngle_pi0_p_cm_MC = new TH1F("hOpeningAngle_pi0_p_cm_MC","",400,0.,200.);
    hOpeningAngle_pi0_p_cm_MC->GetXaxis()->SetTitle("#theta_{p,#pi^{0}} [deg]");
    hOpeningAngle_pi0_p_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hOpeningAngle_pi0_p_cm_MC,hpathMC);

    ////
    hIM_pion_MC = new TH1F("hIM_pion_MC","",500,0.,0.4);
    hIM_pion_MC->GetXaxis()->SetTitle("m_{#gamma_{1}#gamma_{2}} [GeV]");
    hIM_pion_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hIM_pion_MC,hpathMC);

    hMM_deuteron_MC=new TH1F("hMM_deuteron_MC","",500,0.,2.5);
    hMM_deuteron_MC->GetXaxis()->SetTitle("m_{x} [GeV]");
    hMM_deuteron_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hMM_deuteron_MC,hpathMC);

    hMM_proton_MC=new TH1F("hMM_proton_MC","",500,0.,1.25);
    hMM_proton_MC->GetXaxis()->SetTitle("m_{x} [GeV]");
    hMM_proton_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hMM_proton_MC,hpathMC);

    hME_deuteron_MC=new TH1F("hME_deuteron_MC","",500,0.,2.5);
    hME_deuteron_MC->GetXaxis()->SetTitle("E_{x} [GeV]");
    hME_deuteron_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hME_deuteron_MC,hpathMC);

    hME_proton_MC=new TH1F("hME_proton_MC","",500,0.,1.25);
    hME_proton_MC->GetXaxis()->SetTitle("E_{x} [GeV]");
    hME_proton_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hME_proton_MC,hpathMC);

    hMMvsME_nucleon_MC = new TH2F("hMMvsME_nucleon_MC","",500,0.,2.5,500,0.,2.5);
    hMMvsME_nucleon_MC->GetXaxis()->SetTitle("m_{x} [GeV]");
    hMMvsME_nucleon_MC->GetYaxis()->SetTitle("E_{x} [GeV]");
    gHistoManager->Add(hMMvsME_nucleon_MC,hpathMC);

//////////////////////////////////////////RECONSTRUCTED/////////////////////////////////////////

    ////DATA level 0: trigger #10////

    hNeutralTracksCD[0][0] = new TH1F("hNeutralTracksCD_lev0","",11,-0.5,10.5);
    hNeutralTracksCD[0][0]->GetXaxis()->SetTitle("Neutral Tracks in CD");
    hNeutralTracksCD[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hNeutralTracksCD[0][0],hpath[0][0]);

    hChargedTracksCD[0][0] = new TH1F("hChargedTracksCD_lev0","",11,-0.5,10.5);
    hChargedTracksCD[0][0]->GetXaxis()->SetTitle("Charged Tracks in CD");
    hChargedTracksCD[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hChargedTracksCD[0][0],hpath[0][0]);

    hChargedTracksFD[0][0] = new TH1F("hChargedTracksFD_lev0","",11,-0.5,10.5);
    hChargedTracksFD[0][0]->GetXaxis()->SetTitle("Charged Tracks in FD");
    hChargedTracksFD[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hChargedTracksFD[0][0],hpath[0][0]);

    hp_beam[0][0] = new TH1F("hp_beam_lev0","",500,1.41,1.65);
    hp_beam[0][0]->GetXaxis()->SetTitle("p_{beam} [GeV/c]");
    hp_beam[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_beam[0][0],hpath[0][0]);

    hQ[0][0][0][0] = new TH1F("hQ_lev0","",40,-70.,30.);
    hQ[0][0][0][0]->GetXaxis()->SetTitle("Q [MeV]");
    hQ[0][0][0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hQ[0][0][0][0],hpath[0][0]);

    //FD
    hEdepFWC1vsFRH1[0][0] = new TH2F("hEdepFWC1vsFRH1_lev0","",1000,0.,0.02,1000,0.,0.5);
    hEdepFWC1vsFRH1[0][0]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
    hEdepFWC1vsFRH1[0][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
    gHistoManager->Add(hEdepFWC1vsFRH1[0][0],hpath[0][0]);

    hEdepFWC2vsFRH1[0][0] = new TH2F("hEdepFWC2vsFRH1_lev0","",1000,0.,0.02,1000,0.,0.5);
    hEdepFWC2vsFRH1[0][0]->GetXaxis()->SetTitle("Edep(FWC2) [GeV]");
    hEdepFWC2vsFRH1[0][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
    gHistoManager->Add(hEdepFWC2vsFRH1[0][0],hpath[0][0]);

    hEdepFTH1vsFRH1[0][0] = new TH2F("hEdepFTH1vsFRH1_lev0","",1000,0.,0.05,1000,0.,0.5);
    hEdepFTH1vsFRH1[0][0]->GetXaxis()->SetTitle("Edep(FTH1) [GeV]");
    hEdepFTH1vsFRH1[0][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
    gHistoManager->Add(hEdepFTH1vsFRH1[0][0],hpath[0][0]);

    hEdepFRH1vsFRH2[0][0] = new TH2F("hEdepFRH1vsFRH2_lev0","",1000,0.,0.5,1000,0.,0.5);
    hEdepFRH1vsFRH2[0][0]->GetXaxis()->SetTitle("Edep(FRH1) [GeV]");
    hEdepFRH1vsFRH2[0][0]->GetYaxis()->SetTitle("Edep(FRH2) [GeV]");
    gHistoManager->Add(hEdepFRH1vsFRH2[0][0],hpath[0][0]);

    hEdepFRH2vsFRH3[0][0] = new TH2F("hEdepFRH2vsFRH3_lev0","",1000,0.,0.5,1000,0.,0.5);
    hEdepFRH2vsFRH3[0][0]->GetXaxis()->SetTitle("Edep(FRH2) [GeV]");
    hEdepFRH2vsFRH3[0][0]->GetYaxis()->SetTitle("Edep(FRH3) [GeV]");
    gHistoManager->Add(hEdepFRH2vsFRH3[0][0],hpath[0][0]);

    hEdepFWC1vsFRH1FRH2FRH3[0][0] = new TH2F("hEdepFWC1vsFRH1FRH2FRH3_lev0","",1000,0.,0.02,1000,0.,1.);
    hEdepFWC1vsFRH1FRH2FRH3[0][0]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
    hEdepFWC1vsFRH1FRH2FRH3[0][0]->GetYaxis()->SetTitle("Edep(FRH1+FRH2+FRH3) [GeV]");
    gHistoManager->Add(hEdepFWC1vsFRH1FRH2FRH3[0][0],hpath[0][0]);

    hTheta_FDC[0][0] = new TH1F("hTheta_FDC_lev0","",250,0.,25.);
    hTheta_FDC[0][0]->GetXaxis()->SetTitle("#theta^{lab}_{FD} [deg]");
    hTheta_FDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_FDC[0][0],hpath[0][0]);

    hPhi_FDC[0][0] = new TH1F("hPhi_FDC_lev0","",360,-180.,180.);
    hPhi_FDC[0][0]->GetXaxis()->SetTitle("#phi^{lab}_{FD} [deg]");
    hPhi_FDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_FDC[0][0],hpath[0][0]);

    hTime_FDC[0][0] = new TH1F("hTime_FDC_lev0","",1000,0.,2500.);
    hTime_FDC[0][0]->GetXaxis()->SetTitle("Time_{FD} [ns]");
    hTime_FDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTime_FDC[0][0],hpath[0][0]);

    //charged in CD
    hEdepPSBvsSEC[0][0] = new TH2F("hEdepPSBvsSEC_lev0","",1000,0.,0.75,1000,0.,0.025);
    hEdepPSBvsSEC[0][0]->GetXaxis()->SetTitle("E_{dep(SEC)} [GeV]");
    hEdepPSBvsSEC[0][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
    gHistoManager->Add(hEdepPSBvsSEC[0][0],hpath[0][0]);

    hEdepPSBvsSigMom[0][0] = new TH2F("hEdepPSBvsSigMom_lev0","",1000,-2.5,2.5,1000,0.,0.025);
    hEdepPSBvsSigMom[0][0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
    hEdepPSBvsSigMom[0][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
    gHistoManager->Add(hEdepPSBvsSigMom[0][0],hpath[0][0]);

    hEdepSECvsSigMom[0][0] = new TH2F("hEdepSECvsSigMom_lev0","",1000,-2.5,2.5,1000,0.,0.75);
    hEdepSECvsSigMom[0][0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
    hEdepSECvsSigMom[0][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
    gHistoManager->Add(hEdepSECvsSigMom[0][0],hpath[0][0]);

    hMom_CDC[0][0] = new TH1F("hMom_CDC_lev0","",750,0.,2.5);
    hMom_CDC[0][0]->GetXaxis()->SetTitle("p_{CD} [GeV/c]");
    hMom_CDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hMom_CDC[0][0],hpath[0][0]);

    hTheta_CDC[0][0] = new TH1F("hTheta_CDC_lev0","",360,0.,180.);
    hTheta_CDC[0][0]->GetXaxis()->SetTitle("#theta_{CD} [deg]");
    hTheta_CDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_CDC[0][0],hpath[0][0]);

    hPhi_CDC[0][0] = new TH1F("hPhi_CDC_lev0","",360,-180.,180.);
    hPhi_CDC[0][0]->GetXaxis()->SetTitle("#phi_{CD} [deg]");
    hPhi_CDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_CDC[0][0],hpath[0][0]);

    hTime_CDC[0][0] = new TH1F("hTime_CDC_lev0","",1000,0.,2500.);
    hTime_CDC[0][0]->GetXaxis()->SetTitle("Time_{CD} [ns]");
    hTime_CDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTime_CDC[0][0],hpath[0][0]);

    //neutral in CD
    hMom_CDN[0][0] = new TH1F("hMom_CDN_lev0","",500,0.,0.8);
    hMom_CDN[0][0]->GetXaxis()->SetTitle("p_{CD} [GeV/c]");
    hMom_CDN[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hMom_CDN[0][0],hpath[0][0]);

    hTheta_CDN[0][0]= new TH1F("hTheta_CDN_lev0","",360,0.,180.);
    hTheta_CDN[0][0]->GetXaxis()->SetTitle("#theta_{CD} [deg]");
    hTheta_CDN[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_CDN[0][0],hpath[0][0]);

    hPhi_CDN[0][0] = new TH1F("hPhi_CDN_lev0","",360,-180.,180.);
    hPhi_CDN[0][0]->GetXaxis()->SetTitle("#phi_{CD} [deg]");
    hPhi_CDN[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_CDN[0][0],hpath[0][0]);

    hTime_CDN[0][0] = new TH1F("hTime_CDN_lev0","",1000,0.,2500.);
    hTime_CDN[0][0]->GetXaxis()->SetTitle("Time_{CD} [ns]");
    hTime_CDN[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTime_CDN[0][0],hpath[0][0]);

    ////DATA level 1: level 0 + 1chFD1chCD2ntCD////

    hNeutralTracksCD[1][0] = new TH1F("hNeutralTracksCD_lev1","",11,-0.5,10.5);
    hNeutralTracksCD[1][0]->GetXaxis()->SetTitle("Neutral Tracks in CD");
    hNeutralTracksCD[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hNeutralTracksCD[1][0],hpath[1][0]);

    hChargedTracksCD[1][0] = new TH1F("hChargedTracksCD_lev1","",11,-0.5,10.5);
    hChargedTracksCD[1][0]->GetXaxis()->SetTitle("Charged Tracks in CD");
    hChargedTracksCD[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hChargedTracksCD[1][0],hpath[1][0]);

    hChargedTracksFD[1][0] = new TH1F("hChargedTracksFD_lev1","",11,-0.5,10.5);
    hChargedTracksFD[1][0]->GetXaxis()->SetTitle("Charged Tracks in FD");
    hChargedTracksFD[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hChargedTracksFD[1][0],hpath[1][0]);

    hQ[1][0][0][0] = new TH1F("hQ_lev1","",40,-70.,30.);
    hQ[1][0][0][0]->GetXaxis()->SetTitle("Q [MeV]");
    hQ[1][0][0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hQ[1][0][0][0],hpath[1][0]);

    //FD
    hEdepFWC1vsFRH1[1][0] = new TH2F("hEdepFWC1vsFRH1_lev1","",1000,0.,0.02,1000,0.,0.5);
    hEdepFWC1vsFRH1[1][0]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
    hEdepFWC1vsFRH1[1][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
    gHistoManager->Add(hEdepFWC1vsFRH1[1][0],hpath[1][0]);

    hEdepFWC2vsFRH1[1][0] = new TH2F("hEdepFWC2vsFRH1_lev1","",1000,0.,0.02,1000,0.,0.5);
    hEdepFWC2vsFRH1[1][0]->GetXaxis()->SetTitle("Edep(FWC2) [GeV]");
    hEdepFWC2vsFRH1[1][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
    gHistoManager->Add(hEdepFWC2vsFRH1[1][0],hpath[1][0]);

    hEdepFTH1vsFRH1[1][0] = new TH2F("hEdepFTH1vsFRH1_lev1","",1000,0.,0.05,1000,0.,0.5);
    hEdepFTH1vsFRH1[1][0]->GetXaxis()->SetTitle("Edep(FTH1) [GeV]");
    hEdepFTH1vsFRH1[1][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
    gHistoManager->Add(hEdepFTH1vsFRH1[1][0],hpath[1][0]);

    hEdepFRH1vsFRH2[1][0] = new TH2F("hEdepFRH1vsFRH2_lev1","",1000,0.,0.5,1000,0.,0.5);
    hEdepFRH1vsFRH2[1][0]->GetXaxis()->SetTitle("Edep(FRH1) [GeV]");
    hEdepFRH1vsFRH2[1][0]->GetYaxis()->SetTitle("Edep(FRH2) [GeV]");
    gHistoManager->Add(hEdepFRH1vsFRH2[1][0],hpath[1][0]);

    hEdepFRH2vsFRH3[1][0] = new TH2F("hEdepFRH2vsFRH3_lev1","",1000,0.,0.5,1000,0.,0.5);
    hEdepFRH2vsFRH3[1][0]->GetXaxis()->SetTitle("Edep(FRH2) [GeV]");
    hEdepFRH2vsFRH3[1][0]->GetYaxis()->SetTitle("Edep(FRH3) [GeV]");
    gHistoManager->Add(hEdepFRH2vsFRH3[1][0],hpath[1][0]);

    hEdepFWC1vsFRH1FRH2FRH3[1][0] = new TH2F("hEdepFWC1vsFRH1FRH2FRH3_lev1","",1000,0.,0.02,1000,0.,1.);
    hEdepFWC1vsFRH1FRH2FRH3[1][0]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
    hEdepFWC1vsFRH1FRH2FRH3[1][0]->GetYaxis()->SetTitle("Edep(FRH1+FRH2+FRH3) [GeV]");
    gHistoManager->Add(hEdepFWC1vsFRH1FRH2FRH3[1][0],hpath[1][0]);

    hTheta_FDC[1][0] = new TH1F("hTheta_FDC_lev1","",250,0.,25.);
    hTheta_FDC[1][0]->GetXaxis()->SetTitle("#theta^{lab}_{FD} [deg]");
    hTheta_FDC[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_FDC[1][0],hpath[1][0]);

    hPhi_FDC[1][0] = new TH1F("hPhi_FDC_lev1","",360,-180.,180.);
    hPhi_FDC[1][0]->GetXaxis()->SetTitle("#phi^{lab}_{FD} [deg]");
    hPhi_FDC[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_FDC[1][0],hpath[1][0]);

    hTime_FDC[1][0] = new TH1F("hTime_FDC_lev1","",1000,0.,2500.);
    hTime_FDC[1][0]->GetXaxis()->SetTitle("Time_{FD} [ns]");
    hTime_FDC[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTime_FDC[1][0],hpath[1][0]);

    //charged in CD
    hEdepPSBvsSEC[1][0] = new TH2F("hEdepPSBvsSEC_lev1","",1000,0.,0.6,1000,0.,0.015);
    hEdepPSBvsSEC[1][0]->GetXaxis()->SetTitle("E_{dep(SEC)} [GeV]");
    hEdepPSBvsSEC[1][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
    gHistoManager->Add(hEdepPSBvsSEC[1][0],hpath[1][0]);

    hEdepPSBvsSigMom[1][0] = new TH2F("hEdepPSBvsSigMom_lev1","",1000,-2.,2.,1000,0.,0.020);
    hEdepPSBvsSigMom[1][0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
    hEdepPSBvsSigMom[1][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
    gHistoManager->Add(hEdepPSBvsSigMom[1][0],hpath[1][0]);

    hEdepSECvsSigMom[1][0] = new TH2F("hEdepSECvsSigMom_lev1","",1000,-2.,2.,1000,0.,0.6);
    hEdepSECvsSigMom[1][0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
    hEdepSECvsSigMom[1][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
    gHistoManager->Add(hEdepSECvsSigMom[1][0],hpath[1][0]);

    hMom_CDC[1][0] = new TH1F("hMom_CDC_lev1","",750,0.,2.5);
    hMom_CDC[1][0]->GetXaxis()->SetTitle("p_{CD} [GeV/c]");
    hMom_CDC[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hMom_CDC[1][0],hpath[1][0]);

    hTheta_CDC[1][0] = new TH1F("hTheta_CDC_lev1","",360,0.,180.);
    hTheta_CDC[1][0]->GetXaxis()->SetTitle("#theta_{CD} [deg]");
    hTheta_CDC[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_CDC[1][0],hpath[1][0]);

    hPhi_CDC[1][0] = new TH1F("hPhi_CDC_lev1","",360,-180.,180.);
    hPhi_CDC[1][0]->GetXaxis()->SetTitle("#phi_{CD} [deg]");
    hPhi_CDC[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_CDC[1][0],hpath[1][0]);

    hTime_CDC[1][0] = new TH1F("hTime_CDC_lev1","",1000,0.,2500.);
    hTime_CDC[1][0]->GetXaxis()->SetTitle("Time_{CD} [ns]");
    hTime_CDC[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTime_CDC[1][0],hpath[1][0]);

    //gamma #1//
    hp_g1_lab[1][0] = new TH1F("hp_g1_lab_lev1","",500,0.,0.8);
    hp_g1_lab[1][0]->GetXaxis()->SetTitle("p^{lab}_{#gamma_{1}} [GeV/c]");
    hp_g1_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_g1_lab[1][0],hpathG[1][0]);

    hTheta_g1_lab[1][0] = new TH1F("hTheta_g1_lab_lev1","",360,0.,180.);
    hTheta_g1_lab[1][0]->GetXaxis()->SetTitle("#theta^{lab}_{#gamma_{1}} [deg]");
    hTheta_g1_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_g1_lab[1][0],hpathG[1][0]);

    hPhi_g1_lab[1][0] = new TH1F("hPhi_g1_lab_lev1","",360,-180.,180.);
    hPhi_g1_lab[1][0]->GetXaxis()->SetTitle("#phi^{lab}_{#gamma_{1}} [deg]");
    hPhi_g1_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_g1_lab[1][0],hpathG[1][0]);

    hEkin_vs_Theta_g1_lab[1][0] = new TH2F("hEkin_vs_Theta_g1_lab_lev1","",500,0.,0.8,500,0.,180.);
    hEkin_vs_Theta_g1_lab[1][0]->GetXaxis()->SetTitle("E^{kin}_{#gamma_{1}} [GeV]");
    hEkin_vs_Theta_g1_lab[1][0]->GetYaxis()->SetTitle("#theta^{lab}_{#gamma_{1}} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_g1_lab[1][0],hpathG[1][0]);

    //gamma#2//
    hp_g2_lab[1][0] = new TH1F("hp_g2_lab_lev1","",500,0.,0.8);
    hp_g2_lab[1][0]->GetXaxis()->SetTitle("p^{lab}_{#gamma_{2}} [GeV/c]");
    hp_g2_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_g2_lab[1][0],hpathG[1][0]);

    hTheta_g2_lab[1][0] = new TH1F("hTheta_g2_lab_lev1","",360,0.,180.);
    hTheta_g2_lab[1][0]->GetXaxis()->SetTitle("#theta^{lab}_{#gamma_{2}} [deg]");
    hTheta_g2_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_g2_lab[1][0],hpathG[1][0]);

    hPhi_g2_lab[1][0] = new TH1F("hPhi_g2_lab_lev1","",360,-180.,180.);
    hPhi_g2_lab[1][0]->GetXaxis()->SetTitle("#phi^{lab}_{#gamma_{2}} [deg]");
    hPhi_g2_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_g2_lab[1][0],hpathG[1][0]);

    hEkin_vs_Theta_g2_lab[1][0] = new TH2F("hEkin_vs_Theta_g2_lab_lev1","",500,0.,0.8,500,0.,180.);
    hEkin_vs_Theta_g2_lab[1][0]->GetXaxis()->SetTitle("E^{kin}_{#gamma_{2}} [GeV]");
    hEkin_vs_Theta_g2_lab[1][0]->GetYaxis()->SetTitle("#theta^{lab}_{#gamma_{2}} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_g2_lab[1][0],hpathG[1][0]);
    //
    hOpeningAngle_g1_g2_lab[1][0] = new TH1F("hOpeningAngle_g1_g2_lab_lev1","",360,0.,180.);
    hOpeningAngle_g1_g2_lab[1][0]->GetXaxis()->SetTitle("#theta^{lab}_{#gamma,#gamma} [deg]");
    hOpeningAngle_g1_g2_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hOpeningAngle_g1_g2_lab[1][0],hpathG[1][0]);

    hTheta_g1_vs_Theta_g2_lab[1][0] = new TH2F("hTheta_g1_vs_Theta_g2_lab_lev1","",540,0.,180.,540,0.,180.);
    hTheta_g1_vs_Theta_g2_lab[1][0]->GetXaxis()->SetTitle("#theta^{lab}_{#gamma_{1}} [deg]");
    hTheta_g1_vs_Theta_g2_lab[1][0]->GetYaxis()->SetTitle("#theta^{lab}_{#gamma_{2}} [deg]");
    gHistoManager->Add(hTheta_g1_vs_Theta_g2_lab[1][0],hpathG[1][0]);

    //pion//
    hEkin_pi0_lab[1][0] = new TH1F("hEkin_pi0_lab_lev1","",500,0.,0.75);
    hEkin_pi0_lab[1][0]->GetXaxis()->SetTitle("E^{kin}_{#pi^{0}} [GeV]");
    hEkin_pi0_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_pi0_lab[1][0],hpathPi[1][0]);

    hE_pi0_lab[1][0] = new TH1F("hE_pi0_lab_lev1","",500,0.,1.);
    hE_pi0_lab[1][0]->GetXaxis()->SetTitle("E^{lab}_{#pi^{0}} [GeV]");
    hE_pi0_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_pi0_lab[1][0],hpathPi[1][0]);

    hp_pi0_lab[1][0] = new TH1F("hp_pi0_lab_lev1","",500,0.,0.8);
    hp_pi0_lab[1][0]->GetXaxis()->SetTitle("p^{lab}_{#pi^{0}} [GeV/c]");
    hp_pi0_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_pi0_lab[1][0],hpathPi[1][0]);

    hTheta_pi0_lab[1][0] = new TH1F("hTheta_pi0_lab_lev1","",360,0.,180.);
    hTheta_pi0_lab[1][0]->GetXaxis()->SetTitle("#theta^{lab}_{#pi^{0}} [deg]");
    hTheta_pi0_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_pi0_lab[1][0],hpathPi[1][0]);

    hPhi_pi0_lab[1][0] = new TH1F("hPhi_pi0_lab_lev1","",360,-180.,180.);
    hPhi_pi0_lab[1][0]->GetXaxis()->SetTitle("#phi^{lab}_{#pi^{0}} [deg]");
    hPhi_pi0_lab[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_pi0_lab[1][0],hpathPi[1][0]);

    hEkin_vs_Theta_pi0_lab[1][0] = new TH2F("hEkin_vs_Theta_pi0_lab_lev1","",500,0.,0.8,500,0.,180.);
    hEkin_vs_Theta_pi0_lab[1][0]->GetXaxis()->SetTitle("E^{kin}_{#pi^{0}} [GeV]");
    hEkin_vs_Theta_pi0_lab[1][0]->GetYaxis()->SetTitle("#theta^{lab}_{#pi^{0}} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_pi0_lab[1][0],hpathPi[1][0]);

    ////
    hEnergy_additional_gammas = new TH1F("hEnergy_additional_gammas","",400,0.,0.2);
    hEnergy_additional_gammas->GetXaxis()->SetTitle("p_{#gamma} [GeV/c]");
    hEnergy_additional_gammas->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEnergy_additional_gammas,hpath[1][0]);


    hEnergy_additional_gammas_cut = new TH1F("hEnergy_additional_gammas_cut","",400,0.,0.2);
    hEnergy_additional_gammas_cut->GetXaxis()->SetTitle("p_{#gamma} [GeV/c]");
    hEnergy_additional_gammas_cut->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEnergy_additional_gammas_cut,hpath[1][0]);

    ////DATA level 2: level 1 + CUTG1[PSBvsSEC] | DATA level 3: level 2 + positively charged in CD////

    for (Int_t l = 2; l < 4; l++) {
        for (Int_t k = 0; k < 7; k++) {

            hQ[l][k][0][0] = new TH1F(Form("hQ_lev%d_cut%d",l,k),"",40,-70.,30.);
            hQ[l][k][0][0]->GetXaxis()->SetTitle("Q [MeV]");
            hQ[l][k][0][0]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hQ[l][k][0][0],hpath[l][k]);

            //FD
            hEdepFWC1vsFRH1[l][k] =new TH2F(Form("hEdepFWC1vsFRH1_lev%d_cut%d",l,k),"",1000,0.,0.02,1000,0.,0.5);
            hEdepFWC1vsFRH1[l][k]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
            hEdepFWC1vsFRH1[l][k]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
            gHistoManager->Add(hEdepFWC1vsFRH1[l][k],hpath[l][k]);

            hEdepFWC2vsFRH1[l][k] =new TH2F(Form("hEdepFWC2vsFRH1_lev%d_cut%d",l,k),"",1000,0.,0.02,1000,0.,0.5);
            hEdepFWC2vsFRH1[l][k]->GetXaxis()->SetTitle("Edep(FWC2) [GeV]");
            hEdepFWC2vsFRH1[l][k]->GetYaxis()->SetTitle("Edep(FRH1) [Ge0V]");
            gHistoManager->Add(hEdepFWC2vsFRH1[l][k],hpath[l][k]);

            hEdepFTH1vsFRH1[l][k] =new TH2F(Form("hEdepFTH1vsFRH1_lev%d_cut%d",l,k),"",1000,0.,0.05,1000,0.,0.5);
            hEdepFTH1vsFRH1[l][k]->GetXaxis()->SetTitle("Edep(FTH1) [GeV]");
            hEdepFTH1vsFRH1[l][k]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
            gHistoManager->Add(hEdepFTH1vsFRH1[l][k],hpath[l][k]);

            hEdepFRH1vsFRH2[l][k] =new TH2F(Form("hEdepFRH1vsFRH2_lev%d_cut%d",l,k),"",1000,0.,0.5,1000,0.,0.5);
            hEdepFRH1vsFRH2[l][k]->GetXaxis()->SetTitle("Edep(FRH1) [GeV]");
            hEdepFRH1vsFRH2[l][k]->GetYaxis()->SetTitle("Edep(FRH2) [GeV]");
            gHistoManager->Add(hEdepFRH1vsFRH2[l][k],hpath[l][k]);

            hEdepFRH2vsFRH3[l][k] =new TH2F(Form("hEdepFRH2vsFRH3_lev%d_cut%d",l,k),"",1000,0.,0.5,1000,0.,0.5);
            hEdepFRH2vsFRH3[l][k]->GetXaxis()->SetTitle("Edep(FRH2) [GeV]");
            hEdepFRH2vsFRH3[l][k]->GetYaxis()->SetTitle("Edep(FRH3) [GeV]");
            gHistoManager->Add(hEdepFRH2vsFRH3[l][k],hpath[l][k]);

            hEdepFWC1vsFRH1FRH2FRH3[l][k] =new TH2F(Form("hEdepFWC1vsFRH1FRH2FRH3_lev%d_cut%d",l,k),"",1000,0.,0.02,1000,0.,1.);
            hEdepFWC1vsFRH1FRH2FRH3[l][k]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
            hEdepFWC1vsFRH1FRH2FRH3[l][k]->GetYaxis()->SetTitle("Edep(FRH1+FRH2+FRH3) [GeV]");
            gHistoManager->Add(hEdepFWC1vsFRH1FRH2FRH3[l][k],hpath[l][k]);

            //CD
            hEdepPSBvsSEC[l][k] = new TH2F(Form("hEdepPSBvsSEC_lev%d_cut%d",l,k),"",1000,0.,0.6,1000,0.,0.015);
            hEdepPSBvsSEC[l][k]->GetXaxis()->SetTitle("E_{dep(SEC)} [GeV]");
            hEdepPSBvsSEC[l][k]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
            gHistoManager->Add(hEdepPSBvsSEC[l][k],hpath[l][k]);

            hEdepPSBvsSigMom[l][k] = new TH2F(Form("hEdepPSBvsSigMom_lev%d_cut%d",l,k),"",1000,-2.,2.,1000,0.,0.02);
            hEdepPSBvsSigMom[l][k]->GetXaxis()->SetTitle("Momentum [GeV/c]");
            hEdepPSBvsSigMom[l][k]->GetYaxis()->SetTitle("E_{dep(SEC)} [GeV]");
            gHistoManager->Add(hEdepPSBvsSigMom[l][k],hpath[l][k]);

            hEdepSECvsSigMom[l][k] = new TH2F(Form("hEdepSECvsSigMom_lev%d_cut%d",l,k),"",1000,-2.,2.,1000,0.,0.6);
            hEdepSECvsSigMom[l][k]->GetXaxis()->SetTitle("Momentum [GeV/c]");
            hEdepSECvsSigMom[l][k]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
            gHistoManager->Add(hEdepSECvsSigMom[l][k],hpath[l][k]);

            ////
            hOpeningAngle_pi0_p_lab[l][k] = new TH1F(Form("hOpeningAngle_pi0_p_lab_lev%d_cut%d",l,k),"",1000,0.,200.);
            hOpeningAngle_pi0_p_lab[l][k]->GetXaxis()->SetTitle("#theta^{lab}_{p,#pi^{0}} [deg]");
            hOpeningAngle_pi0_p_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hOpeningAngle_pi0_p_lab[l][k],hpath[l][k]);

            hOpeningAngle_pi0_p_cm[l][k] = new TH1F(Form("hOpeningAngle_pi0_p_cm_lev%d_cut%d",l,k),"",1000,0.,200.);
            hOpeningAngle_pi0_p_cm[l][k]->GetXaxis()->SetTitle("#theta^{cm}_{p,#pi^{0}} [deg]");
            hOpeningAngle_pi0_p_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hOpeningAngle_pi0_p_cm[l][k],hpath[l][k]);

            ////
            hTheta_p_vs_Theta_d_lab[l][k] = new TH2F(Form("hTheta_p_vs_Theta_d_lab_lev%d_cut%d",l,k),"",500,0.,20.,500,0.,180.);
            hTheta_p_vs_Theta_d_lab[l][k]->GetXaxis()->SetTitle("#theta^{lab}_{d} [deg]");
            hTheta_p_vs_Theta_d_lab[l][k]->GetYaxis()->SetTitle("#theta^{lab}_{p} [deg]");
            gHistoManager->Add(hTheta_p_vs_Theta_d_lab[l][k],hpath[l][k]);

            hTheta_p_vs_Theta_d_cm[l][k] = new TH2F(Form("hTheta_p_vs_Theta_d_cm_lev%d_cut%d",l,k),"",360,0.,180.,360,0.,180.);
            hTheta_p_vs_Theta_d_cm[l][k]->GetXaxis()->SetTitle("#theta^{cm}_{d} [deg]");
            hTheta_p_vs_Theta_d_cm[l][k]->GetYaxis()->SetTitle("#theta^{cm}_{p} [deg]");
            gHistoManager->Add(hTheta_p_vs_Theta_d_cm[l][k],hpath[l][k]);

            ////
            hIM_pion[l][k] = new TH1F(Form("hIM_pion_lev%d_cut%d",l,k),"",1000,0.,0.4);
            hIM_pion[l][k]->GetXaxis()->SetTitle("m_{#gamma_{1}#gamma_{2}} [GeV]");
            hIM_pion[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hIM_pion[l][k],hpath[l][k]);

            hMM_nucleon[l][k] = new TH1F(Form("hMM_nucleon_lev%d_cut%d",l,k),"",1000,0.,2.5);
            hMM_nucleon[l][k]->GetXaxis()->SetTitle("m_{x} [GeV]");
            hMM_nucleon[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hMM_nucleon[l][k],hpath[l][k]);

            hME_nucleon[l][k] = new TH1F(Form("hME_nucleon_lev%d_cut%d",l,k),"",1000,0.,2.5);
            hME_nucleon[l][k]->GetXaxis()->SetTitle("E_{x} [GeV]");
            hME_nucleon[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hME_nucleon[l][k],hpath[l][k]);

            hMMvsME_nucleon[l][k] = new TH2F(Form("hMMvsME_nucleon_lev%d_cut%d",l,k),"",500,0.,2.5,500,0.,3.);
            hMMvsME_nucleon[l][k]->GetXaxis()->SetTitle("m_{x} [GeV]");
            hMMvsME_nucleon[l][k]->GetYaxis()->SetTitle("E_{x} [GeV]");
            gHistoManager->Add(hMMvsME_nucleon[l][k],hpath[l][k]);

            hIMvsMM[l][k]  = new TH2F(Form("hIMvsMM_lev%d_cut%d",l,k),"",500,0.,0.5,500,0.,2.5);
            hIMvsMM[l][k]->GetXaxis()->SetTitle("m_{#gamma_{1}#gamma_{2}} [GeV]");
            hIMvsMM[l][k]->GetYaxis()->SetTitle("m_{x} [GeV]");
            gHistoManager->Add(hIMvsMM[l][k],hpath[l][k]);

            ////
            //proton//
            hEkin_p_lab[l][k] = new TH1F(Form("hEkin_p_lab_lev%d_cut%d",l,k),"",500,0.,1.0);
            hEkin_p_lab[l][k]->GetXaxis()->SetTitle("E^{kin}_{p} [GeV]");
            hEkin_p_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hEkin_p_lab[l][k],hpathP[l][k]);

            hEkin_p_cm[l][k] = new TH1F(Form("hEkin_p_cm_lev%d_cut%d",l,k),"",500,0.,0.5);
            hEkin_p_cm[l][k]->GetXaxis()->SetTitle("E^{kin,cm}_{p} [GeV]");
            hEkin_p_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hEkin_p_cm[l][k],hpathP[l][k]);

            hE_p_lab[l][k] = new TH1F(Form("hE_p_lab_lev%d_cut%d",l,k),"",500,0.9,1.5);
            hE_p_lab[l][k]->GetXaxis()->SetTitle("E^{lab}_{p} [GeV]");
            hE_p_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hE_p_lab[l][k],hpathP[l][k]);

            hE_p_cm[l][k] = new TH1F(Form("hE_p_cm_lev%d_cut%d",l,k),"",500,0.9,1.25);
            hE_p_cm[l][k]->GetXaxis()->SetTitle("E^{cm}_{p} [GeV]");
            hE_p_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hE_p_cm[l][k],hpathP[l][k]);

            hp_p_lab[l][k] = new TH1F(Form("hp_p_lab_lev%d_cut%d",l,k),"",500,0.,1.5);
            hp_p_lab[l][k]->GetXaxis()->SetTitle("p^{lab}_{p} [GeV/c]");
            hp_p_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hp_p_lab[l][k],hpathP[l][k]);

            hp_p_cm[l][k] = new TH1F(Form("hp_p_cm_lev%d_cut%d",l,k),"",500,0.,1.);
            hp_p_cm[l][k]->GetXaxis()->SetTitle("p^{cm}_{p} [GeV/c]");
            hp_p_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hp_p_cm[l][k],hpathP[l][k]);

            hTheta_p_lab[l][k] = new TH1F(Form("hTheta_p_lab_lev%d_cut%d",l,k),"",360,0.,180.);
            hTheta_p_lab[l][k]->GetXaxis()->SetTitle("#theta^{lab}_{p} [deg]");
            hTheta_p_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_p_lab[l][k],hpathP[l][k]);

            hTheta_p_cm[l][k] = new TH1F(Form("hTheta_p_cm_lev%d_cut%d",l,k),"",360,0.,180.);
            hTheta_p_cm[l][k]->GetXaxis()->SetTitle("#theta^{cm}_{p} [deg]");
            hTheta_p_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_p_cm[l][k],hpathP[l][k]);

            hPhi_p_lab[l][k] = new TH1F(Form("hPhi_p_lab_lev%d_cut%d",l,k),"",360,-180.,180.);
            hPhi_p_lab[l][k]->GetXaxis()->SetTitle("#phi^{lab}_{p} [deg]");
            hPhi_p_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hPhi_p_lab[l][k],hpathP[l][k]);

            hPhi_p_cm[l][k] = new TH1F(Form("hPhi_p_cm_lev%d_cut%d",l,k),"",360,-180.,180.);
            hPhi_p_cm[l][k]->GetXaxis()->SetTitle("#phi^{cm}_{p} [deg]");
            hPhi_p_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hPhi_p_cm[l][k],hpathP[l][k]);

            hEkin_vs_Theta_p_lab[l][k] = new TH2F(Form("hEkin_vs_Theta_p_lab_lev%d_cut%d",l,k),"",500,0.,1.,500,0.,180.);
            hEkin_vs_Theta_p_lab[l][k]->GetXaxis()->SetTitle("E^{kin}_{p} [GeV]");
            hEkin_vs_Theta_p_lab[l][k]->GetYaxis()->SetTitle("#theta^{lab}_{p} [deg]");
            gHistoManager->Add(hEkin_vs_Theta_p_lab[l][k],hpathP[l][k]);

            hEkin_vs_Theta_p_cm[l][k] = new TH2F(Form("hEkin_vs_Theta_p_cm_lev%d_cut%d",l,k),"",500,0.,0.5,500,0.,180.);
            hEkin_vs_Theta_p_cm[l][k]->GetXaxis()->SetTitle("E^{kin,cm}_{p} [GeV]");
            hEkin_vs_Theta_p_cm[l][k]->GetYaxis()->SetTitle("#theta^{cm}_{p} [deg]");
            gHistoManager->Add(hEkin_vs_Theta_p_cm[l][k],hpathP[l][k]);

            //deuteron//
            hEkin_d_lab[l][k] = new TH1F(Form("hEkin_d_lab_lev%d_cut%d",l,k),"",500,0.,1.0);
            hEkin_d_lab[l][k]->GetXaxis()->SetTitle("E^{kin}_{d} [GeV]");
            hEkin_d_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hEkin_d_lab[l][k],hpathD[l][k]);

            hEkin_d_cm[l][k] = new TH1F(Form("hEkin_d_cm_lev%d_cut%d",l,k),"",500,0.,0.1);
            hEkin_d_cm[l][k]->GetXaxis()->SetTitle("E^{kin,cm}_{d} [GeV]");
            hEkin_d_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hEkin_d_cm[l][k],hpathD[l][k]);

            hE_d_lab[l][k] = new TH1F(Form("hE_d_lab_lev%d_cut%d",l,k),"",500,1.8,2.6);
            hE_d_lab[l][k]->GetXaxis()->SetTitle("E^{lab}_{d} [GeV]");
            hE_d_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hE_d_lab[l][k],hpathD[l][k]);

            hE_d_cm[l][k] = new TH1F(Form("hE_d_cm_lev%d_cut%d",l,k),"",500,1.8,2.);
            hE_d_cm[l][k]->GetXaxis()->SetTitle("E^{cm}_{d} [GeV]");
            hE_d_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hE_d_cm[l][k],hpathD[l][k]);

            hp_d_lab[l][k] = new TH1F(Form("hp_d_lab_lev%d_cut%d",l,k),"",500,0.,2.);
            hp_d_lab[l][k]->GetXaxis()->SetTitle("p^{lab}_{d} [GeV/c]");
            hp_d_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hp_d_lab[l][k],hpathD[l][k]);

            hp_d_cm[l][k] = new TH1F(Form("hp_d_cm_lev%d_cut%d",l,k),"",500,0.,1.);
            hp_d_cm[l][k]->GetXaxis()->SetTitle("p^{cm}_{d} [GeV/c]");
            hp_d_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hp_d_cm[l][k],hpathD[l][k]);

            hTheta_d_lab[l][k] = new TH1F(Form("hTheta_d_lab_lev%d_cut%d",l,k),"",720,0.,60.);
            hTheta_d_lab[l][k]->GetXaxis()->SetTitle("#theta^{lab}_{d} [deg]");
            hTheta_d_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_d_lab[l][k],hpathD[l][k]);

            hTheta_d_cm[l][k] = new TH1F(Form("hTheta_d_cm_lev%d_cut%d",l,k),"",360,0.,180.);
            hTheta_d_cm[l][k]->GetXaxis()->SetTitle("#theta^{cm}_{d} [deg]");
            hTheta_d_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_d_cm[l][k],hpathD[l][k]);

            hPhi_d_lab[l][k] = new TH1F(Form("hPhi_d_lab_lev%d_cut%d",l,k),"",360,-180.,180.);
            hPhi_d_lab[l][k]->GetXaxis()->SetTitle("#phi^{lab}_{d} [deg]");
            hPhi_d_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hPhi_d_lab[l][k],hpathD[l][k]);

            hPhi_d_cm[l][k] = new TH1F(Form("hPhi_d_cm_lev%d_cut%d",l,k),"",360,-180.,180.);
            hPhi_d_cm[l][k]->GetXaxis()->SetTitle("#phi^{cm}_{d} [deg]");
            hPhi_d_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hPhi_d_cm[l][k],hpathD[l][k]);

            hEkin_vs_Theta_d_lab[l][k] = new TH2F(Form("hEkin_vs_Theta_d_lab_lev%d_cut%d",l,k),"",500,0.,1.,500,0.,60.);
            hEkin_vs_Theta_d_lab[l][k]->GetXaxis()->SetTitle("E^{kin}_{d} [GeV]");
            hEkin_vs_Theta_d_lab[l][k]->GetYaxis()->SetTitle("#theta^{lab}_{d} [deg]");
            gHistoManager->Add(hEkin_vs_Theta_d_lab[l][k],hpathD[l][k]);

            hEkin_vs_Theta_d_cm[l][k] = new TH2F(Form("hEkin_vs_Theta_d_cm_lev%d_cut%d",l,k),"",500,0.,0.1,500,0.,180.);
            hEkin_vs_Theta_d_cm[l][k]->GetXaxis()->SetTitle("E^{kin,cm}_{d} [GeV]");
            hEkin_vs_Theta_d_cm[l][k]->GetYaxis()->SetTitle("#theta^{cm}_{d} [deg]");
            gHistoManager->Add(hEkin_vs_Theta_d_cm[l][k],hpathD[l][k]);

            //gamma #1//
            hp_g1_lab[l][k] = new TH1F(Form("hp_g1_lab_lev%d_cut%d",l,k),"",500,0.,0.8);
            hp_g1_lab[l][k]->GetXaxis()->SetTitle("p^{lab}_{#gamma_{1}} [GeV/c]");
            hp_g1_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hp_g1_lab[l][k],hpathG[l][k]);

            hp_g1_cm[l][k] = new TH1F(Form("hp_g1_cm_lev%d_cut%d",l,k),"",500,0.,0.8);
            hp_g1_cm[l][k]->GetXaxis()->SetTitle("p^{cm}_{#gamma_{1}} [GeV/c]");
            hp_g1_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hp_g1_cm[l][k],hpathG[l][k]);

            hTheta_g1_lab[l][k] = new TH1F(Form("hTheta_g1_lab_lev%d_cut%d",l,k),"",360,0.,180.0);
            hTheta_g1_lab[l][k]->GetXaxis()->SetTitle("#theta^{lab}_{#gamma_{1}} [deg]");
            hTheta_g1_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_g1_lab[l][k],hpathG[l][k]);

            hTheta_g1_cm[l][k] = new TH1F(Form("hTheta_g1_cm_lev%d_cut%d",l,k),"",360,0.,180.0);
            hTheta_g1_cm[l][k]->GetXaxis()->SetTitle("#theta^{cm}_{#gamma_{1}} [deg]");
            hTheta_g1_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_g1_cm[l][k],hpathG[l][k]);

            hPhi_g1_lab[l][k] = new TH1F(Form("hPhi_g1_lab_lev%d_cut%d",l,k),"",360,-180.,180.);
            hPhi_g1_lab[l][k]->GetXaxis()->SetTitle("#phi^{lab}_{#gamma_{1}} [deg]");
            hPhi_g1_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hPhi_g1_lab[l][k],hpathG[l][k]);

            hPhi_g1_cm[l][k] = new TH1F(Form("hPhi_g1_cm_lev%d_cut%d",l,k),"",360,-180.,180.);
            hPhi_g1_cm[l][k]->GetXaxis()->SetTitle("#phi^{cm}_{#gamma_{1}} [deg]");
            hPhi_g1_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hPhi_g1_cm[l][k],hpathG[l][k]);

            hEkin_vs_Theta_g1_lab[l][k] = new TH2F(Form("hEkin_vs_Theta_g1_lab_lev%d_cut%d",l,k),"",500,0.,0.8,500,0.,180.);
            hEkin_vs_Theta_g1_lab[l][k]->GetXaxis()->SetTitle("E^{kin}_{#gamma_{1}} [GeV]");
            hEkin_vs_Theta_g1_lab[l][k]->GetYaxis()->SetTitle("#theta^{lab}_{#gamma_{1}} [deg]");
            gHistoManager->Add(hEkin_vs_Theta_g1_lab[l][k],hpathG[l][k]);

            hEkin_vs_Theta_g1_cm[l][k] = new TH2F(Form("hEkin_vs_Theta_g1_cm_lev%d_cut%d",l,k),"",500,0.,0.8,500,0.,180.);
            hEkin_vs_Theta_g1_cm[l][k]->GetXaxis()->SetTitle("E^{kin,cm}_{#gamma_{1}} [GeV]");
            hEkin_vs_Theta_g1_cm[l][k]->GetYaxis()->SetTitle("#theta^{cm}_{#gamma_{1}} [deg]");
            gHistoManager->Add(hEkin_vs_Theta_g1_cm[l][k],hpathG[l][k]);

            //gamma #2//
            hp_g2_lab[l][k] = new TH1F(Form("hp_g2_lab_lev%d_cut%d",l,k),"",500,0.,0.8);
            hp_g2_lab[l][k]->GetXaxis()->SetTitle("p^{lab}_{#gamma_{2}} [GeV/c]");
            hp_g2_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hp_g2_lab[l][k],hpathG[l][k]);

            hp_g2_cm[l][k] = new TH1F(Form("hp_g2_cm_lev%d_cut%d",l,k),"",500,0.,0.8);
            hp_g2_cm[l][k]->GetXaxis()->SetTitle("p^{cm}_{#gamma_{2}} [GeV/c]");
            hp_g2_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hp_g2_cm[l][k],hpathG[l][k]);

            hTheta_g2_lab[l][k] = new TH1F(Form("hTheta_g2_lab_lev%d_cut%d",l,k),"",360,0.,180.0);
            hTheta_g2_lab[l][k]->GetXaxis()->SetTitle("#theta^{lab}_{#gamma_{2}} [deg]");
            hTheta_g2_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_g2_lab[l][k],hpathG[l][k]);

            hTheta_g2_cm[l][k] = new TH1F(Form("hTheta_g2_cm_lev%d_cut%d",l,k),"",360,0.,180.0);
            hTheta_g2_cm[l][k]->GetXaxis()->SetTitle("#theta^{cm}_{#gamma_{2}} [deg]");
            hTheta_g2_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_g2_cm[l][k],hpathG[l][k]);

            hPhi_g2_lab[l][k] = new TH1F(Form("hPhi_g2_lab_lev%d_cut%d",l,k),"",360,-180.,180.);
            hPhi_g2_lab[l][k]->GetXaxis()->SetTitle("#phi^{lab}_{#gamma_{2}} [deg]");
            hPhi_g2_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hPhi_g2_lab[l][k],hpathG[l][k]);

            hPhi_g2_cm[l][k] = new TH1F(Form("hPhi_g2_cm_lev%d_cut%d",l,k),"",360,-180.,180.);
            hPhi_g2_cm[l][k]->GetXaxis()->SetTitle("#phi^{cm}_{#gamma_{2}} [deg]");
            hPhi_g2_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hPhi_g2_cm[l][k],hpathG[l][k]);

            hEkin_vs_Theta_g2_lab[l][k] = new TH2F(Form("hEkin_vs_Theta_g2_lab_lev%d_cut%d",l,k),"",500,0.,0.8,500,0.,180.);
            hEkin_vs_Theta_g2_lab[l][k]->GetXaxis()->SetTitle("E^{kin}_{#gamma_{2}} [GeV]");
            hEkin_vs_Theta_g2_lab[l][k]->GetYaxis()->SetTitle("#theta^{lab}_{#gamma_{2}} [deg]");
            gHistoManager->Add(hEkin_vs_Theta_g2_lab[l][k],hpathG[l][k]);

            hEkin_vs_Theta_g2_cm[l][k] = new TH2F(Form("hEkin_vs_Theta_g2_cm_lev%d_cut%d",l,k),"",500,0.,0.8,500,0.,180.);
            hEkin_vs_Theta_g2_cm[l][k]->GetXaxis()->SetTitle("E^{kin,cm}_{#gamma_{2}} [GeV]");
            hEkin_vs_Theta_g2_cm[l][k]->GetYaxis()->SetTitle("#theta^{cm}_{#gamma_{2}} [deg]");
            gHistoManager->Add(hEkin_vs_Theta_g2_cm[l][k],hpathG[l][k]);

            ////
            hOpeningAngle_g1_g2_lab[l][k] =new TH1F(Form("hOpeningAngle_g1_g2_lab_lev%d_cut%d",l,k),"",360,0.,180.);
            hOpeningAngle_g1_g2_lab[l][k]->GetXaxis()->SetTitle("#theta^{lab}_{#gamma,#gamma} [deg]");
            hOpeningAngle_g1_g2_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hOpeningAngle_g1_g2_lab[l][k],hpathG[l][k]);

            hOpeningAngle_g1_g2_cm[l][k] =new TH1F(Form("hOpeningAngle_g1_g2_cm_lev%d_cut%d",l,k),"",360,0.,180.);
            hOpeningAngle_g1_g2_cm[l][k]->GetXaxis()->SetTitle("#theta^{cm}_{#gamma,#gamma} [deg]");
            hOpeningAngle_g1_g2_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hOpeningAngle_g1_g2_cm[l][k],hpathG[l][k]);

            hTheta_g1_vs_Theta_g2_lab[l][k] = new TH2F(Form("hTheta_g1_vs_Theta_g2_lab_lev%d_cut%d",l,k),"",360,0.,180.,360,0.,180.);
            hTheta_g1_vs_Theta_g2_lab[l][k]->GetXaxis()->SetTitle("#theta^{lab}_{#gamma_{1}} [deg]");
            hTheta_g1_vs_Theta_g2_lab[l][k]->GetYaxis()->SetTitle("#theta^{lab}_{#gamma_{2}} [deg]");
            gHistoManager->Add(hTheta_g1_vs_Theta_g2_lab[l][k],hpathG[l][k]);

            hTheta_g1_vs_Theta_g2_cm[l][k] = new TH2F(Form("hTheta_g1_vs_Theta_g2_cm_lev%d_cut%d",l,k),"",360,0.,180.,360,0.,180.);
            hTheta_g1_vs_Theta_g2_cm[l][k]->GetXaxis()->SetTitle("#theta^{cm}_{#gamma_{1}} [deg]");
            hTheta_g1_vs_Theta_g2_cm[l][k]->GetYaxis()->SetTitle("#theta^{cm}_{#gamma_{2}} [deg]");
            gHistoManager->Add(hTheta_g1_vs_Theta_g2_cm[l][k],hpathG[l][k]);

            //pion//
            hEkin_pi0_lab[l][k] = new TH1F(Form("hEkin_pi0_lab_lev%d_cut%d",l,k),"",500,0.,1.);
            hEkin_pi0_lab[l][k]->GetXaxis()->SetTitle("E^{kin}_{#pi^{0}} [GeV]");
            hEkin_pi0_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hEkin_pi0_lab[l][k],hpathPi[l][k]);

            hEkin_pi0_cm[l][k] = new TH1F(Form("hEkin_pi0_cm_lev%d_cut%d",l,k),"",500,0.,1.);
            hEkin_pi0_cm[l][k]->GetXaxis()->SetTitle("E^{kin}_{#pi^{0}} [GeV]");
            hEkin_pi0_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hEkin_pi0_cm[l][k],hpathPi[l][k]);

            hE_pi0_lab[l][k] = new TH1F(Form("hE_pi0_lab_lev%d_cut%d",l,k),"",500,0.,1.);
            hE_pi0_lab[l][k]->GetXaxis()->SetTitle("E^{lab}_{#pi^{0}} [GeV]");
            hE_pi0_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hE_pi0_lab[l][k],hpathPi[l][k]);

            hE_pi0_cm[l][k] = new TH1F(Form("hE_pi0_cm_lev%d_cut%d",l,k),"",500,0.,1.);
            hE_pi0_cm[l][k]->GetXaxis()->SetTitle("E^{cm}_{#pi^{0}} [GeV]");
            hE_pi0_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hE_pi0_cm[l][k],hpathPi[l][k]);

            hp_pi0_lab[l][k] = new TH1F(Form("hp_pi0_lab_lev%d_cut%d",l,k),"",500,0.,1.);
            hp_pi0_lab[l][k]->GetXaxis()->SetTitle("p^{lab}_{#pi^{0}} [GeV/c]");
            hp_pi0_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hp_pi0_lab[l][k],hpathPi[l][k]);

            hp_pi0_cm[l][k] = new TH1F(Form("hp_pi0_cm_lev%d_cut%d",l,k),"",500,0.,1.);
            hp_pi0_cm[l][k]->GetXaxis()->SetTitle("p^{cm}_{#pi^{0}} [GeV/c]");
            hp_pi0_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hp_pi0_cm[l][k],hpathPi[l][k]);

            hTheta_pi0_lab[l][k] = new TH1F(Form("hTheta_pi0_lab_lev%d_cut%d",l,k),"",360,0.,180.0);
            hTheta_pi0_lab[l][k]->GetXaxis()->SetTitle("#theta^{lab}_{#pi^{0}} [deg]");
            hTheta_pi0_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_pi0_lab[l][k],hpathPi[l][k]);

            hTheta_pi0_cm[l][k] = new TH1F(Form("hTheta_pi0_cm_lev%d_cut%d",l,k),"",360,0.,180.0);
            hTheta_pi0_cm[l][k]->GetXaxis()->SetTitle("#theta^{cm}_{#pi^{0}} [deg]");
            hTheta_pi0_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_pi0_cm[l][k],hpathPi[l][k]);

            hPhi_pi0_lab[l][k] = new TH1F(Form("hPhi_pi0_lab_lev%d_cut%d",l,k),"",360,-180.,180.);
            hPhi_pi0_lab[l][k]->GetXaxis()->SetTitle("#phi^{lab}_{#pi^{0}} [deg]");
            hPhi_pi0_lab[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hPhi_pi0_lab[l][k],hpathPi[l][k]);

            hPhi_pi0_cm[l][k] = new TH1F(Form("hPhi_pi0_cm_lev%d_cut%d",l,k),"",360,-180.,180.);
            hPhi_pi0_cm[l][k]->GetXaxis()->SetTitle("#phi^{cm}_{#pi^{0}} [deg]");
            hPhi_pi0_cm[l][k]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hPhi_pi0_cm[l][k],hpathPi[l][k]);

            hEkin_vs_Theta_pi0_lab[l][k] = new TH2F(Form("hEkin_vs_Theta_pi0_lab_lev%d_cut%d",l,k),"",500,0.,0.8,500,0.,180.);
            hEkin_vs_Theta_pi0_lab[l][k]->GetXaxis()->SetTitle("E^{kin}_{#pi^{0}} [GeV]");
            hEkin_vs_Theta_pi0_lab[l][k]->GetYaxis()->SetTitle("#theta^{lab}_{#pi^{0}} [deg]");
            gHistoManager->Add(hEkin_vs_Theta_pi0_lab[l][k],hpathPi[l][k]);

            hEkin_vs_Theta_pi0_cm[l][k] = new TH2F(Form("hEkin_vs_Theta_pi0_cm_lev%d_cut%d",l,k),"",500,0.,0.8,500,0.,180.);
            hEkin_vs_Theta_pi0_cm[l][k]->GetXaxis()->SetTitle("E^{kin,cm}_{#pi^{0}} [GeV]");
            hEkin_vs_Theta_pi0_cm[l][k]->GetYaxis()->SetTitle("#theta^{cm}_{#pi^{0}} [deg]");
            gHistoManager->Add(hEkin_vs_Theta_pi0_cm[l][k],hpathPi[l][k]);

        }

    }

    ////////////////////////
    for (Int_t q = 3; q < 5; q++) {
        for (Int_t r = 1; r < 3; r++) {
            hQ[q][r][4][0] = new TH1F(Form("hQ_lev%d_%d_cut4",q-1,r),"",40,-70.,30.);
            hQ[q][r][4][0]->GetXaxis()->SetTitle("Q [MeV]");
            hQ[q][r][4][0]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hQ[q][r][4][0],"Systematcs");

            hQ[q][r][5][0] = new TH1F(Form("hQ_lev%d_%d_cut5",q-1,r),"",40,-70.,30.);
            hQ[q][r][5][0]->GetXaxis()->SetTitle("Q [MeV]");
            hQ[q][r][5][0]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hQ[q][r][5][0],"Systematcs");

        }
    }

    for (Int_t l = 2; l < 4; l++) {
        for (Int_t m = 4; m < 6; m++) {
            for (Int_t n = 1; n < 5; n++) {
                for (Int_t p = 1; p < 3; p++) {
                    hQ[l][m][n][p] = new TH1F(Form("hQ_lev%d_cut%d_%d_%d",l,m,n,p),"",40,-70.,30.);
                    hQ[l][m][n][p]->GetXaxis()->SetTitle("Q [MeV]");
                    hQ[l][m][n][p]->GetYaxis()->SetTitle("counts");
                    gHistoManager->Add(hQ[l][m][n][p],"Systematcs");
                }
            }
        }
    }

    ////statistics [h_st]/////

    hStatistics[0] = new TH1F("Statistics_trigger","",11,-0.5,10.5);
    hStatistics[0]->GetXaxis()->SetTitle("cut");
    hStatistics[0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hStatistics[0],h_st);

    hStatistics[1] = new TH1F("Statistics_level_2","",11,-0.5,10.5);
    hStatistics[1]->GetXaxis()->SetTitle("cut");
    hStatistics[1]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hStatistics[1],h_st);

    hStatistics[2] = new TH1F("Statistics_level_3","",11,-0.5,10.5);
    hStatistics[2]->GetXaxis()->SetTitle("cut");
    hStatistics[2]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hStatistics[2],h_st);

    return;

}   //02//

void eventselection::Clear(Option_t *option){
    fProcessed = kFALSE;
    return;
}

void eventselection::Print(Option_t *option){
    return;
}

void eventselection::UserCommand(CCommand * command){
    return;
}
