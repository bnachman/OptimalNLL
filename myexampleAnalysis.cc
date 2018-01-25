#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>


#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "myexampleAnalysis.h"
#include "myTools.h"

#include "myFastJetBase.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "IterativeSoftDrop.h"
//#include "fastjet/contrib/IteratedSoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"

#include "Pythia8/Pythia.h"

using namespace std;

double Rval = 0.8;

// Constructor 
myexampleAnalysis::myexampleAnalysis(){
    if(fDebug) cout << "myexampleAnalysis::myexampleAnalysis Start " << endl;
    ftest = 0;
    fDebug = false;
    fOutName = "test.root";
    tool = new myTools();

    // jet def 
    m_jet_def               = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rval);
    m_jet_def_ca  = new fastjet::JetDefinition(fastjet::cambridge_algorithm, Rval);

    if(fDebug) cout << "myexampleAnalysis::myexampleAnalysis End " << endl;
}

// Destructor 
myexampleAnalysis::~myexampleAnalysis(){
    delete tool;
    delete m_jet_def;
    delete m_jet_def_ca;
}

// Begin method
void myexampleAnalysis::Begin(){
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("EventTree", "Event Tree for myexample");
   int nbins = 80;
   lund = new TH2F("lund","lund",nbins,0,10,nbins,0,10);
   lundg = new TH2F("lundg","lundg",nbins,0,10,nbins,0,10);
   lundq = new TH2F("lundq","lundq",nbins,0,10,nbins,0,10);
   lund_p = new TH2F("lund_p","lund_p",nbins,0,10,nbins,0,10);
   lundg_p = new TH2F("lundg_p","lundg_p",nbins,0,10,nbins,0,10);
   lundq_p = new TH2F("lundq_p","lundq_p",nbins,0,10,nbins,0,10);
   integrals = new TH1F("integrals","integrals",6,0.5,6.5);
   DeclareBranches();
   ResetBranches();   

   return;
}

// End
void myexampleAnalysis::End(){
    
    tT->Write();
    lund->Write();
    lundg->Write();
    lundq->Write();
    lund_p->Write();
    lundg_p->Write();
    lundq_p->Write();
    integrals->Write();
    tF->Close();
    return;
}

// Analyze
void myexampleAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8){
    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent Begin " << endl;

    // -------------------------
    if (!pythia8->next()) return;
    //pythia8->event.list();

    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent Event Number " << ievt << endl;

    // reset branches 
    ResetBranches();
    
    // new event-----------------------
    fTEventNumber = ievt;
    std::vector <fastjet::PseudoJet>           particlesForJets;
    std::vector <fastjet::PseudoJet>           HS;
    
    // Particle loop -----------------------------------------------------------
    for (unsigned int ip=0; ip<pythia8->event.size(); ++ip){

      fastjet::PseudoJet p(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e());
      p.set_user_info(new MyUserInfo(pythia8->event[ip].id(),ip,pythia8->event[ip].charge()));

      if (fabs(pythia8->event[ip].id()) < 40){
	HS.push_back(p);
      }
      
      // particles for jets --------------
      if (!pythia8->event[ip].isFinal() )      continue;
      if (fabs(pythia8->event[ip].id())  ==12) continue;
      if (fabs(pythia8->event[ip].id())  ==13) continue;
      if (fabs(pythia8->event[ip].id())  ==14) continue;
      if (fabs(pythia8->event[ip].id())  ==16) continue;

      particlesForJets.push_back(p);
      
    } // end particle loop -----------------------------------------------

    fastjet::ClusterSequence cs(particlesForJets, *m_jet_def);
    vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(500.0));
    for (unsigned int ij = 0; ij < jets.size(); ij++){ 

      fTJsmallPt[fTNJetsSmallRFilled] = jets[ij].perp();
      fTJsmallN[fTNJetsSmallRFilled] = 0;
      fTJsmallSquirrel[fTNJetsSmallRFilled] = 0.;
      fTNJetsSmallRFilled+=1;

      double Emax = -1;
      int flav = 0;
      for (unsigned int i=0; i < HS.size(); i++){
	if (HS[i].delta_R(jets[ij]) > m_jet_def->R()) continue;
	if (HS[i].e() > Emax){
	  flav = HS[i].user_info<MyUserInfo>().pdg_id();
	  Emax = HS[i].e();
	}
      }

      fTJsmalltype[fTNJetsSmallRFilled-1] = flav;
      fTJsmallTrackMult[fTNJetsSmallRFilled-1]=0;
      for (int k=0; k<jets[ij].constituents().size(); k++){
	//std::cout << k << " " << jets[ij].constituents()[k].user_info<MyUserInfo>().charge() << std::endl;
	if (jets[ij].constituents()[k].perp() > 0.5 && fabs(jets[ij].constituents()[k].user_info<MyUserInfo>().charge()) > 0) fTJsmallTrackMult[fTNJetsSmallRFilled-1]+=1;//note this is nonsense for parton level!
      }

      for (int k=0; k<15; k++){
	fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][k]=0.;
      }
      
      fastjet::ClusterSequence cs2(jets[ij].constituents(), *m_jet_def_ca);
      vector<fastjet::PseudoJet> jets2 = fastjet::sorted_by_pt(cs2.inclusive_jets(0.0));
      fastjet::PseudoJet j = jets2[0];
      fastjet::PseudoJet jj, j1, j2;
      jj = j;

      fastjet::contrib::IterativeSoftDrop isd( 0.007, -1, 0, 0.8 );
      int sdemissions = isd.result(j);
      fTJsmallISD[fTNJetsSmallRFilled-1] = sdemissions;

      integrals->Fill(1);
      if (flav==21) integrals->Fill(2);
      else integrals->Fill(3);

      int mycounter = 0;
      while (jj.has_parents(j1,j2)) {
        // make sure j1 is always harder branch
        if (j1.pt2() < j2.pt2()) swap(j1,j2);
	
	mycounter+=1;

        // collect info and fill in the histogram
        double delta_R = j1.delta_R(j2);
        double delta_R_norm = delta_R / m_jet_def_ca->R();
        double z = j2.pt()/(j1.pt() + j2.pt());
        double y = log(1.0 / delta_R_norm);

        // there is an ambiguity here: can use z or j2.pt() / j.pt()
        double lnpt_rel = log(z * delta_R_norm);
        double lnpt_abs = log(j2.pt()/j.pt() * delta_R_norm);
        
        //hists_2d["lund-zrel"].add_entry(y, lnpt_rel, evwgt);
        //hists_2d["lund-zabs"].add_entry(y, lnpt_abs, evwgt);

        double lnpt = log(j2.pt() * delta_R);
        //hists_2d["lund-lnpt"].add_entry(y, lnpt, evwgt);
        // follow harder branch
	lund->Fill(y,log(1./z));
	if (flav==21) lundg->Fill(y,log(1./z));
	else lundq->Fill(y,log(1./z));
	//std::cout << log(1./z) << " " << y << std::endl;

	//std::cout << "    " << mycounter << " " << jets[ij].constituents().size() << " " << z << " " << delta_R_norm << " " << z*delta_R_norm << " " << 2./(j1.pt() + j2.pt()) << std::endl;
	if (z*delta_R_norm > 2./(j1.pt() + j2.pt())){
	  fTJsmallN[fTNJetsSmallRFilled-1]+=1.;
	  fTJsmallSquirrel[fTNJetsSmallRFilled-1]+=(1.-1.5*z*z/log(4./9.));

	  /*
	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][0] += (1.-4.*1.5*z*z/log(4./9.));
	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][1] += (1.-2.*1.5*z*z/log(4./9.));
	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][2] += (1.-0.5*1.5*z*z/log(4./9.));
	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][3] += (1.-0.25*1.5*z*z/log(4./9.));
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][4] += (1.+0.25*1.5*z*z/log(4./9.));
	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][5] += (1.+0.5*1.5*z*z/log(4./9.));
	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][6] += (1.+1.*1.5*z*z/log(4./9.));
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][7] += (1.+2.*1.5*z*z/log(4./9.));
	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][8] += (1.+4.*1.5*z*z/log(4./9.));
	  
	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][9] += (1. - 0.5*z);
	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][10] += (1. - 0.25*z);
	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][11] += (1. - 0.125*z);
	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][12] += (1. + 0.5*z);
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][13] += (1. + 0.25*z);
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][14] += (1. + 0.125*z);
	  */

	  fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][0] += (1.-50.*1.5*z*z/log(4./9.));
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][1] += (1.-20.*1.5*z*z/log(4./9.));
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][2] += (1.-10.*1.5*z*z/log(4./9.));
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][3] += (1.-5.*1.5*z*z/log(4./9.));
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][4] += (1.+1.*1.5*z*z/log(4./9.));
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][5] += (1.+5.*1.5*z*z/log(4./9.));
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][6] += (1.+10.*1.5*z*z/log(4./9.));
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][7] += (1.+20.*1.5*z*z/log(4./9.));
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][8] += (0.+50.*1.5*z*z/log(4./9.));

          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][9] += (1. - 100.*z);
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][10] += (1. - 10.*z);
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][11] += (1. - 1.*z);
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][12] += (1. + 1.*z);
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][13] += (1. + 10.*z);
          fTJsmallSquirrel_many[fTNJetsSmallRFilled-1][14] += (0. + 100.*z);

	  lund_p->Fill(y,log(1./z));
	  if (flav==21) lundg_p->Fill(y,log(1./z));
	  else lundq_p->Fill(y,log(1./z));
	}

        jj = j1;
      }
      //std::cout << "blah blah " << ij << " " << sdemissions << " " << mycounter << std::endl;
    }

    tT->Fill();

    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent End " << endl;
    return;
}



// declate branches
void myexampleAnalysis::DeclareBranches(){
   
   // Event Properties 
   tT->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
   
   // smallR jets
   tT->Branch("NJetsFilledSmallR",         &fTNJetsSmallRFilled,       "NJetsFilledSmallR/I");
   tT->Branch("JsmallPt",                  &fTJsmallPt,                "JsmallPt[NJetsFilledSmallR]/F");
   tT->Branch("JsmallN",                  &fTJsmallN,                "JsmallN[NJetsFilledSmallR]/F");
   tT->Branch("JsmallSquirrel",                  &fTJsmallSquirrel,                "JsmallSquirrel[NJetsFilledSmallR]/F");
   tT->Branch("Jsmalltype",                  &fTJsmalltype,                "Jsmalltype[NJetsFilledSmallR]/I");
   tT->Branch("JsmallTrackMult",                  &fTJsmallTrackMult,                "JsmallTrackMult[NJetsFilledSmallR]/I");
   tT->Branch("JsmallISD",                  &fTJsmallISD,                "JsmallISD[NJetsFilledSmallR]/I");
   tT->Branch("JsmallSquirrel_many",                  &fTJsmallSquirrel_many,                "JsmallSquirrel_many[NJetsFilledSmallR][15]/F");

   tT->GetListOfBranches()->ls();
    
   return;
}


// resets vars
void myexampleAnalysis::ResetBranches(){
      // reset branches 
      fTEventNumber                 = -999;

      fTNJetsSmallRFilled=0;
      for (int iP=0; iP < MaxNJetSmallR; ++iP){
          fTJsmallPt      [iP]= -999;
	  fTJsmallN[iP] = -999;
	  fTJsmallSquirrel[iP] = -999;
	  fTJsmalltype[iP] = -999;
	  fTJsmallTrackMult[iP] = -999;
	  fTJsmallISD[iP] = -999;
	  for (int j=0; j<15; j++){
	    fTJsmallSquirrel_many[iP][j] = -999;
	  }
      }
}
