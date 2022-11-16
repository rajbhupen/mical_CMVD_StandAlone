#ifndef MULTISIMDIGI_H
#define MULTISIMDIGI_H 1
#include <vector>
using std::vector;

#include "ParameterMessenger.hh"
#include "micalDetectorParameterDef.hh"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TTree.h"
#include "TFile.h"
// #include "evetree.h"
#include "RPCEve.h"
#include "globals.hh"
#include "Hits.h"
#include "HitPos.h"
#include "TProfile.h"
#include <iostream>
#include <fstream>
using namespace std;

struct vectGr{
  float x;
  float y;
  float z;
  float dx;
  float dy;
  float dz;
};

class MultiSimAnalysisDigi {
public:
  MultiSimAnalysisDigi();
  ~MultiSimAnalysisDigi();
  static MultiSimAnalysisDigi* AnPointer;

  void OpenInputRootFiles(char* inf);
  void OpenOutputRootFiles(char* outf);
  void OpenCollatedRootFile();
  void CloseInputRootFiles();
  void CloseOutputRootFiles();
  void SaveGenVisFile();

  void SetCorrTimeError(G4double val);
  void SetUnCorrTimeError(G4double val);
  void SetTimeToDigiConvVal(G4double val);
  void SetSignalSpeedVal(G4double val);
  void SetPhotonSpeedVal(G4double val);
  void SetCMVadctons(G4double val);
  double GetCorrTimeError() {return CorrTimeError;}
  double GetUnCorrTimeError() {return UnCorrTimeError;}
  double GetTimeToDigiConvVal() {return TimeToDigiConv;}
  double GetSignalSpeedVal() {return SignalSpeed;}
  double GetPhotonSpeedVal() {return PhotonSpeed;}
  double GetCMVadctons() {return CMVadctons;}
  int isInOut;
  int isVisOut;
  int isXtermOut;
  int collatedIn;

  TH1F *strpXtime;
  TH1F *strpYtime;
  TH1F *strpXtimeCorr;
  TH1F *strpYtimeCorr;
  TH1F *hitXtime;
  TH1F *hitYtime;
  TH1D* smagFieldX;
  TH1D* smagFieldY;
  TH2D* smag2dX;
  TH2D* smag2dY;
  TH1D* rmagFieldX;
  TH1D* rmagFieldY;
  TH2D* rmag2dX;
  TH2D* rmag2dY;
  TH2D* smag2dXYpixel_iron;
  TH2D* smag2dXYpixel_air;
  TH2D* rmag2dXYpixel_iron;
  TH2D* rmag2dXYpixel_air;
  TH2D* xyvsbxin;
  TH2D* xyvsbyin;
  TH2D* xyvsbxdiff;
  TH2D* xyvsbydiff;
  TH2D* xyvsbxindiff;
  TH2D* xyvsbyindiff;
  TH2D* xyvsbxout;
  TH2D* xyvsbyout;
  // Collated Histograms
  TH2D* inefficiency_corx[20];
  TH2D* inefficiency_uncx[20];
  TH2D* inefficiency_uncy[20];
  TH2D* triggereffi_xevt[20];
  TH2D* triggereffi_yevt[20];
  TH2D* strp_xmulsim_cor[20];
  TH2D* strp_ymulsim_cor[20];

  TH2D* block_xmulsim[20][16][16];
  TH2D* block_ymulsim[20][16][16];







  TH1F *DeadStripX;
  TH1F *NoisyStripX;
  TH1F *DeadStripY;
  TH1F *NoisyStripY;
  TH1F *DiffTime;


  TH2F   *RC;
  TH1F   *DGap;
  TH1D* hdifftime1[20];
  TH1D* hdifftime2[20];
  TH1D* hxtime_ext[20];
  TH1D* hytime_ext[20];

  TH1D* hxpos_ext[20];
  TH1D* hypos_ext[20];
  TH1D* hxpos_ext_kalman[20];
  TH1D* hypos_ext_kalman[20];
  TH1D* h_hit_time_ext[20];

  TH1D* xtdc_minus_ref[20][8];
  TH1D* ytdc_minus_ref[20][8];

  TH1D* tshift_xtdc_minus_ref[20][8];
  TH1D* tshift_ytdc_minus_ref[20][8];

  Hits *H;
  HitPos *Hp;
  int EveCnt;
  int nloops;
  TH1F  *pdedz[20];
  TH1F  *hitDist;    //asm
  TH1F  *TrkDist;    //asm
  TH1F  *EffDist;    //asm
  TH1F  *InoTrack_listsize;   //asm
  TFile *pRootFile;
  TFile *inputRootFile;
  TFile *pVisFile;
  TFile* collatedRootFile;

  TTree *pEventTree;
  TTree *inputEventTree;
  TTree *visTree;
  RPCEve *data_event;
  // evetree *data_event;

  TH1F* ShwXw;
  TH1F* ShwYw;
  TH1D* trk_gap;
  TH2D* trk_edge;
  TH1F* pPosX;
  TH1F* pPosY;
  TH1F* pPosZ;
  TH2F* pPosXX;
  TH2F* pPosYY;
  TH2F* pPosZZ;


  static const int nhistmx=1000;
  int   ihist;
  TH3F* gens_list[6][nhistmx]; // Not used
  vector<vectGr> gens_vect[6];

  // Common input read and output store for ICALsim and miniICALsim
  static const unsigned int ngenmx=50;
  UInt_t          irun;
  UInt_t          ievt;
  UInt_t          ievt2;
  UInt_t          ievt3;
  UInt_t          ngent;
  Int_t           pidin[ngenmx];   //[ngent]
  Float_t         ievt_wt;
  Int_t           intxn_id;
  Float_t         momin[ngenmx];   //[ngent]
  Float_t         thein[ngenmx];   //[ngent]
  Float_t         phiin[ngenmx];   //[ngent]
  Float_t         posxin[ngenmx];   //[ngent]
  Float_t         posyin[ngenmx];   //[ngent]
  Float_t         poszin[ngenmx];   //[ngent]

  // ICALsim Root Files Input Data Read (SIM)
  static const unsigned int nsimhtmx=4000;
  UInt_t          nsimht;
  UInt_t          detid[nsimhtmx];   //[nsimht]
  Int_t           simpdgid[nsimhtmx];   //[nsimht]
  Float_t         simtime[nsimhtmx];   //[nsimht]
  Float_t         simenr[nsimhtmx];   //[nsimht]
  Float_t         simvx[nsimhtmx];   //[nsimht]
  Float_t         simvy[nsimhtmx];   //[nsimht]
  Float_t         simvz[nsimhtmx];   //[nsimht]
  Float_t         simpx[nsimhtmx];   //[nsimht]
  Float_t         simpy[nsimhtmx];   //[nsimht]
  Float_t         simpz[nsimhtmx];   //[nsimht]
  Float_t         simlocvx[nsimhtmx];   //[nsimht]
  Float_t         simlocvy[nsimhtmx];   //[nsimht]
  Float_t         simlocvz[nsimhtmx];   //[nsimht]
  // ICALsim Root Files Input Data Read
  static const unsigned int ndigihtmx=5000;
  UInt_t          ndigiht;
  Int_t           trigx;
  Int_t           trigy;
  UInt_t          stripid[ndigihtmx];   //[ndigiht]
  Int_t           digipdgid[ndigihtmx];   //[ndigiht]
  Int_t           digitime[ndigihtmx];   //[ndigiht]
  Int_t           digitruetime[ndigihtmx];   //[ndigiht]
  Float_t         digienr[ndigihtmx];   //[ndigiht]
  Float_t         digivx[ndigihtmx];   //[ndigiht]
  Float_t         digivy[ndigihtmx];   //[ndigiht]
  Float_t         digivz[ndigihtmx];   //[ndigiht]
  Float_t         digipx[ndigihtmx];   //[ndigiht]
  Float_t         digipy[ndigihtmx];   //[ndigiht]
  Float_t         digipz[ndigihtmx];   //[ndigiht]
  int diginoise[ndigihtmx];

  // Reco Output store for both sim.
  static const unsigned int nvishtmx=5000;
  Float_t fitposxx[nvishtmx];
  Float_t fitposyy[nvishtmx];
  Float_t fitposzz[nvishtmx];
  Float_t fitlayzz[nvishtmx];
  Float_t fitlayx2[nvishtmx];
  Float_t fitlayx3[nvishtmx];
  Float_t fitlayx4[nvishtmx];
  Float_t fitlaymom[nvishtmx];
  Float_t fitlaythe[nvishtmx];
  Float_t fitlayphi[nvishtmx];
  Float_t extrapolxx[nvishtmx];
  Float_t extrapolyy[nvishtmx];
  Float_t extrapolmom[nvishtmx];

  Float_t momdiff1;
  Float_t radialdiff1;

  unsigned int nvisht;
  Float_t clstposxx[nvishtmx];
  Float_t clstposyy[nvishtmx];
  Float_t clstposzz[nvishtmx];
  Int_t clstposzpln[nvishtmx];
  unsigned int nvisclst;

  static  const unsigned int nthtmx=100;

  Int_t ntrecord1x;
  Int_t ntrecord1y;
  Int_t ntrecord2x;
  Int_t ntrecord2y;
  Int_t striprec1x[nthtmx];
  Int_t striprec1y[nthtmx];
  Int_t striprec2x[nthtmx];
  Int_t striprec2y[nthtmx];
  Float_t tdcrec1x[nthtmx];
  Float_t tdcrec1y[nthtmx];
  Float_t tdcrec2x[nthtmx];
  Float_t tdcrec2y[nthtmx];

  Int_t nhits_last;
  Int_t nhits_last_m1;

  Int_t strtnhitsx;
  Int_t strtnhitsy;
  Float_t strtchisqx;
  Float_t strtchisqy;
  Float_t strtintercptx;
  Float_t strtintercpty;
  Float_t strtslopex;
  Float_t strtslopey;

  Float_t simpleradii;
  Float_t simplecurv;
  Float_t simplex0;
  Float_t simplez0;
  Float_t simplechisqpos;
  Float_t simplechisqneg;
  Float_t simplechisqcndn;
  Float_t simpleavgxpos;
  Float_t simpleavgxneg;
  Float_t simpleavgxcndn;
  Float_t simpleavgxmeas;
  Float_t simplenhits;

  Int_t ntdc1x;
  Int_t ntstrp1x;
  Int_t tdcID1x[nthtmx];
  Int_t StrpID1x[nthtmx];
  Float_t TDCval1x[nthtmx];

  Int_t ntdc2x;
  Int_t ntstrp2x;
  Int_t tdcID2x[nthtmx];
  Int_t StrpID2x[nthtmx];
  Float_t TDCval2x[nthtmx];

  Int_t ntdc1y;
  Int_t ntstrp1y;
  Int_t tdcID1y[nthtmx];
  Int_t StrpID1y[nthtmx];
  Float_t TDCval1y[nthtmx];

  Int_t ntdc2y;
  Int_t ntstrp2y;
  Int_t tdcID2y[nthtmx];
  Int_t StrpID2y[nthtmx];
  Float_t TDCval2y[nthtmx];

  Int_t nhits_below;
  Float_t ftime_last;

  static  const unsigned int ntrkmx=20;
  // Int_t           hw_trigx;
  UInt_t          ngenerated;
  UInt_t          naperture;
  UInt_t          triggeracceptance;
  Int_t           hw_trig;
  Int_t           sw_trigx;
  Int_t           sw_trigy;
  UInt_t          ntrkt;
  Int_t           itype[ntrkmx];   //[ntrkt]
  Int_t           nLayer;
  Int_t           nhits[ntrkmx];   //[ntrkt]
  Int_t           nhits_finder[ntrkmx];   //[ntrkt]
  Float_t         chisq[ntrkmx];   //[ntrkt]
  Float_t chisq2[ntrkmx];
  Int_t ndof[ntrkmx];
  Int_t ndof2[ntrkmx];
  
  Float_t         cvalue[ntrkmx];   //[ntrkt]
  Int_t           fc_or_pc[ntrkmx];   //[ntrkt]
  Float_t         trkmm[ntrkmx];   //[ntrkt]
  Float_t         trkth[ntrkmx];   //[ntrkt]
  Float_t         trkph[ntrkmx];   //[ntrkt]
  Float_t         momvx[ntrkmx];   //[ntrkt]
  Float_t         thevx[ntrkmx];   //[ntrkt]

  Float_t         phivx[ntrkmx];   //[ntrkt]
  Float_t         posxvx[ntrkmx];   //[ntrkt]
  Float_t         posyvx[ntrkmx];   //[ntrkt]
  Float_t         poszvx[ntrkmx];   //[ntrkt]

  //ROOF RSA

  Float_t momrf[ntrkmx]; 	  //Measured momentum of reconstrued track
  Float_t therf[ntrkmx];	  //Measured polar angle of track
  Float_t phirf[ntrkmx];        //Measured azimuthal angle of track
  Float_t posxrf[ntrkmx];	  //Starting X-position
  Float_t posyrf[ntrkmx];       //Starting Y-position
  Float_t poszrf[ntrkmx];       //Starting Z-position

  static const unsigned int nlayermx=10;
  //for posresolution SL
  Float_t XPosdev_exclu[nlayermx];
  Float_t YPosdev_exclu[nlayermx];
  

  Float_t XPosdev[nlayermx];
  Float_t YPosdev[nlayermx];


  
  //Position Resolutions JIM
  Float_t inPosX[10];
  Float_t extPosX[10];
  Float_t inPosY[10];
  Float_t extPosY[10];
  Int_t nXStrips[10];
  Int_t nYStrips[10];

  Float_t inPosfX[10];
  Float_t extPosfX[10];
  Float_t inPosfY[10];
  Float_t extPosfY[10];
  Int_t nfXStrips[10];
  Int_t nfYStrips[10];


  
  Float_t XchisqOccu[10];
  Float_t YchisqOccu[10];
  Int_t XndfOccu[10];
  Int_t YndfOccu[10];
  Int_t nxFailOccu[10];
  Int_t nyFailOccu[10];
  

  Float_t         momend[ntrkmx];   //[ntrkt]
  Float_t         theend[ntrkmx];   //[ntrkt]
  Float_t         phiend[ntrkmx];   //[ntrkt]
  Float_t         posxend[ntrkmx];   //[ntrkt]
  Float_t         posyend[ntrkmx];   //[ntrkt]
  Int_t            strpxend[ntrkmx];   //[ntrkt]
  Int_t            strpyend[ntrkmx];   //[ntrkt]
  Float_t         poszend[ntrkmx];   //[ntrkt]
  Float_t         tx_end[ntrkmx];   //[ntrkt]
  Float_t         ty_end[ntrkmx];   //[ntrkt]
  Float_t         momds[ntrkmx];   //[ntrkt]
  Float_t         momrg[ntrkmx];   //[ntrkt]
  Float_t         mcxgnvx[ntrkmx];   //[ntrkt]
  Float_t         mcygnvx[ntrkmx];   //[ntrkt]
  Float_t         momgnvx[ntrkmx];   //[ntrkt]
  Float_t         thegnvx[ntrkmx];   //[ntrkt]
  Float_t         phignvx[ntrkmx];   //[ntrkt]
  Float_t         momgnend[ntrkmx];   //[ntrkt]
  Float_t         thegnend[ntrkmx];   //[ntrkt]
  Float_t         phignend[ntrkmx];   //[ntrkt]
  Int_t           vtxzplane[ntrkmx];   //[ntrkt]
  Int_t           endzplane[ntrkmx];   //[ntrkt]
  Int_t           ntrkcl[ntrkmx];   //[ntrkt]
  Int_t           ntrkst[ntrkmx];   //[ntrkt]
  Int_t           ntotcl;
  Int_t           ntotst;
  Int_t           inohits;
  Int_t           orighits;
  Int_t           inoclust;
  Int_t           origclust;
  Float_t         hPathlength;
  Int_t           x_hits;
  Int_t           y_hits;
  Int_t           inohits_old;
  Int_t           orighits_old;
  Int_t           x_hits_old;
  Int_t           y_hits_old;
  Int_t           hit_wo_ghst;
  Float_t         e_hadron;
  Int_t           nhits_largest_cluster;
  Int_t           orighits_trape;
  Int_t           orighits_cluster;
  Int_t           hit_wogh_orighits;
  Float_t         theta_hadron_shw;
  Float_t         had_eigen_val[3];
  Float_t         phi_hadron_shw;
  Float_t         theta_hadron_in;
  Float_t         phi_hadron_in;
  Float_t         dot_angle_had_shw;
  Int_t           nhits_largest_cluster_selected;
  Float_t         range;
  Float_t         tx[ntrkmx];
  Float_t         ty[ntrkmx];
  Float_t         xxin[ntrkmx];
  Float_t         yyin[ntrkmx];
  Float_t         txin[ntrkmx];
  Float_t         tyin[ntrkmx];
  Float_t         therr[ntrkmx];
  Float_t         pherr[ntrkmx];
  Float_t         xxerr[ntrkmx];
  Float_t         yyerr[ntrkmx];
  Float_t         txerr[ntrkmx];
  Float_t         tyerr[ntrkmx];
  Float_t         qperr[ntrkmx];
  Float_t         xxenderr[ntrkmx];
  Float_t         yyenderr[ntrkmx];
  Float_t         txenderr[ntrkmx];
  Float_t         tyenderr[ntrkmx];
  Float_t         qpenderr[ntrkmx];

  Int_t nmxhit;





  Float_t atimslope[ntrkmx];
  Float_t atiminter[ntrkmx];

  Float_t xxtxerr[ntrkmx];
  Float_t xxtyerr[ntrkmx];
  Float_t yytyerr[ntrkmx];
  Float_t yytxerr[ntrkmx];
  Float_t txtyerr[ntrkmx];
 

  Float_t atim[ntrkmx];


  
  //cmvd extrapol :

  static const unsigned int cmvdnlayermx=32;

  unsigned int CMVDnLayer;
  
  
  Float_t CMVDTruePosX[32];
  Float_t CMVDTruePosY[32];
  Float_t CMVDTruePosZ[32];

  Float_t CMVDRecoPosX[32];
  Float_t CMVDRecoPosY[32];
  Float_t CMVDRecoPosZ[32];

  Float_t CMVDWRecoPosX[32];
  Float_t CMVDWRecoPosY[32];
  Float_t CMVDWRecoPosZ[32];
  
  Float_t CMVDExpPosX[32];
  Float_t CMVDExpPosY[32];
  Float_t CMVDExpPosZ[32];

  Int_t CMVDExpnHit[32];
  Int_t CMVDisExpIn[32];
  Int_t CMVDisHit[32];


  Int_t CMVDpl2[32];
  Int_t CMVDClustsInLay[32];
  Int_t CMVDHitsInLay[32];
  Int_t CMVDClusterSize[32];
  
  
  //

 




  

  static const unsigned int cmv_nhtmx=10000;

  unsigned int cmv_nhit;
  //cmvd cluster hit:
  unsigned int cmv_hitid[ cmv_nhtmx];
  Int_t cmv_hitpdgid[ cmv_nhtmx];
  Float_t cmv_hitLeTim[ cmv_nhtmx];
  Float_t cmv_hitRiTim[ cmv_nhtmx];
  Float_t cmv_hitLePul[ cmv_nhtmx];
  Float_t cmv_hitRiPul[ cmv_nhtmx];
  Float_t cmv_hitTrueposx[ cmv_nhtmx];
  Float_t cmv_hitTrueposy[ cmv_nhtmx];
  Float_t cmv_hitTrueposz[ cmv_nhtmx];
  Float_t cmv_hitRecoposx[ cmv_nhtmx];
  Float_t cmv_hitRecoposy[ cmv_nhtmx];
  Float_t cmv_hitRecoposz[ cmv_nhtmx];


  //  Int_t cmv_hitsiz[ cmv_nhtmx];

  //cmv_ cluster:


  static const unsigned int cmv_nclustmx=10000;
  UInt_t cmv_nclust;
  //cmvd cluster hit:
  UInt_t cmv_clustid[ cmv_nclustmx];
  Int_t cmv_clustpdgid[ cmv_nclustmx];
  Float_t cmv_clustLeTim[ cmv_nclustmx];
  Float_t cmv_clustRiTim[ cmv_nclustmx];
  Float_t cmv_clustLePul[ cmv_nclustmx];
  Float_t cmv_clustRiPul[ cmv_nclustmx];
  Float_t cmv_clustTrueposx[ cmv_nclustmx];
  Float_t cmv_clustTrueposy[ cmv_nclustmx];
  Float_t cmv_clustTrueposz[ cmv_nclustmx];
  Float_t cmv_clustRecoposx[ cmv_nclustmx];
  Float_t cmv_clustRecoposy[ cmv_nclustmx];
  Float_t cmv_clustRecoposz[ cmv_nclustmx];
  Int_t cmv_clustsiz[ cmv_nclustmx];










  //cmvd extrapol :


  static const unsigned int cmv_nexphtmx=1000;

  UInt_t cmv_nexphit;


  UInt_t cmv_expid[cmv_nexphtmx];

  Float_t cmv_Expposx[cmv_nexphtmx];
  Float_t cmv_Expposy[cmv_nexphtmx];
  Float_t cmv_Expposz[cmv_nexphtmx];



  Float_t cmv_DCAposx[cmv_nexphtmx];
  Float_t cmv_DCAposy[cmv_nexphtmx];
  Float_t cmv_DCAposz[cmv_nexphtmx];



  Float_t distofclosapp[ntrkmx];
  Int_t planeedge[ntrkmx];
  Int_t clustersize00[ntrkmx];

  Float_t extrapolposx00[ntrkmx];
  Float_t extrapolposy00[ntrkmx];
  Float_t extrapolposz00[ntrkmx];

  Float_t extrapolposxerr00[ntrkmx];
  Float_t extrapolposyerr00[ntrkmx];
  Float_t extrapolposzerr00[ntrkmx];

  Float_t cmvhitrecoposx00[ntrkmx];
  Float_t cmvhitrecoposy00[ntrkmx];
  Float_t cmvhitrecoposz00[ntrkmx];

  Float_t cmvhittrueposx00[ntrkmx];
  Float_t cmvhittrueposy00[ntrkmx];
  Float_t cmvhittrueposz00[ntrkmx];

  Float_t cmvhitrecoposxerr00[ntrkmx];
  Float_t cmvhitrecoposyerr00[ntrkmx];
  Float_t cmvhitrecoposzerr00[ntrkmx];

  Float_t LeTime00[ntrkmx];
  Float_t RiTime00[ntrkmx];
  Float_t LePulse00[ntrkmx];
  Float_t RiPulse00[ntrkmx];

  //   pEventTree->Branch("LeTime31", LeTime31, "LeTime31[ntrkt]/F");
  // pEventTree->Branch("RiTime31", RiTime31, "RiTime31[ntrkt]/F");

  //   pEventTree->Branch("LePulse31", LePulse31, "LePulse31[ntrkt]/F");
  // pEventTree->Branch("RiPulse31", RiPulse31, "RiPulse31[ntrkt]/F");





  Int_t clustersize01[ntrkmx];

  Float_t extrapolposx01[ntrkmx];
  Float_t extrapolposy01[ntrkmx];
  Float_t extrapolposz01[ntrkmx];

  Float_t extrapolposxerr01[ntrkmx];
  Float_t extrapolposyerr01[ntrkmx];
  Float_t extrapolposzerr01[ntrkmx];

  Float_t cmvhitrecoposx01[ntrkmx];
  Float_t cmvhitrecoposy01[ntrkmx];
  Float_t cmvhitrecoposz01[ntrkmx];

  Float_t cmvhittrueposx01[ntrkmx];
  Float_t cmvhittrueposy01[ntrkmx];
  Float_t cmvhittrueposz01[ntrkmx];

  Float_t cmvhitrecoposxerr01[ntrkmx];
  Float_t cmvhitrecoposyerr01[ntrkmx];
  Float_t cmvhitrecoposzerr01[ntrkmx];

  Float_t LeTime01[ntrkmx];
  Float_t RiTime01[ntrkmx];
  Float_t LePulse01[ntrkmx];
  Float_t RiPulse01[ntrkmx];


  Int_t clustersize02[ntrkmx];
  Float_t extrapolposx02[ntrkmx];
  Float_t extrapolposy02[ntrkmx];
  Float_t extrapolposz02[ntrkmx];

  Float_t extrapolposxerr02[ntrkmx];
  Float_t extrapolposyerr02[ntrkmx];
  Float_t extrapolposzerr02[ntrkmx];

  Float_t cmvhitrecoposx02[ntrkmx];
  Float_t cmvhitrecoposy02[ntrkmx];
  Float_t cmvhitrecoposz02[ntrkmx];

  Float_t cmvhittrueposx02[ntrkmx];
  Float_t cmvhittrueposy02[ntrkmx];
  Float_t cmvhittrueposz02[ntrkmx];

  Float_t cmvhitrecoposxerr02[ntrkmx];
  Float_t cmvhitrecoposyerr02[ntrkmx];
  Float_t cmvhitrecoposzerr02[ntrkmx];




  Float_t LeTime02[ntrkmx];
  Float_t RiTime02[ntrkmx];
  Float_t LePulse02[ntrkmx];
  Float_t RiPulse02[ntrkmx];



  Int_t clustersize03[ntrkmx];






  Float_t extrapolposx03[ntrkmx];
  Float_t extrapolposy03[ntrkmx];
  Float_t extrapolposz03[ntrkmx];

  Float_t extrapolposxerr03[ntrkmx];
  Float_t extrapolposyerr03[ntrkmx];
  Float_t extrapolposzerr03[ntrkmx];

  Float_t cmvhitrecoposx03[ntrkmx];
  Float_t cmvhitrecoposy03[ntrkmx];
  Float_t cmvhitrecoposz03[ntrkmx];

  Float_t cmvhittrueposx03[ntrkmx];
  Float_t cmvhittrueposy03[ntrkmx];
  Float_t cmvhittrueposz03[ntrkmx];

  Float_t cmvhitrecoposxerr03[ntrkmx];
  Float_t cmvhitrecoposyerr03[ntrkmx];
  Float_t cmvhitrecoposzerr03[ntrkmx];

  Float_t LeTime03[ntrkmx];
  Float_t RiTime03[ntrkmx];
  Float_t LePulse03[ntrkmx];
  Float_t RiPulse03[ntrkmx];


  Float_t debug[ntrkmx];


  Int_t clustersize10[ntrkmx];

  Float_t extrapolposx10[ntrkmx];
  Float_t extrapolposy10[ntrkmx];
  Float_t extrapolposz10[ntrkmx];

  Float_t extrapolposxerr10[ntrkmx];
  Float_t extrapolposyerr10[ntrkmx];
  Float_t extrapolposzerr10[ntrkmx];

  Float_t cmvhitrecoposx10[ntrkmx];
  Float_t cmvhitrecoposy10[ntrkmx];
  Float_t cmvhitrecoposz10[ntrkmx];

  Float_t cmvhittrueposx10[ntrkmx];
  Float_t cmvhittrueposy10[ntrkmx];
  Float_t cmvhittrueposz10[ntrkmx];

  Float_t cmvhitrecoposxerr10[ntrkmx];
  Float_t cmvhitrecoposyerr10[ntrkmx];
  Float_t cmvhitrecoposzerr10[ntrkmx];

  Float_t LeTime10[ntrkmx];
  Float_t RiTime10[ntrkmx];
  Float_t LePulse10[ntrkmx];
  Float_t RiPulse10[ntrkmx];

  Int_t clustersize11[ntrkmx];

  Float_t extrapolposx11[ntrkmx];
  Float_t extrapolposy11[ntrkmx];
  Float_t extrapolposz11[ntrkmx];

  Float_t extrapolposxerr11[ntrkmx];
  Float_t extrapolposyerr11[ntrkmx];
  Float_t extrapolposzerr11[ntrkmx];

  Float_t cmvhitrecoposx11[ntrkmx];
  Float_t cmvhitrecoposy11[ntrkmx];
  Float_t cmvhitrecoposz11[ntrkmx];

  Float_t cmvhittrueposx11[ntrkmx];
  Float_t cmvhittrueposy11[ntrkmx];
  Float_t cmvhittrueposz11[ntrkmx];

  Float_t cmvhitrecoposxerr11[ntrkmx];
  Float_t cmvhitrecoposyerr11[ntrkmx];
  Float_t cmvhitrecoposzerr11[ntrkmx];

  Float_t LeTime11[ntrkmx];
  Float_t RiTime11[ntrkmx];
  Float_t LePulse11[ntrkmx];
  Float_t RiPulse11[ntrkmx];


  Int_t clustersize12[ntrkmx];


  Float_t extrapolposx12[ntrkmx];
  Float_t extrapolposy12[ntrkmx];
  Float_t extrapolposz12[ntrkmx];

  Float_t extrapolposxerr12[ntrkmx];
  Float_t extrapolposyerr12[ntrkmx];
  Float_t extrapolposzerr12[ntrkmx];

  Float_t cmvhitrecoposx12[ntrkmx];
  Float_t cmvhitrecoposy12[ntrkmx];
  Float_t cmvhitrecoposz12[ntrkmx];

  Float_t cmvhittrueposx12[ntrkmx];
  Float_t cmvhittrueposy12[ntrkmx];
  Float_t cmvhittrueposz12[ntrkmx];

  Float_t cmvhitrecoposxerr12[ntrkmx];
  Float_t cmvhitrecoposyerr12[ntrkmx];
  Float_t cmvhitrecoposzerr12[ntrkmx];

  Float_t LeTime12[ntrkmx];
  Float_t RiTime12[ntrkmx];
  Float_t LePulse12[ntrkmx];
  Float_t RiPulse12[ntrkmx];



  Int_t clustersize21[ntrkmx];

  Float_t extrapolposx21[ntrkmx];
  Float_t extrapolposy21[ntrkmx];
  Float_t extrapolposz21[ntrkmx];

  Float_t extrapolposxerr21[ntrkmx];
  Float_t extrapolposyerr21[ntrkmx];
  Float_t extrapolposzerr21[ntrkmx];

  Float_t cmvhitrecoposx21[ntrkmx];
  Float_t cmvhitrecoposy21[ntrkmx];
  Float_t cmvhitrecoposz21[ntrkmx];

  Float_t cmvhittrueposx21[ntrkmx];
  Float_t cmvhittrueposy21[ntrkmx];
  Float_t cmvhittrueposz21[ntrkmx];

  Float_t cmvhitrecoposxerr21[ntrkmx];
  Float_t cmvhitrecoposyerr21[ntrkmx];
  Float_t cmvhitrecoposzerr21[ntrkmx];


  Float_t LeTime21[ntrkmx];
  Float_t RiTime21[ntrkmx];
  Float_t LePulse21[ntrkmx];
  Float_t RiPulse21[ntrkmx];

  Int_t clustersize20[ntrkmx];

  Float_t extrapolposx20[ntrkmx];
  Float_t extrapolposy20[ntrkmx];
  Float_t extrapolposz20[ntrkmx];

  Float_t extrapolposxerr20[ntrkmx];
  Float_t extrapolposyerr20[ntrkmx];
  Float_t extrapolposzerr20[ntrkmx];

  Float_t cmvhitrecoposx20[ntrkmx];
  Float_t cmvhitrecoposy20[ntrkmx];
  Float_t cmvhitrecoposz20[ntrkmx];

  Float_t cmvhittrueposx20[ntrkmx];
  Float_t cmvhittrueposy20[ntrkmx];
  Float_t cmvhittrueposz20[ntrkmx];

  Float_t cmvhitrecoposxerr20[ntrkmx];
  Float_t cmvhitrecoposyerr20[ntrkmx];
  Float_t cmvhitrecoposzerr20[ntrkmx];

  Float_t LeTime20[ntrkmx];
  Float_t RiTime20[ntrkmx];
  Float_t LePulse20[ntrkmx];
  Float_t RiPulse20[ntrkmx];

  Int_t clustersize22[ntrkmx];

  Float_t extrapolposx22[ntrkmx];
  Float_t extrapolposy22[ntrkmx];
  Float_t extrapolposz22[ntrkmx];

  Float_t extrapolposxerr22[ntrkmx];
  Float_t extrapolposyerr22[ntrkmx];
  Float_t extrapolposzerr22[ntrkmx];

  Float_t cmvhitrecoposx22[ntrkmx];
  Float_t cmvhitrecoposy22[ntrkmx];
  Float_t cmvhitrecoposz22[ntrkmx];

  Float_t cmvhittrueposx22[ntrkmx];
  Float_t cmvhittrueposy22[ntrkmx];
  Float_t cmvhittrueposz22[ntrkmx];

  Float_t cmvhitrecoposxerr22[ntrkmx];
  Float_t cmvhitrecoposyerr22[ntrkmx];
  Float_t cmvhitrecoposzerr22[ntrkmx];

  Float_t LeTime22[ntrkmx];
  Float_t RiTime22[ntrkmx];
  Float_t LePulse22[ntrkmx];
  Float_t RiPulse22[ntrkmx];

  Int_t clustersize30[ntrkmx];

  Float_t extrapolposx30[ntrkmx];
  Float_t extrapolposy30[ntrkmx];
  Float_t extrapolposz30[ntrkmx];

  Float_t extrapolposxerr30[ntrkmx];
  Float_t extrapolposyerr30[ntrkmx];
  Float_t extrapolposzerr30[ntrkmx];

  Float_t cmvhitrecoposx30[ntrkmx];
  Float_t cmvhitrecoposy30[ntrkmx];
  Float_t cmvhitrecoposz30[ntrkmx];

  Float_t cmvhittrueposx30[ntrkmx];
  Float_t cmvhittrueposy30[ntrkmx];
  Float_t cmvhittrueposz30[ntrkmx];

  Float_t cmvhitrecoposxerr30[ntrkmx];
  Float_t cmvhitrecoposyerr30[ntrkmx];
  Float_t cmvhitrecoposzerr30[ntrkmx];

  Float_t LeTime30[ntrkmx];
  Float_t RiTime30[ntrkmx];
  Float_t LePulse30[ntrkmx];
  Float_t RiPulse30[ntrkmx];



  Int_t clustersize31[ntrkmx];


  Float_t extrapolposx31[ntrkmx];
  Float_t extrapolposy31[ntrkmx];
  Float_t extrapolposz31[ntrkmx];

  Float_t extrapolposxerr31[ntrkmx];
  Float_t extrapolposyerr31[ntrkmx];
  Float_t extrapolposzerr31[ntrkmx];

  Float_t cmvhitrecoposx31[ntrkmx];
  Float_t cmvhitrecoposy31[ntrkmx];
  Float_t cmvhitrecoposz31[ntrkmx];

  Float_t cmvhittrueposx31[ntrkmx];
  Float_t cmvhittrueposy31[ntrkmx];
  Float_t cmvhittrueposz31[ntrkmx];

  Float_t cmvhitrecoposxerr31[ntrkmx];
  Float_t cmvhitrecoposyerr31[ntrkmx];
  Float_t cmvhitrecoposzerr31[ntrkmx];

  Float_t LeTime31[ntrkmx];
  Float_t RiTime31[ntrkmx];
  Float_t LePulse31[ntrkmx];
  Float_t RiPulse31[ntrkmx];


  Int_t clustersize32[ntrkmx];


  Float_t extrapolposx32[ntrkmx];
  Float_t extrapolposy32[ntrkmx];
  Float_t extrapolposz32[ntrkmx];

  Float_t extrapolposxerr32[ntrkmx];
  Float_t extrapolposyerr32[ntrkmx];
  Float_t extrapolposzerr32[ntrkmx];

  Float_t cmvhitrecoposx32[ntrkmx];
  Float_t cmvhitrecoposy32[ntrkmx];
  Float_t cmvhitrecoposz32[ntrkmx];

  Float_t cmvhittrueposx32[ntrkmx];
  Float_t cmvhittrueposy32[ntrkmx];
  Float_t cmvhittrueposz32[ntrkmx];

  Float_t cmvhitrecoposxerr32[ntrkmx];
  Float_t cmvhitrecoposyerr32[ntrkmx];
  Float_t cmvhitrecoposzerr32[ntrkmx];



  Float_t LeTime32[ntrkmx];
  Float_t RiTime32[ntrkmx];
  Float_t LePulse32[ntrkmx];
  Float_t RiPulse32[ntrkmx];


  //11022022

  Int_t clustersize40[ntrkmx];

  Float_t extrapolposx40[ntrkmx];
  Float_t extrapolposy40[ntrkmx];
  Float_t extrapolposz40[ntrkmx];

  Float_t extrapolposxerr40[ntrkmx];
  Float_t extrapolposyerr40[ntrkmx];
  Float_t extrapolposzerr40[ntrkmx];

  Float_t cmvhitrecoposx40[ntrkmx];
  Float_t cmvhitrecoposy40[ntrkmx];
  Float_t cmvhitrecoposz40[ntrkmx];

  Float_t cmvhittrueposx40[ntrkmx];
  Float_t cmvhittrueposy40[ntrkmx];
  Float_t cmvhittrueposz40[ntrkmx];

  Float_t cmvhitrecoposxerr40[ntrkmx];
  Float_t cmvhitrecoposyerr40[ntrkmx];
  Float_t cmvhitrecoposzerr40[ntrkmx];

  Float_t LeTime40[ntrkmx];
  Float_t RiTime40[ntrkmx];
  Float_t LePulse40[ntrkmx];
  Float_t RiPulse40[ntrkmx];



  Int_t clustersize41[ntrkmx];


  Float_t extrapolposx41[ntrkmx];
  Float_t extrapolposy41[ntrkmx];
  Float_t extrapolposz41[ntrkmx];

  Float_t extrapolposxerr41[ntrkmx];
  Float_t extrapolposyerr41[ntrkmx];
  Float_t extrapolposzerr41[ntrkmx];

  Float_t cmvhitrecoposx41[ntrkmx];
  Float_t cmvhitrecoposy41[ntrkmx];
  Float_t cmvhitrecoposz41[ntrkmx];

  Float_t cmvhittrueposx41[ntrkmx];
  Float_t cmvhittrueposy41[ntrkmx];
  Float_t cmvhittrueposz41[ntrkmx];

  Float_t cmvhitrecoposxerr41[ntrkmx];
  Float_t cmvhitrecoposyerr41[ntrkmx];
  Float_t cmvhitrecoposzerr41[ntrkmx];

  Float_t LeTime41[ntrkmx];
  Float_t RiTime41[ntrkmx];
  Float_t LePulse41[ntrkmx];
  Float_t RiPulse41[ntrkmx];


  Int_t clustersize42[ntrkmx];


  Float_t extrapolposx42[ntrkmx];
  Float_t extrapolposy42[ntrkmx];
  Float_t extrapolposz42[ntrkmx];

  Float_t extrapolposxerr42[ntrkmx];
  Float_t extrapolposyerr42[ntrkmx];
  Float_t extrapolposzerr42[ntrkmx];

  Float_t cmvhitrecoposx42[ntrkmx];
  Float_t cmvhitrecoposy42[ntrkmx];
  Float_t cmvhitrecoposz42[ntrkmx];

  Float_t cmvhittrueposx42[ntrkmx];
  Float_t cmvhittrueposy42[ntrkmx];
  Float_t cmvhittrueposz42[ntrkmx];

  Float_t cmvhitrecoposxerr42[ntrkmx];
  Float_t cmvhitrecoposyerr42[ntrkmx];
  Float_t cmvhitrecoposzerr42[ntrkmx];



  Float_t LeTime42[ntrkmx];
  Float_t RiTime42[ntrkmx];
  Float_t LePulse42[ntrkmx];
  Float_t RiPulse42[ntrkmx];


  Int_t clustersize50[ntrkmx];

  Float_t extrapolposx50[ntrkmx];
  Float_t extrapolposy50[ntrkmx];
  Float_t extrapolposz50[ntrkmx];

  Float_t extrapolposxerr50[ntrkmx];
  Float_t extrapolposyerr50[ntrkmx];
  Float_t extrapolposzerr50[ntrkmx];

  Float_t cmvhitrecoposx50[ntrkmx];
  Float_t cmvhitrecoposy50[ntrkmx];
  Float_t cmvhitrecoposz50[ntrkmx];

  Float_t cmvhittrueposx50[ntrkmx];
  Float_t cmvhittrueposy50[ntrkmx];
  Float_t cmvhittrueposz50[ntrkmx];

  Float_t cmvhitrecoposxerr50[ntrkmx];
  Float_t cmvhitrecoposyerr50[ntrkmx];
  Float_t cmvhitrecoposzerr50[ntrkmx];

  Float_t LeTime50[ntrkmx];
  Float_t RiTime50[ntrkmx];
  Float_t LePulse50[ntrkmx];
  Float_t RiPulse50[ntrkmx];



  Int_t clustersize51[ntrkmx];


  Float_t extrapolposx51[ntrkmx];
  Float_t extrapolposy51[ntrkmx];
  Float_t extrapolposz51[ntrkmx];

  Float_t extrapolposxerr51[ntrkmx];
  Float_t extrapolposyerr51[ntrkmx];
  Float_t extrapolposzerr51[ntrkmx];

  Float_t cmvhitrecoposx51[ntrkmx];
  Float_t cmvhitrecoposy51[ntrkmx];
  Float_t cmvhitrecoposz51[ntrkmx];

  Float_t cmvhittrueposx51[ntrkmx];
  Float_t cmvhittrueposy51[ntrkmx];
  Float_t cmvhittrueposz51[ntrkmx];

  Float_t cmvhitrecoposxerr51[ntrkmx];
  Float_t cmvhitrecoposyerr51[ntrkmx];
  Float_t cmvhitrecoposzerr51[ntrkmx];

  Float_t LeTime51[ntrkmx];
  Float_t RiTime51[ntrkmx];
  Float_t LePulse51[ntrkmx];
  Float_t RiPulse51[ntrkmx];


  Int_t clustersize52[ntrkmx];


  Float_t extrapolposx52[ntrkmx];
  Float_t extrapolposy52[ntrkmx];
  Float_t extrapolposz52[ntrkmx];

  Float_t extrapolposxerr52[ntrkmx];
  Float_t extrapolposyerr52[ntrkmx];
  Float_t extrapolposzerr52[ntrkmx];

  Float_t cmvhitrecoposx52[ntrkmx];
  Float_t cmvhitrecoposy52[ntrkmx];
  Float_t cmvhitrecoposz52[ntrkmx];

  Float_t cmvhittrueposx52[ntrkmx];
  Float_t cmvhittrueposy52[ntrkmx];
  Float_t cmvhittrueposz52[ntrkmx];

  Float_t cmvhitrecoposxerr52[ntrkmx];
  Float_t cmvhitrecoposyerr52[ntrkmx];
  Float_t cmvhitrecoposzerr52[ntrkmx];



  Float_t LeTime52[ntrkmx];
  Float_t RiTime52[ntrkmx];
  Float_t LePulse52[ntrkmx];
  Float_t RiPulse52[ntrkmx];


  Int_t clustersize60[ntrkmx];

  Float_t extrapolposx60[ntrkmx];
  Float_t extrapolposy60[ntrkmx];
  Float_t extrapolposz60[ntrkmx];

  Float_t extrapolposxerr60[ntrkmx];
  Float_t extrapolposyerr60[ntrkmx];
  Float_t extrapolposzerr60[ntrkmx];

  Float_t cmvhitrecoposx60[ntrkmx];
  Float_t cmvhitrecoposy60[ntrkmx];
  Float_t cmvhitrecoposz60[ntrkmx];

  Float_t cmvhittrueposx60[ntrkmx];
  Float_t cmvhittrueposy60[ntrkmx];
  Float_t cmvhittrueposz60[ntrkmx];

  Float_t cmvhitrecoposxerr60[ntrkmx];
  Float_t cmvhitrecoposyerr60[ntrkmx];
  Float_t cmvhitrecoposzerr60[ntrkmx];

  Float_t LeTime60[ntrkmx];
  Float_t RiTime60[ntrkmx];
  Float_t LePulse60[ntrkmx];
  Float_t RiPulse60[ntrkmx];



  Int_t clustersize61[ntrkmx];


  Float_t extrapolposx61[ntrkmx];
  Float_t extrapolposy61[ntrkmx];
  Float_t extrapolposz61[ntrkmx];

  Float_t extrapolposxerr61[ntrkmx];
  Float_t extrapolposyerr61[ntrkmx];
  Float_t extrapolposzerr61[ntrkmx];

  Float_t cmvhitrecoposx61[ntrkmx];
  Float_t cmvhitrecoposy61[ntrkmx];
  Float_t cmvhitrecoposz61[ntrkmx];

  Float_t cmvhittrueposx61[ntrkmx];
  Float_t cmvhittrueposy61[ntrkmx];
  Float_t cmvhittrueposz61[ntrkmx];

  Float_t cmvhitrecoposxerr61[ntrkmx];
  Float_t cmvhitrecoposyerr61[ntrkmx];
  Float_t cmvhitrecoposzerr61[ntrkmx];

  Float_t LeTime61[ntrkmx];
  Float_t RiTime61[ntrkmx];
  Float_t LePulse61[ntrkmx];
  Float_t RiPulse61[ntrkmx];


  Int_t clustersize62[ntrkmx];


  Float_t extrapolposx62[ntrkmx];
  Float_t extrapolposy62[ntrkmx];
  Float_t extrapolposz62[ntrkmx];

  Float_t extrapolposxerr62[ntrkmx];
  Float_t extrapolposyerr62[ntrkmx];
  Float_t extrapolposzerr62[ntrkmx];

  Float_t cmvhitrecoposx62[ntrkmx];
  Float_t cmvhitrecoposy62[ntrkmx];
  Float_t cmvhitrecoposz62[ntrkmx];

  Float_t cmvhittrueposx62[ntrkmx];
  Float_t cmvhittrueposy62[ntrkmx];
  Float_t cmvhittrueposz62[ntrkmx];

  Float_t cmvhitrecoposxerr62[ntrkmx];
  Float_t cmvhitrecoposyerr62[ntrkmx];
  Float_t cmvhitrecoposzerr62[ntrkmx];



  Float_t LeTime62[ntrkmx];
  Float_t RiTime62[ntrkmx];
  Float_t LePulse62[ntrkmx];
  Float_t RiPulse62[ntrkmx];



  //11

  Float_t extra_diff1[ntrkmx];
  Float_t extra_diff2[ntrkmx];
  Float_t extra_diff3[ntrkmx];



  Float_t ellip_diff00[ntrkmx];
  Float_t ellip_diff01[ntrkmx];
  Float_t ellip_diff02[ntrkmx];
  Float_t ellip_diff03[ntrkmx];

  Float_t ellip_diff10[ntrkmx];
  Float_t ellip_diff11[ntrkmx];
  Float_t ellip_diff12[ntrkmx];

  Float_t ellip_diff20[ntrkmx];
  Float_t ellip_diff21[ntrkmx];
  Float_t ellip_diff22[ntrkmx];

  Float_t ellip_diff30[ntrkmx];
  Float_t ellip_diff31[ntrkmx];
  Float_t ellip_diff32[ntrkmx];


  //cmv

  //simulation
  static const unsigned int cmv_nsimhtmx=1000;
  UInt_t cmv_nsimhit;
  UInt_t cmv_detid[cmv_nsimhtmx];
  Int_t  cmv_simpdgid[cmv_nsimhtmx];
  Float_t cmv_simtime[cmv_nsimhtmx];
  Float_t cmv_simenr[cmv_nsimhtmx];
  Float_t cmv_simposx[cmv_nsimhtmx];
  Float_t cmv_simposy[cmv_nsimhtmx];
  Float_t cmv_simposz[cmv_nsimhtmx];

  Float_t cmv_simpx[cmv_nsimhtmx];
  Float_t cmv_simpy[cmv_nsimhtmx];
  Float_t cmv_simpz[cmv_nsimhtmx];

  Float_t cmv_simlocx[cmv_nsimhtmx];
  Float_t cmv_simlocy[cmv_nsimhtmx];
  Float_t cmv_simlocz[cmv_nsimhtmx];
  //digitization
  // for data only one 16 bit integer + one 32bit integer arrays will be there
  static const unsigned int cmv_ndigihtmx=1000;

  UInt_t cmv_ndigihit;
  Int_t  cmv_digipdgid[cmv_ndigihtmx];
  UInt_t cmv_sipmid[cmv_ndigihtmx];
  // 3 bit for CMVD
  // 2 bit for Layer
  // 7 bit for Strip
  // 2 bit for SiPM
  Int_t cmv_digitimpul[cmv_ndigihtmx]; // 12bit for energy and 20 bit for time
  Int_t cmv_digitime[cmv_ndigihtmx];
  Int_t cmv_digipul[cmv_ndigihtmx];
  Float_t cmv_digiposx[cmv_ndigihtmx];
  Float_t cmv_digiposy[cmv_ndigihtmx];
  Float_t cmv_digiposz[cmv_ndigihtmx];

  Float_t cmv_digimom[cmv_ndigihtmx];
  Float_t cmv_digithe[cmv_ndigihtmx];
  Float_t cmv_digiphi[cmv_ndigihtmx];

  Float_t cmv_digilocx[cmv_ndigihtmx];
  Float_t cmv_digilocy[cmv_ndigihtmx];
  Float_t cmv_digilocz[cmv_ndigihtmx];
  //cmv sipm noise
  TFile* sipmnoise;
  TH1F* noise_hist[10][4];




private:
  micalDetectorParameterDef *paradef;
  ParameterMessenger *CardFile;

  double CorrTimeError;
  double UnCorrTimeError;
  double TimeToDigiConv;
  double SignalSpeed;

  int numberInLA;
  double parirlay[3];
  double parlay[3];

  double h2dDimX;
  double h2dDimY;
  int nbinxMag2d;
  int nbinyMag2d;
  double magZmax;
  int nbinxMagZ;


  int nUp;
  int nDown;
  //  int DetectorType;
  double PhotonSpeed;
  double CMVadctons;
};

#endif


// exp([0]*x+[1])+[2]
