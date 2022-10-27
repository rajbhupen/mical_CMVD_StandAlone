#include "MultiSimAnalysisDigi.hh"

MultiSimAnalysisDigi *MultiSimAnalysisDigi::AnPointer;

MultiSimAnalysisDigi::MultiSimAnalysisDigi() {
  AnPointer = this;
  paradef = micalDetectorParameterDef::AnPointer;
  CardFile = ParameterMessenger::AnPointer;
  //  DetectorType = CardFile->GetDetectorType();

 h2dDimX = paradef->GetParino(0);
  h2dDimY = paradef->GetParino(1);
   numberInLA = paradef->GetnLayer();
  cout<<"numberInLA "<<numberInLA<<endl;
  nbinxMag2d = (int)(h2dDimX/10);
  nbinyMag2d = (int)(h2dDimY/10);

  magZmax = 2*(paradef->GetParirlay(2) + paradef->GetParlay(2));
  nbinxMagZ = (int)(magZmax);
  CardFile->PrintParameters();
  isInOut = CardFile->GetInputOutput();
  isVisOut = CardFile->GetVisualOutput();
  isXtermOut = CardFile->GetXTermOutput();
  collatedIn = CardFile->GetCollatedIn();
  PhotonSpeed = CardFile->GetPhotonSpeed();
  CMVadctons = CardFile->GetCMVadctons();

  nmxhit	= 10000;
  ievt		= 0;
  ihist		= 0;
  ievt_wt=0;
  intxn_id=0;

  CorrTimeError = CardFile->GetCorrTimeSmr();
  UnCorrTimeError = CardFile->GetUnCorrTimeSmr();
  TimeToDigiConv = CardFile->GetTimeToDigiConv();
  // cout<<"Hellllllllllll"<<TimeToDigiConv<<endl;
  SignalSpeed = CardFile->GetSignalSpeed();

  H	=0;
  Hp	=0;

  pRootFile	=0;
  inputRootFile=0;
  pVisFile = 0;
  pEventTree	=0;
  inputEventTree=0;
  visTree = 0;

  // double pPosLimXY = 52.5;
  // double pPosLimZ = 2.5;

  double pPosLimXY = 105;
  double pPosLimZ = 20.5;
  cout<<"checck2"<<endl;
  pPosX = new TH1F("deltaX", "#Delta X (cm)", 105, -pPosLimXY, pPosLimXY);
  pPosY = new TH1F("deltaY", "#Delta Y (cm)", 105, -pPosLimXY, pPosLimXY);
  pPosZ = new TH1F("deltaZ", "#Delta Z (cm)", 100, -pPosLimZ, pPosLimZ);

  char title[100];
  for (int ij=0; ij<20; ij++) {
    sprintf(title, "dedz_%i", ij);
      pdedz[ij] = new TH1F(title, title, 4000, -600., 600.);
  }

  pPosXX = new TH2F("deltaXX", "#Delta XX (cm)", 100, -2500, 2500, 100, -pPosLimXY, pPosLimXY);
  pPosYY = new TH2F("deltaYY", "#Delta YY (cm)", 100, -500, 500, 100, -pPosLimXY, pPosLimXY);
  pPosZZ = new TH2F("deltaZZ", "#Delta ZZ (cm)", 100, -1000, 1000, 100, -pPosLimZ, pPosLimZ);

  hitDist  = new TH1F("HitDist","This is the total rawhit distribution across layers",182,0,182);     //asm
  TrkDist  = new TH1F("TrkDist","This is the total track-hit distribution across layers",182,0,182);     //asm
  EffDist  = new TH1F("EffDist","This is the efficiency distribution across layers",182,0,182);     //asm

  InoTrack_listsize =new TH1F("TrksizeDist","Distribution of tracks per event ",10,0,10); //asm

  ShwXw	= new TH1F("ShwXw ","This is a distribution for X position resolution in shower region " ,80 , -3.14 , 3.14);  //asm
  ShwYw	= new TH1F("ShwYw ","This is a distribution for Y position resolution in shower region " ,80 , -3.14 , 3.14);  //asm
  RC		= new TH2F("RCvsff","This is a distribution for raduis of curvature at the edges " ,1000 , 0 , 1000,100,-50,50);  //asm
  DGap	= new TH1F("DetGap"," Plot to see the performance of code at detgap ",10,0,10);  //asm
  //  nhtcal1 = 0;
  //  nhtcal2 = 0;

  cout <<"booooooooooooooooooooooooooooooooooooooooooooooked "<<endl;

  DeadStripX = new TH1F("DeadStripX","Plot to see dead strips on X-plane",64,-0.5,63.5);
  NoisyStripX = new TH1F("NoisyStripX","Plot to see noisy strips on X-Plane",64,-0.5,63.5);
  DeadStripY = new TH1F("DeadStripY","Plot to see dead strips on Y-Plane",64,-0.5,63.5);
  NoisyStripY = new TH1F("NoisyStripY","Plot to see noisy strips on Y-Plane",64,-0.5,63.5);
  DiffTime = new TH1F("DiffTime","DiffTime",100,-100.,100.);

  strpXtime = new TH1F("strpXtime","xsmrtime - xtruetime",100,-40,40.);
  strpYtime = new TH1F("strpYtime","ysmrtime - ytruetime",100,-40,40.);
  strpXtimeCorr = new TH1F("strpXtimeCorr","xsmrtimeCorr - xtruetime",100,-20,20.);
  strpYtimeCorr = new TH1F("strpYtimeCorr","ysmrtimeCorr - ytruetime",100,-20,20.);
  hitXtime = new TH1F("hitXtime","hittime - xtruetime",100,-20,20.);
  hitYtime = new TH1F("hitYtime","hittime - ytruetime",100,-20,20.);

  smagFieldX = new TH1D("smagfieldinx","magfieldinx simu",nbinxMagZ,0,magZmax);
  smagFieldY = new TH1D("smagfieldiny","magfieldiny simu",nbinxMagZ,0,magZmax);
  smag2dX = new TH2D("smag2dinx","magfieldinx 2D simu",nbinxMagZ,0,magZmax,180,-2.0,2.0);
  smag2dY = new TH2D("smag2diny","magfieldiny 2D simu",nbinxMagZ,0,magZmax,180,-2.0,2.0);

  rmagFieldX = new TH1D("rmagfieldinx","magfieldinx in reco",nbinxMagZ,0,magZmax);
  rmagFieldY = new TH1D("rmagfieldiny","magfieldiny in reco",nbinxMagZ,0,magZmax);
  rmag2dX = new TH2D("rmag2dinx","magfieldinx 2D in reco",nbinxMagZ,0,magZmax,180,-2.0,2.0);
  rmag2dY = new TH2D("rmag2diny","magfieldiny 2D in reco",nbinxMagZ,0,magZmax,180,-2.0,2.0);
  // int nbin2d = 1600;
  smag2dXYpixel_iron = new TH2D("smag2dxypixeliron","magfield in xy pixel in sim in iron",nbinxMag2d,-h2dDimX,h2dDimX,nbinyMag2d,-h2dDimY,h2dDimY);
  smag2dXYpixel_air = new TH2D("smag2dxypixelair","magfield in xy pixel in sim in air",nbinxMag2d,-h2dDimX,h2dDimX,nbinyMag2d,-h2dDimY,h2dDimY);
  rmag2dXYpixel_iron = new TH2D("rmag2dxypixeliron","magfield in xy pixel in reco in iron",nbinxMag2d,-h2dDimX,h2dDimX,nbinyMag2d,-h2dDimY,h2dDimY);
  rmag2dXYpixel_air = new TH2D("rmag2dxypixelair","magfield in xy pixel in reco in air",nbinxMag2d,-h2dDimX,h2dDimX,nbinyMag2d,-h2dDimY,h2dDimY);

  xyvsbxin = new TH2D("xyvsbxin", "xyvsbxin", 321, -8025, 8025, 321, -8025, 8025);
  xyvsbyin = new TH2D("xyvsbyin", "xyvsbyin", 321, -8025, 8025, 321, -8025, 8025);
  xyvsbxdiff = new TH2D("xyvsbxdiff", "xyvsbxdiff", 321, -8025, 8025, 321, -8025, 8025);
  xyvsbydiff = new TH2D("xyvsbydiff", "xyvsbydiff", 321, -8025, 8025, 321, -8025, 8025);
  // xyvsbxindiff = new TH2D("xyvsbxindiff", "xyvsbxindiff", 319, -7975, 7975, 319, -7975, 7975);
  // xyvsbyindiff = new TH2D("xyvsbyindiff", "xyvsbyindiff", 319, -7975, 7975, 319, -7975, 7975);

   xyvsbxout = new TH2D("xyvsbxout", "xyvsbxout", 969, -24225, 24225, 321, -8025, 8025);
   xyvsbyout =  new TH2D("xyvsbyout", "xyvsbyout", 969, -24225, 24225, 321, -8025, 8025);

  //G4int hit_wo_ghst;//SSE
  //hprof  = new TProfile("hprof","Profile",2,0,100,0.,1000.);//SSE
 // hh_E 	= new TH2D("hh_E","E_had vs orighits;orighits;E_had (GeV)",25, 0. ,50.,100, 0., 100.);  //SSE
  //hh_woghst_E 	= new TH2D("hh_woghst_E","E_had vs hits;hits;E_had (GeV)",50, 0. ,500.,100, 0., 100.);  //SSE
  //hist_nhits_LargestCluster_E 	= new TH2D("hist_nhits_LargestCluster_E","E_had vs no of hits in largest cluster;hits;E_had (GeV)",50, 0. ,500.,100, 0., 100.);  //SSE
  //hist_orighits_new_E= new TH2D("hist_orighits_new_E","E_had vs orighits (new) ;orighits_new;E_had (GeV)",50, 0. ,100.,100, 0., 100.);//SSE 09/15
 //hist_orighits_mod_E= new TH2D("hist_orighits_mod_E","E_had vs orighits (modified) ;orighits_mod;E_had (GeV)",25, 0. ,50.,100, 0., 100.);//SSE 09/15
 //hist_wogh_orighits_E= new TH2D("hist_wogh_orighits_E","E_had vs wogh_orighits ;Ghosthit_removed_orighit;E_had (GeV)",25, 0. ,50.,100, 0., 100.);//SSE 09/15
   cout<<"checck2"<<endl;


}

MultiSimAnalysisDigi::~MultiSimAnalysisDigi() {
  cout<<"Deleting MultiSimAnalysisDigi class pointer ..."<<endl;
}

void MultiSimAnalysisDigi::OpenInputRootFiles(char* infile) {

  char inRootFile[400];

  cout << " opening root file " << infile << " inout " <<  isInOut << endl;

  if(isInOut==0 ||isInOut==3 || isInOut==1 ) { ////inout3
    //Input is sim file from MC for sim to reco(3) or sim to digi(0)

    if(isInOut==3 || isInOut==0){
    sprintf(inRootFile,"%s",infile);
    inputRootFile = new TFile(inRootFile, "read");
    cout<< "Data is read from simulation output rootfile : "<< inRootFile <<endl;
    inputEventTree = (TTree*)inputRootFile->Get("T3"); //T1
    }
    else if(isInOut==1){
   sprintf(inRootFile,"%s",infile);
    inputRootFile = new TFile(inRootFile, "read");
    cout<< "Data is read from digitization output file : "<< inRootFile  <<endl;
    inputEventTree = (TTree*)inputRootFile->Get("T1"); //T2
    }

    
    inputEventTree->SetBranchAddress("irun",&irun);
    inputEventTree->SetBranchAddress("ievt",&ievt);
    inputEventTree->SetBranchAddress("ngent",&ngent);
    inputEventTree->SetBranchAddress("pidin",pidin);
    inputEventTree->SetBranchAddress("momin",momin);
    inputEventTree->SetBranchAddress("thein",thein);
    inputEventTree->SetBranchAddress("phiin",phiin);
    inputEventTree->SetBranchAddress("posxin",posxin);
    inputEventTree->SetBranchAddress("posyin",posyin);
    inputEventTree->SetBranchAddress("poszin",poszin);
    inputEventTree->SetBranchAddress("nsimht", &nsimht);
    inputEventTree->SetBranchAddress("detid", detid);
    inputEventTree->SetBranchAddress("simpdgid", simpdgid);
    inputEventTree->SetBranchAddress("simtime", simtime);
    inputEventTree->SetBranchAddress("simenr", simenr);
    inputEventTree->SetBranchAddress("simvx", simvx);
    inputEventTree->SetBranchAddress("simvy", simvy);
    inputEventTree->SetBranchAddress("simvz", simvz);
    inputEventTree->SetBranchAddress("simpx", simpx);
    inputEventTree->SetBranchAddress("simpy", simpy);
    inputEventTree->SetBranchAddress("simpz", simpz);
    inputEventTree->SetBranchAddress("simlocvx", simlocvx);
    inputEventTree->SetBranchAddress("simlocvy", simlocvy);
    inputEventTree->SetBranchAddress("ngenerated", &ngenerated);
    inputEventTree->SetBranchAddress("naperture", &naperture);
    inputEventTree->SetBranchAddress("trigx", &trigx);
    inputEventTree->SetBranchAddress("trigy", &trigy);

    if(CardFile->GetCMVD()==1){
  //cmv
    inputEventTree->SetBranchAddress("cmv_nsimhit", &cmv_nsimhit);
     inputEventTree->SetBranchAddress("cmv_detid", cmv_detid) ;
    inputEventTree->SetBranchAddress("cmv_simpdgid", cmv_simpdgid);
    inputEventTree->SetBranchAddress("cmv_simtime", cmv_simtime);
    inputEventTree->SetBranchAddress("cmv_simenr", cmv_simenr);
    inputEventTree->SetBranchAddress("cmv_simposx", cmv_simposx);
    inputEventTree->SetBranchAddress("cmv_simposy", cmv_simposy);
    inputEventTree->SetBranchAddress("cmv_simposz", cmv_simposz);
    inputEventTree->SetBranchAddress("cmv_simpx", cmv_simpx);
    inputEventTree->SetBranchAddress("cmv_simpy", cmv_simpy);
    inputEventTree->SetBranchAddress("cmv_simpz", cmv_simpz);
    inputEventTree->SetBranchAddress("cmv_simlocx", cmv_simlocx);
    inputEventTree->SetBranchAddress("cmv_simlocy", cmv_simlocy);
 inputEventTree->SetBranchAddress("cmv_simlocz", cmv_simlocz);
    }// if(CardFile->GetCMVD()==1){

    if(isInOut==1) {
    inputEventTree->SetBranchAddress("ndigiht", &ndigiht);
    inputEventTree->SetBranchAddress("stripid", stripid);
    inputEventTree->SetBranchAddress("digipdgid", digipdgid);
    inputEventTree->SetBranchAddress("digitime", digitime);
    inputEventTree->SetBranchAddress("digitruetime", digitruetime);
    inputEventTree->SetBranchAddress("digienr", digienr);
    inputEventTree->SetBranchAddress("digivx", digivx);
    inputEventTree->SetBranchAddress("digivy", digivy);
    inputEventTree->SetBranchAddress("digivz", digivz);
    inputEventTree->SetBranchAddress("digipx", digipx);
    inputEventTree->SetBranchAddress("digipy", digipy);
    inputEventTree->SetBranchAddress("digipz", digipz);
    // inputEventTree->SetBranchAddress("diginoise", diginoise);
       //cmv
     if(CardFile->GetCMVD()==1){
    inputEventTree->SetBranchAddress("cmv_ndigihit", &cmv_ndigihit);
    inputEventTree->SetBranchAddress("cmv_sipmid", cmv_sipmid);
    inputEventTree->SetBranchAddress("cmv_digipdgid", cmv_digipdgid);
    inputEventTree->SetBranchAddress("cmv_digitimpul", cmv_digitimpul);
    inputEventTree->SetBranchAddress("cmv_digitime", cmv_digitime);
    inputEventTree->SetBranchAddress("cmv_digipul", cmv_digipul);
    inputEventTree->SetBranchAddress("cmv_digiposx", cmv_digiposx);
    inputEventTree->SetBranchAddress("cmv_digiposy", cmv_digiposy);
    inputEventTree->SetBranchAddress("cmv_digiposz", cmv_digiposz);
    inputEventTree->SetBranchAddress("cmv_digimom", cmv_digimom);
    inputEventTree->SetBranchAddress("cmv_digithe", cmv_digithe);
    inputEventTree->SetBranchAddress("cmv_digiphi", cmv_digiphi);
    inputEventTree->SetBranchAddress("cmv_digilocx", cmv_digilocx);
    inputEventTree->SetBranchAddress("cmv_digilocy", cmv_digilocy);
    inputEventTree->SetBranchAddress("cmv_digilocz", cmv_digilocz);
    //cmv
     }// if(CardFile->GetCMVD()==1){

    }
  } else if(isInOut==2) { // datareading
    sprintf(inRootFile,"%s",infile);
    inputRootFile = new TFile(inRootFile, "read");
    cout<< "Data is read from data file : "<< inRootFile  <<endl;

    inputEventTree = (TTree*)inputRootFile->Get("evetree");
    data_event = new RPCEve(inputEventTree);
    // data_event = new evetree(inputEventTree);
    data_event->Loop();

  }

}

void MultiSimAnalysisDigi::OpenOutputRootFiles(char* outfile) {

  cout << " OpenOutputRootFiles " << outfile << " inout " <<  isInOut << endl;

  char outRootFile[400];
  sprintf(outRootFile,"%s.root",outfile);


    pRootFile = new TFile(outRootFile, "RECREATE");
    if (!pRootFile) {
      cout << "Error opening output root file !" << endl;
      exit(-1);
    } else {
      cout<< "Output stored in root file: "<< outRootFile <<endl;
    }

    if (isInOut==0 ){
      pEventTree = new TTree("T2", "INODIGI");
    }else if (isInOut==1||isInOut==3 ){
      pEventTree = new TTree("T1", "INORECO");
    } else if(isInOut==2) {
      pEventTree = new TTree("T4", "INODATA");
 }

    //These are common to all:

    pEventTree->Branch("irun",&irun,"irun/i"); //VALGRIND
    pEventTree->Branch("ievt",&ievt2,"ievt/i");
    pEventTree->Branch("ngent",&ngent,"ngent/i");
    pEventTree->Branch("pidin",pidin,"pidin[ngent]/I");
    pEventTree->Branch("ievt_wt",&ievt_wt,"ievt_wt/F");
    pEventTree->Branch("intxn_id",&intxn_id,"intxn_id/I");
    pEventTree->Branch("momin",momin,"momin[ngent]/F");
    pEventTree->Branch("thein",thein,"thein[ngent]/F");
    pEventTree->Branch("phiin",phiin,"phiin[ngent]/F");
    pEventTree->Branch("posxin",posxin,"posxin[ngent]/F");
    pEventTree->Branch("posyin",posyin,"posyin[ngent]/F");
    pEventTree->Branch("poszin",poszin,"poszin[ngent]/F");
    pEventTree->Branch("ngenerated",&ngenerated,"ngenerated/i");
    pEventTree->Branch("naperture",&naperture,"naperture/i");


    //DigiOutput
     if (isInOut==0 || isInOut==3 || isInOut==1 ){


       pEventTree->Branch("nsimht", &nsimht, "nsimht/i");
       pEventTree->Branch("detid", detid, "detid[nsimht]/i");
       pEventTree->Branch("simpdgid", simpdgid, "simpdgid[nsimht]/I");
       pEventTree->Branch("simtime", simtime, "simtime[nsimht]/F");
       pEventTree->Branch("simenr", simenr, "simenr[nsimht]/F");
       pEventTree->Branch("simvx", simvx, "simvx[nsimht]/F");
       pEventTree->Branch("simvy", simvy, "simvy[nsimht]/F");
       pEventTree->Branch("simvz", simvz, "simvz[nsimht]/F");
       pEventTree->Branch("simpx", simpx, "simpx[nsimht]/F");
       pEventTree->Branch("simpy", simpy, "simpy[nsimht]/F");
       pEventTree->Branch("simpz", simpz, "simpz[nsimht]/F");
       pEventTree->Branch("simlocvx", simlocvx, "simlocvx[nsimht]/F");
       pEventTree->Branch("simlocvy", simlocvy, "simlocvy[nsimht]/F");
       pEventTree->Branch("simlocvz", simlocvz, "simlocvz[nsimht]/F");

       
       
    pEventTree->Branch("ndigiht", &ndigiht, "ndigiht/i");
    pEventTree->Branch("trigx", &trigx, "trigx/I");
    pEventTree->Branch("trigy", &trigy, "trigy/I");
    pEventTree->Branch("triggeracceptance",&triggeracceptance,"triggeracceptance/i");
    pEventTree->Branch("stripid", stripid, "stripid[ndigiht]/i");
    pEventTree->Branch("digipdgid", digipdgid, "digipdgid[ndigiht]/I");
    pEventTree->Branch("digitime", digitime, "digitime[ndigiht]/I");
    pEventTree->Branch("digitruetime", digitruetime, "digitruetime[ndigiht]/I");
    pEventTree->Branch("digienr", digienr, "digienr[ndigiht]/F");
    pEventTree->Branch("digivx", digivx, "digivx[ndigiht]/F");
    pEventTree->Branch("digivy", digivy, "digivy[ndigiht]/F");
    pEventTree->Branch("digivz", digivz, "digivz[ndigiht]/F");
    pEventTree->Branch("digipx", digipx, "digipx[ndigiht]/F");
    pEventTree->Branch("digipy", digipy, "digipy[ndigiht]/F");
    pEventTree->Branch("digipz", digipz, "digipz[ndigiht]/F");

       cout<<"check4"<<endl;
     if(CardFile->GetCMVD()==1){

       pEventTree->Branch("cmv_nsimhit", &cmv_nsimhit, "cmv_nsimhit/i");
       pEventTree->Branch("cmv_detid", cmv_detid, "cmv_detid[cmv_nsimhit]/I") ;
       pEventTree->Branch("cmv_simpdgid", cmv_simpdgid, "cmv_simpdgid[cmv_nsimhit]/I");
       pEventTree->Branch("cmv_simtime", cmv_simtime, "cmv_simtime[cmv_nsimhit]/F");
       pEventTree->Branch("cmv_simenr", cmv_simenr, "cmv_simenr[cmv_nsimhit]/F");
       pEventTree->Branch("cmv_simposx", cmv_simposx, "cmv_simposx[cmv_nsimhit]/F");
       pEventTree->Branch("cmv_simposy", cmv_simposy, "cmv_simposy[cmv_nsimhit]/F");
       pEventTree->Branch("cmv_simposz", cmv_simposz, "cmv_simposz[cmv_nsimhit]/F");
       
       pEventTree->Branch("cmv_simpx", cmv_simpx, "cmv_simpx[cmv_nsimhit]/F");
       pEventTree->Branch("cmv_simpy", cmv_simpy, "cmv_simpy[cmv_nsimhit]/F");
       pEventTree->Branch("cmv_simpz", cmv_simpz, "cmv_simpz[cmv_nsimhit]/F");
       
     pEventTree->Branch("cmv_simlocx", cmv_simlocx, "cmv_simlocx[cmv_nsimhit]/F");
     pEventTree->Branch("cmv_simlocy", cmv_simlocy, "cmv_simlocy[cmv_nsimhit]/F");
     pEventTree->Branch("cmv_simlocz", cmv_simlocz, "cmv_simlocz[cmv_nsimhit]/F");
     
       //cmv
    pEventTree->Branch("cmv_ndigihit", &cmv_ndigihit, "cmv_ndigihit/i");
    pEventTree->Branch("cmv_sipmid", cmv_sipmid, "cmv_sipmid[cmv_ndigihit]/i") ;
    pEventTree->Branch("cmv_digipdgid", cmv_digipdgid, "cmv_digipdgid[cmv_ndigihit]/I");
    pEventTree->Branch("cmv_digitimpul", cmv_digitimpul, "cmv_digitimpul[cmv_ndigihit]/i");
   pEventTree->Branch("cmv_digitime", cmv_digitime, "cmv_digitime[cmv_ndigihit]/i");
    pEventTree->Branch("cmv_digipul", cmv_digipul, "cmv_digipul[cmv_ndigihit]/i");
    pEventTree->Branch("cmv_digiposx", cmv_digiposx, "cmv_digiposx[cmv_ndigihit]/F");
    pEventTree->Branch("cmv_digiposy", cmv_digiposy, "cmv_digiposy[cmv_ndigihit]/F");
    pEventTree->Branch("cmv_digiposz", cmv_digiposz, "cmv_digiposz[cmv_ndigihit]/F");

    pEventTree->Branch("cmv_digimom", cmv_digimom, "cmv_digimom[cmv_ndigihit]/F");
    pEventTree->Branch("cmv_digithe", cmv_digithe, "cmv_digithe[cmv_ndigihit]/F");
    pEventTree->Branch("cmv_digiphi", cmv_digiphi, "cmv_digiphi[cmv_ndigihit]/F");

    pEventTree->Branch("cmv_digilocx", cmv_digilocx, "cmv_digilocx[cmv_ndigihit]/F");
    pEventTree->Branch("cmv_digilocy", cmv_digilocy, "cmv_digilocy[cmv_ndigihit]/F");
    pEventTree->Branch("cmv_digilocz", cmv_digilocz, "cmv_digilocz[cmv_ndigihit]/F");
     }// if(CardFile->GetCMVD()==1){

     }//     if (isInOut==0 || isInOut==3  ){


        if (isInOut){

    pEventTree->Branch("hw_trig",&hw_trig,"hw_trig/I");
    // pEventTree->Branch("hw_trigy",&hw_trigy,"hw_trigy/I");
    pEventTree->Branch("sw_trigx",&sw_trigx,"sw_trigx/I");
    pEventTree->Branch("sw_trigy",&sw_trigy,"sw_trigy/I");
    pEventTree->Branch("ntrkt",&ntrkt,"ntrkt/i");
    pEventTree->Branch("itype",itype,"itype[ntrkt]/I");
    pEventTree->Branch("nLayer",&nLayer,"nLayer/I");
    pEventTree->Branch("nhits",nhits,"nhits[ntrkt]/I");
    pEventTree->Branch("nhits_finder",nhits_finder,"nhits_finder[ntrkt]/I");
    pEventTree->Branch("chisq",chisq,"chisq[ntrkt]/F");
    pEventTree->Branch("chisq2",chisq2,"chisq2[ntrkt]/F");
    pEventTree->Branch("ndof",ndof,"ndof[ntrkt]/i");
    pEventTree->Branch("ndof2",ndof2,"ndof2[ntrkt]/i");
    pEventTree->Branch("cvalue",cvalue,"cvalue[ntrkt]/F");
    pEventTree->Branch("trkmm",trkmm,"trkmm[ntrkt]/F");
    pEventTree->Branch("trkth",trkth,"trkth[ntrkt]/F");
    pEventTree->Branch("trkph",trkph,"trkph[ntrkt]/F");
    pEventTree->Branch("momvx",momvx,"momvx[ntrkt]/F");
    pEventTree->Branch("thevx",thevx,"thevx[ntrkt]/F");
    pEventTree->Branch("phivx",phivx,"phivx[ntrkt]/F");
    pEventTree->Branch("posxvx",posxvx,"posxvx[ntrkt]/F");
    pEventTree->Branch("posyvx",posyvx,"posyvx[ntrkt]/F");
    pEventTree->Branch("poszvx",poszvx,"poszvx[ntrkt]/F");

    pEventTree->Branch("momrf",momrf,"momrf[ntrkt]/F");
    pEventTree->Branch("therf",therf,"therf[ntrkt]/F");
    pEventTree->Branch("phirf",phirf,"phirf[ntrkt]/F");
    pEventTree->Branch("posxrf",posxrf,"posxrf[ntrkt]/F");
    pEventTree->Branch("posyrf",posyrf,"posyrf[ntrkt]/F");
    pEventTree->Branch("poszrf",poszrf,"poszrf[ntrkt]/F");



    //added for straightline fit:raj
    pEventTree->Branch("XPosdev_exclu",XPosdev_exclu,"XPosdev_exclu[nLayer]/F");
    pEventTree->Branch("YPosdev_exclu",YPosdev_exclu,"YPosdev_exclu[nLayer]/F");

    pEventTree->Branch("XPosdev",XPosdev,"XPosdev[10]/F");
    pEventTree->Branch("YPosdev",YPosdev,"YPosdev[10]/F");



    pEventTree->Branch("momend",momend,"momend[ntrkt]/F");
    pEventTree->Branch("theend",theend,"theend[ntrkt]/F");
    pEventTree->Branch("phiend",phiend,"phiend[ntrkt]/F");
    pEventTree->Branch("posxend",posxend,"posxend[ntrkt]/F");
    pEventTree->Branch("posyend",posyend,"posyend[ntrkt]/F");
    pEventTree->Branch("strpxend",strpxend,"strpxend[ntrkt]/I");
    pEventTree->Branch("strpyend",strpyend,"strpyend[ntrkt]/I");
    pEventTree->Branch("poszend",poszend,"poszend[ntrkt]/F");
    pEventTree->Branch("momds"   ,momds   ,"momds[ntrkt]/F");
    pEventTree->Branch("momrg"   ,momrg   ,"momrg[ntrkt]/F");
    pEventTree->Branch("vtxzplane",vtxzplane,"vtxzplane[ntrkt]/I");
    pEventTree->Branch("endzplane",endzplane,"endzplane[ntrkt]/I");
    pEventTree->Branch("nvisht", &nvisht, "nvisht/i");
    pEventTree->Branch("fitposxx",fitposxx,"fitposxx[nvisht]/F");
    pEventTree->Branch("fitposyy",fitposyy,"fitposyy[nvisht]/F");
    pEventTree->Branch("fitposzz",fitposzz,"fitposzz[nvisht]/F");
    pEventTree->Branch("fitlayzz",fitlayzz,"fitlayzz[nvisht]/F");
    // pEventTree->Branch("fitlayx2",fitlayx2,"fitlayx2[nvisht]/F");
    // pEventTree->Branch("fitlayx3",fitlayx3,"fitlayx3[nvisht]/F");
    // pEventTree->Branch("fitlayx4",fitlayx4,"fitlayx4[nvisht]/F");
    pEventTree->Branch("fitlaymom",fitlaymom,"fitlaymom[nvisht]/F");
    pEventTree->Branch("extrapolxx",extrapolxx,"extrapolxx[nvisht]/F");
    pEventTree->Branch("extrapolyy",extrapolyy,"extrapolyy[nvisht]/F");
    pEventTree->Branch("extrapolmom",extrapolmom,"extrapolmom[nvisht]/F");
    pEventTree->Branch("momdiff1",&momdiff1,"momdiff1/F");
    pEventTree->Branch("radialdiff1",&radialdiff1,"radialdiff1/F");
    // pEventTree->Branch("fitlaythe",fitlaythe,"fitlaythe[nvisht]/F");
    // pEventTree->Branch("fitlayphi",fitlayphi,"fitlayphi[nvisht]/F");
    pEventTree->Branch("nvisclst", &nvisclst, "nvisclst/i");
    pEventTree->Branch("clstposxx",clstposxx,"clstposxx[nvisclst]/F");
    pEventTree->Branch("clstposyy",clstposyy,"clstposyy[nvisclst]/F");
    pEventTree->Branch("clstposzz",clstposzz,"clstposzz[nvisclst]/F");
    pEventTree->Branch("clstposzpln",clstposzpln,"clstposzpln[nvisclst]/I");

    pEventTree->Branch("strtnhitsx",&strtnhitsx,"strtnhitsx/I");
    pEventTree->Branch("strtchisqx",&strtchisqx,"strtchisqx/F");
    pEventTree->Branch("strtintercptx",&strtintercptx,"strtintercptx/F");
    pEventTree->Branch("strtslopex",&strtslopex,"strtslopex/F");
    pEventTree->Branch("strtnhitsy",&strtnhitsy,"strtnhitsy/I");
    pEventTree->Branch("strtchisqy",&strtchisqy,"strtchisqy/F");
    pEventTree->Branch("strtintercpty",&strtintercpty,"strtintercpty/F");
    pEventTree->Branch("strtslopey",&strtslopey,"strtslopey/F");

    pEventTree->Branch("simpleradii",&simpleradii,"simpleradii/F");
    pEventTree->Branch("simplecurv",&simplecurv,"simplecurv/F");
    pEventTree->Branch("simplex0",&simplex0,"simplex0/F");
    pEventTree->Branch("simplez0",&simplez0,"simplez0/F");
    pEventTree->Branch("simplechisqpos",&simplechisqpos,"simplechisqpos/F");
    pEventTree->Branch("simplechisqneg",&simplechisqneg,"simplechisqneg/F");
    pEventTree->Branch("simplechisqcndn",&simplechisqcndn,"simplechisqcndn/F");
    pEventTree->Branch("simpleavgxpos",&simpleavgxpos,"simpleavgxpos/F");
    pEventTree->Branch("simpleavgxneg",&simpleavgxneg,"simpleavgxneg/F");
    pEventTree->Branch("simpleavgxcndn",&simpleavgxcndn,"simpleavgxcndn/F");
    pEventTree->Branch("simpleavgxmeas",&simpleavgxmeas,"simpleavgxmeas/F");
    pEventTree->Branch("simplenhits",&simplenhits,"simplenhits/I");
    pEventTree->Branch("nhits_last",&nhits_last,"nhits_last/I");
    pEventTree->Branch("nhits_last_m1",&nhits_last_m1,"nhits_last_m1/I");
    pEventTree->Branch("nhits_below",&nhits_below,"nhits_below/I");
    pEventTree->Branch("ftime_last",&ftime_last,"ftime_last/F");

    //Extrapolated and measured state vector informations JIM
    pEventTree->Branch("inPosX", inPosX, "inPosX[10]/F");
    pEventTree->Branch("inPosY", inPosY, "inPosY[10]/F");
    pEventTree->Branch("extPosX", extPosX, "extPosX[10]/F");
    pEventTree->Branch("extPosY", extPosY, "extPosY[10]/F");
    pEventTree->Branch("nXStrips", nXStrips, "nXStrips[10]/I");
    pEventTree->Branch("nYStrips", nYStrips, "nYStrips[10]/I");

    pEventTree->Branch("inPosfX", inPosfX, "inPosfX[10]/F");
    pEventTree->Branch("inPosfY", inPosfY, "inPosfY[10]/F");
    pEventTree->Branch("extPosfX", extPosfX, "extPosfX[10]/F");
    pEventTree->Branch("extPosfY", extPosfY, "extPosfY[10]/F");
    pEventTree->Branch("nfXStrips", nfXStrips, "nfXStrips[10]/I");
    pEventTree->Branch("nfYStrips", nfYStrips, "nfYStrips[10]/I");





    
    pEventTree->Branch("XchisqOccu", XchisqOccu, "XchisqOccu[10]/F");
    pEventTree->Branch("YchisqOccu", YchisqOccu, "YchisqOccu[10]/F");
    pEventTree->Branch("XndfOccu", XndfOccu, "XndfOccu[10]/I");
    pEventTree->Branch("YndfOccu", YndfOccu, "YndfOccu[10]/I");
    pEventTree->Branch("nxFailOccu", nxFailOccu, "nxFailOccu[10]/I");
    pEventTree->Branch("nyFailOccu", nyFailOccu, "nyFailOccu[10]/I");

    //cmvd

 if(CardFile->GetCMVD()==1){



  pEventTree->Branch("CMVDnLayer", &CMVDnLayer, "CMVDnLayer/i");

  pEventTree->Branch("CMVDTruePosX", CMVDTruePosX, "CMVDTruePosX[32]/F");
  pEventTree->Branch("CMVDTruePosY", CMVDTruePosY, "CMVDTruePosY[32]/F");
  pEventTree->Branch("CMVDTruePosZ", CMVDTruePosZ, "CMVDTruePosZ[32]/F");

  pEventTree->Branch("CMVDRecoPosX", CMVDRecoPosX, "CMVDRecoPosX[32]/F");
  pEventTree->Branch("CMVDRecoPosY", CMVDRecoPosY, "CMVDRecoPosY[32]/F");
  pEventTree->Branch("CMVDRecoPosZ", CMVDRecoPosZ, "CMVDRecoPosZ[32]/F");


  pEventTree->Branch("CMVDWRecoPosX", CMVDWRecoPosX, "CMVDWRecoPosX[32]/F");
  pEventTree->Branch("CMVDWRecoPosY", CMVDWRecoPosY, "CMVDWRecoPosY[32]/F");
  pEventTree->Branch("CMVDWRecoPosZ", CMVDWRecoPosZ, "CMVDWRecoPosZ[32]/F");


  
  pEventTree->Branch("CMVDExpPosX", CMVDExpPosX, "CMVDExpPosX[32]/F");
  pEventTree->Branch("CMVDExpPosY", CMVDExpPosY, "CMVDExpPosY[32]/F");
  pEventTree->Branch("CMVDExpPosZ", CMVDExpPosZ, "CMVDExpPosZ[32]/F");

  pEventTree->Branch("CMVDExpnHit", CMVDExpnHit, "CMVDExpnHit[32]/I");

  pEventTree->Branch("CMVDisExpIn", CMVDisExpIn, "CMVDisExpIn[32]/I");
  pEventTree->Branch("CMVDisHit", CMVDisHit, "CMVDisHit[32]/I");

  pEventTree->Branch("CMVDpl2", CMVDpl2, "CMVDpl2[32]/I");
  pEventTree->Branch("CMVDClusterSize", CMVDClusterSize, "CMVDClusterSize[32]/I");

   pEventTree->Branch("CMVDClustsInLay", CMVDClustsInLay, "CMVDClustsInLay[32]/I");
   pEventTree->Branch("CMVDHitsInLay", CMVDHitsInLay, "CMVDHitsInLay[32]/I");
     



   

    pEventTree->Branch("cmv_nhit", &cmv_nhit, "cmv_nhit/i");
    pEventTree->Branch("cmv_hitid", cmv_hitid, "cmv_hitid[cmv_nhit]/I") ;
    pEventTree->Branch("cmv_hitpdgid", cmv_hitpdgid, "cmv_hitpdgid[cmv_nhit]/I");
    pEventTree->Branch("cmv_hitLeTim", cmv_hitLeTim, "cmv_hitLeTim[cmv_nhit]/F");
    pEventTree->Branch("cmv_hitRiTim", cmv_hitLeTim, "cmv_hitRiTim[cmv_nhit]/F");
    pEventTree->Branch("cmv_hitLePul", cmv_hitLePul, "cmv_hitLePul[cmv_nhit]/F");
    pEventTree->Branch("cmv_hitRiPul", cmv_hitLePul, "cmv_hitRiPul[cmv_nhit]/F");
    pEventTree->Branch("cmv_hitTrueposx", cmv_hitTrueposx, "cmv_hitTrueposx[cmv_nhit]/F");
    pEventTree->Branch("cmv_hitTrueposy", cmv_hitTrueposy, "cmv_hitTrueposy[cmv_nhit]/F");
    pEventTree->Branch("cmv_hitTrueposz", cmv_hitTrueposz, "cmv_hitTrueposz[cmv_nhit]/F");
    // pEventTree->Branch("cmv_hitmom", cmv_hitmom, "cmv_hitmom[cmv_nhit]/F");
    // pEventTree->Branch("cmv_hitthe", cmv_hitthe, "cmv_hitthe[cmv_nhit]/F");
    // pEventTree->Branch("cmv_hitphi", cmv_hitphi, "cmv_hitphi[cmv_nhit]/F");
    pEventTree->Branch("cmv_hitRecoposx", cmv_hitRecoposx, "cmv_hitRecoposx[cmv_nhit]/F");
    pEventTree->Branch("cmv_hitRecoposy", cmv_hitRecoposy, "cmv_hitRecoposy[cmv_nhit]/F");
    pEventTree->Branch("cmv_hitRecoposz", cmv_hitRecoposz, "cmv_hitRecoposz[cmv_nhit]/F");




    //cmv clust branches
    pEventTree->Branch("cmv_nclust", &cmv_nclust, "cmv_nclust/i");
    pEventTree->Branch("cmv_clustid", cmv_clustid, "cmv_clustid[cmv_nclust]/I") ;
    pEventTree->Branch("cmv_clustpdgid", cmv_clustpdgid, "cmv_clustpdgid[cmv_nclust]/I");
    pEventTree->Branch("cmv_clustLeTim", cmv_clustLeTim, "cmv_clustLeTim[cmv_nclust]/F");
    pEventTree->Branch("cmv_clustRiTim", cmv_clustLeTim, "cmv_clustRiTim[cmv_nclust]/F");
    pEventTree->Branch("cmv_clustLePul", cmv_clustLePul, "cmv_clustLePul[cmv_nclust]/F");
    pEventTree->Branch("cmv_clustRiPul", cmv_clustLePul, "cmv_clustRiPul[cmv_nclust]/F");
    pEventTree->Branch("cmv_clustTrueposx", cmv_clustTrueposx, "cmv_clustTrueposx[cmv_nclust]/F");
    pEventTree->Branch("cmv_clustTrueposy", cmv_clustTrueposy, "cmv_clustTrueposy[cmv_nclust]/F");
    pEventTree->Branch("cmv_clustTrueposz", cmv_clustTrueposz, "cmv_clustTrueposz[cmv_nclust]/F");
    // pEventTree->Branch("cmv_clustmom", cmv_clustmom, "cmv_clustmom[cmv_nclust]/F");
    // pEventTree->Branch("cmv_clustthe", cmv_clustthe, "cmv_clustthe[cmv_nclust]/F");
    // pEventTree->Branch("cmv_clustphi", cmv_clustphi, "cmv_clustphi[cmv_nclust]/F");
    pEventTree->Branch("cmv_clustRecoposx", cmv_clustRecoposx, "cmv_clustRecoposx[cmv_nclust]/F");
    pEventTree->Branch("cmv_clustRecoposy", cmv_clustRecoposy, "cmv_clustRecoposy[cmv_nclust]/F");
    pEventTree->Branch("cmv_clustRecoposz", cmv_clustRecoposz, "cmv_clustRecoposz[cmv_nclust]/F");
    pEventTree->Branch("cmv_clustsiz", cmv_clustsiz, "cmv_clustsiz[cmv_nclust]/I") ;








    // //extrapol brances
    // pEventTree->Branch("cmv_nexphit", &cmv_nexphit, "cmv_nexphit/i");
    // pEventTree->Branch("cmv_expid", cmv_expid, "cmv_expid[cmv_nexphit]/I") ;
    // pEventTree->Branch("cmv_Expposx", cmv_Expposx, "cmv_Expposx[cmv_nexphit]/F");
    // pEventTree->Branch("cmv_Expposy", cmv_Expposy, "cmv_Expposy[cmv_nexphit]/F");
    // pEventTree->Branch("cmv_Expposz", cmv_Expposz, "cmv_Expposz[cmv_nexphit]/F");
    //  pEventTree->Branch("cmv_DCAposx", cmv_DCAposx, "cmv_DCAposx[cmv_nexphit]/F");
    // pEventTree->Branch("cmv_DCAposy", cmv_DCAposy, "cmv_DCAposy[cmv_nexphit]/F");
    // pEventTree->Branch("cmv_DCAposz", cmv_DCAposz, "cmv_DCAposz[cmv_nexphit]/F");


    // pEventTree->Branch("distofclosapp", distofclosapp, "distofclosapp[cmv_nexphit]/F");
    // pEventTree->Branch("planeedge", planeedge, "planeedge[cmv_nexphit]/I") ;


 }// if(CardFile->GetCMVD()==1){

   
    pEventTree->Branch("tx_end",tx_end,"tx_end[ntrkt]/F");
    pEventTree->Branch("ty_end",ty_end,"ty_end[ntrkt]/F");



    pEventTree->Branch("mcxgnvx" ,mcxgnvx ,"mcxgnvx[ntrkt]/F");
    pEventTree->Branch("mcygnvx" ,mcygnvx ,"mcygnvx[ntrkt]/F");
    pEventTree->Branch("momgnvx" ,momgnvx ,"momgnvx[ntrkt]/F");
    pEventTree->Branch("thegnvx" ,thegnvx ,"thegnvx[ntrkt]/F");
    pEventTree->Branch("phignvx" ,phignvx ,"phignvx[ntrkt]/F");
    pEventTree->Branch("momgnend",momgnend,"momgnend[ntrkt]/F");
    pEventTree->Branch("thegnend",thegnend,"thegnend[ntrkt]/F");
    pEventTree->Branch("phignend",phignend,"phignend[ntrkt]/F");

    pEventTree->Branch("ntrkcl",ntrkcl,"ntrkcl[ntrkt]/I");
    pEventTree->Branch("ntrkst",ntrkst,"ntrkst[ntrkt]/I");
    pEventTree->Branch("range", &range, "range/F");

    pEventTree->Branch("tx", &tx, "tx/F");
    pEventTree->Branch("ty", &ty, "ty/F");
    pEventTree->Branch("xxin", &xxin, "xxin/F");
    pEventTree->Branch("yyin", &yyin, "yyin/F");
    pEventTree->Branch("txin", &txin, "txin/F");
    pEventTree->Branch("tyin", &tyin, "tyin/F");

    pEventTree->Branch("therr", &therr, "therr/F");
    pEventTree->Branch("pherr", &pherr, "pherr/F");

    pEventTree->Branch("atimslope", &atimslope, "atimslope/F");
    pEventTree->Branch("atiminter", &atiminter, "atiminter/F");


    pEventTree->Branch("xxerr", &xxerr, "xxerr/F");
    pEventTree->Branch("yyerr", &yyerr, "yyerr/F");
    pEventTree->Branch("txerr", &txerr, "txerr/F");
    pEventTree->Branch("tyerr", &tyerr, "tyerr/F");
    pEventTree->Branch("qperr", &qperr, "qperr/F");


   pEventTree->Branch("xxtxerr", &xxtxerr, "xxtxerr/F");
   pEventTree->Branch("xxtyerr", &xxtyerr, "xxtyerr/F");
    pEventTree->Branch("yytyerr", &yytyerr, "yytyerr/F");
    pEventTree->Branch("yytxerr", &yytxerr, "yytxerr/F");
    pEventTree->Branch("txtyerr", &txtyerr, "txtyerr/F");

    pEventTree->Branch("xxenderr", &xxenderr, "xxenderr/F");
    pEventTree->Branch("yyenderr", &yyenderr, "yyenderr/F");
    pEventTree->Branch("txenderr", &txenderr, "txenderr/F");
    pEventTree->Branch("tyenderr", &tyenderr, "tyenderr/F");
    pEventTree->Branch("qpenderr", &qpenderr, "qpenderr/F");

   }




    if(isInOut==2) {
      pEventTree->Branch("ntdc1x",&ntdc1x,"ntdc1x/i");
      pEventTree->Branch("tdcID1x",tdcID1x,"tdcID1x[ntdc1x]/I");
      pEventTree->Branch("TDCval1x",TDCval1x,"TDCval1x[ntdc1x]/F");
      pEventTree->Branch("ntstrp1x",&ntstrp1x,"ntstrp1x/i");
      pEventTree->Branch("StrpID1x",StrpID1x,"StrpID1x[ntstrp1x]/I");
      pEventTree->Branch("ntdc2x",&ntdc2x,"ntdc2x/i");
      pEventTree->Branch("tdcID2x",tdcID2x,"tdcID2x[ntdc2x]/I");
      pEventTree->Branch("TDCval2x",TDCval2x,"TDCval2x[ntdc2x]/F");
      pEventTree->Branch("ntstrp2x",&ntstrp2x,"ntstrp2x/i");
      pEventTree->Branch("StrpID2x",StrpID2x,"StrpID2x[ntstrp2x]/I");

      pEventTree->Branch("ntdc1y",&ntdc1y,"ntdc1y/i");
      pEventTree->Branch("tdcID1y",tdcID1y,"tdcID1y[ntdc1y]/I");
      pEventTree->Branch("TDCval1y",TDCval1y,"TDCval1y[ntdc1y]/F");
      pEventTree->Branch("ntstrp1y",&ntstrp1y,"ntstrp1y/i");
      pEventTree->Branch("StrpID1y",StrpID1y,"StrpID1y[ntstrp1y]/I");
      pEventTree->Branch("ntdc2y",&ntdc2y,"ntdc2y/i");
      pEventTree->Branch("tdcID2y",tdcID2y,"tdcID2y[ntdc2y]/I");
      pEventTree->Branch("TDCval2y",TDCval2y,"TDCval2y[ntdc2y]/F");
      pEventTree->Branch("ntstrp2y",&ntstrp2y,"ntstrp2y/i");
      pEventTree->Branch("StrpID2y",StrpID2y,"StrpID2y[ntstrp2y]/I");
      pEventTree->Branch("ntrecord1x",&ntrecord1x,"ntrecord1x/i");
      pEventTree->Branch("striprec1x",striprec1x,"striprec1x[ntrecord1x]/I");
      pEventTree->Branch("tdcrec1x",tdcrec1x,"tdcrec1x[ntrecord1x]/F");
      pEventTree->Branch("ntrecord1y",&ntrecord1y,"ntrecord1y/i");
      pEventTree->Branch("striprec1y",striprec1y,"striprec1y[ntrecord1y]/I");
      pEventTree->Branch("tdcrec1y",tdcrec1y,"tdcrec1y[ntrecord1y]/F");
      pEventTree->Branch("ntrecord2x",&ntrecord2x,"ntrecord2x/i");
      pEventTree->Branch("striprec2x",striprec2x,"striprec2x[ntrecord2x]/I");
      pEventTree->Branch("tdcrec2x",tdcrec2x,"tdcrec2x[ntrecord2x]/F");
      pEventTree->Branch("ntrecord2y",&ntrecord2y,"ntrecord2y/i");
      pEventTree->Branch("striprec2y",striprec2y,"striprec2y[ntrecord2y]/I");
      pEventTree->Branch("tdcrec2y",tdcrec2y,"tdcrec2y[ntrecord2y]/F");
    }

    if(isVisOut==1) {
      char outVisFile[300];
      sprintf(outVisFile,"%s.inh",outfile);
      pVisFile = new TFile(outVisFile, "RECREATE"); //VALGRIND
      if (!pVisFile) {
	cout << "Error opening .inh file !" << endl;
	exit(-1);
      } else {
	cout << "Vis Hits stored in: " << outVisFile << endl;
      }
      pVisFile->cd();
      EveCnt=0;//Event Counter
            if(!H) H = new Hits(); //VALGRIND
       if(!Hp) Hp= new HitPos();
      nloops=0;// number of tree fills
      if(!visTree){visTree = new TTree("Hitstree","Geant Hits File");}
      visTree->Branch("Hits_Branch","Hits",&H,1600000,2);
      H->Clear();
      Hp->Clear();
      H->ClearTracks();
      H->ENum=0;
    }


    char namex[200];
    for(int iaj=0; iaj<numberInLA; iaj++) {
      sprintf(namex,"hdifftime1_xy_%i",iaj);
      hdifftime1[iaj] = new TH1D(namex,namex,120,-25.,25.);
      sprintf(namex,"hdifftime2_xy_%i",iaj);
      hdifftime2[iaj] = new TH1D(namex,namex,120,-5.,5.);

      sprintf(namex,"hxtime_ext_%i",iaj);
      hxtime_ext[iaj] = new TH1D(namex,namex,120,-25.,25.);
      sprintf(namex,"hytime_ext_%i",iaj);
      hytime_ext[iaj] = new TH1D(namex,namex,120,-25.,25.);

      sprintf(namex,"hxpos_ext_%i",iaj);
      hxpos_ext[iaj] = new TH1D(namex,namex,120,-6.,6.);
      sprintf(namex,"hypos_ext_%i",iaj);
      hypos_ext[iaj] = new TH1D(namex,namex,120,-6.,6.);
      sprintf(namex,"hxpos_ext_kalman_%i",iaj);
      hxpos_ext_kalman[iaj] = new TH1D(namex,namex,120,-6.,6.);
      sprintf(namex,"hypos_ext_kalman_%i",iaj);
      hypos_ext_kalman[iaj] = new TH1D(namex,namex,120,-6.,6.);

      sprintf(namex,"h_hit_time_ext_%i",iaj);
      h_hit_time_ext[iaj] = new TH1D(namex,namex,120,-25.,25.);
      // for(int jak=0; jak<8; jak++) {
      // 	sprintf(namex,"xtdc_minus_ref_l%i_%i",iaj,jak);
      // 	xtdc_minus_ref[iaj][jak] = new TH1D(namex,namex,120,-20000.,200000);
      // 	sprintf(namex,"ytdc_minus_ref_l%i_%i",iaj,jak);
      // 	ytdc_minus_ref[iaj][jak] = new TH1D(namex,namex,120,-20000.,200000);
      // 	sprintf(namex,"tshift_xtdc_minus_ref%i_%i",iaj,jak);
      // 	tshift_xtdc_minus_ref[iaj][jak] = new TH1D(namex,namex,120,-20000.,200000);
      // 	sprintf(namex,"tshift_ytdc_minus_ref%i_%i",iaj,jak);
      // 	tshift_ytdc_minus_ref[iaj][jak] = new TH1D(namex,namex,120,-20000.,200000);
      // }
    }
    // DGap = new TH1D("DetGap"," Plot to see the performance of code at detgap ",10,0,10);  //asm
    // ShwXw	= new TH1D("ShwXw ","This is a distribution for X position resolution in shower region " ,80 , -3.14 , 3.14);  //asm
    // ShwYw	= new TH1D("ShwYw ","This is a distribution for Y position resolution in shower region " ,80 , -3.14 , 3.14);  //asm
    trk_edge = new TH2D("theta_dmn at trk_edge","theta_dmn at trk_edge" ,20, 0 ,20,200,0,3);  //asm
    trk_gap = new TH1D("difference in energy lost in the  dead space and the tracklength - )"," dpp-m*dxx + c" ,100, -5 ,5);  //asm
    // if(!pPosX)pPosX = new TH1D("deltaX", "#Delta X (cm)", 100, -52.5, 52.5);
    // if(!pPosY)pPosY = new TH1D("deltaY", "#Delta Y (cm)", 100, -52.5, 52.5);
    // if(!pPosZ)pPosZ = new TH1D("deltaZ", "#Delta Z (cm)", 100, -2.5, 2.5);
    // if(!pPosXX)pPosXX = new TH2D("deltaXX", "#Delta XX (cm)", 100, -2500, 2500, 100, -52.5, 52.5);
    // if(!pPosYY)pPosYY = new TH2D("deltaYY", "#Delta YY (cm)", 100, -500, 500, 100, -52.5, 52.5);
    // if(!pPosZZ)pPosZZ = new TH2D("deltaZZ", "#Delta ZZ (cm)", 100, -1000, 1000, 100, -2.5, 2.5);
  }


void MultiSimAnalysisDigi::OpenCollatedRootFile() {
  cout<<"Reading Collated Root Input file for digitization."<<endl;
  cout<<"collatedin "<<collatedIn<<endl;
  if(collatedIn) {
    collatedRootFile = new TFile("Collated_evtraw_20181120_20210813ai_5of8.root","read");
    if(!collatedRootFile) {
      cout << "Error opening collated file !" << endl;
      exit(-1);
    } else{
      cout<<"Collated file open successful..."<<endl;
      // char namex[200];
      // for(int iki=0; iki<numberInLA; iki++) {
      // 	sprintf(namex,"inefficiency_corx_l%i",iki);
      // 	inefficiency_corx[iki] = (TH2D*)collatedRootFile->Get(namex);
      // 	sprintf(namex,"inefficiency_uncx_l%i",iki);
      // 	inefficiency_uncx[iki] = (TH2D*)collatedRootFile->Get(namex);
      // 	// cout<<"inefficiency_uncx["<<iki<<"] "<<inefficiency_uncx[iki]<<endl;
      // 	sprintf(namex,"inefficiency_uncy_l%i",iki);
      // 	inefficiency_uncy[iki] = (TH2D*)collatedRootFile->Get(namex);
      // 	// cout<<"inefficiency_uncy["<<iki<<"] "<<inefficiency_uncy[iki]<<endl;
      // 	sprintf(namex,"triggereffi_xevt_l%i",iki);
      // 	triggereffi_xevt[iki] = (TH2D*)collatedRootFile->Get(namex);
      // 	sprintf(namex,"triggereffi_yevt_l%i",iki);
      // 	triggereffi_yevt[iki] = (TH2D*)collatedRootFile->Get(namex);
      // 	sprintf(namex,"strp_xmulsim_cor_l%i",iki);
      // 	strp_xmulsim_cor[iki] = (TH2D*)collatedRootFile->Get(namex);
      // 	sprintf(namex,"strp_ymulsim_cor_l%i",iki);
      // 	strp_ymulsim_cor[iki] = (TH2D*)collatedRootFile->Get(namex);
      // }

 char namex[200];
    for(int iki=0; iki<numberInLA; iki++) {
      sprintf(namex,"inefficiency_cor_m0_xr0_yr0_l%i",iki);
      inefficiency_corx[iki] = (TH2D*)collatedRootFile->Get(namex);
      sprintf(namex,"inefficiency_unc_m0_xr0_yr0_x%i",iki);
      inefficiency_uncx[iki] = (TH2D*)collatedRootFile->Get(namex);
      cout<<"inefficiency_uncx["<<iki<<"] "<<inefficiency_uncx[iki]<<endl;
      sprintf(namex,"inefficiency_unc_m0_xr0_yr0_y%i",iki);
      inefficiency_uncy[iki] = (TH2D*)collatedRootFile->Get(namex);
      cout<<"inefficiency_uncy["<<iki<<"] "<<inefficiency_uncy[iki]<<endl;
      sprintf(namex,"triggereffi_evt_m0_xr0_yr0_x%i",iki);
      triggereffi_xevt[iki] = (TH2D*)collatedRootFile->Get(namex);
      sprintf(namex,"triggereffi_evt_m0_xr0_yr0_y%i",iki);
      triggereffi_yevt[iki] = (TH2D*)collatedRootFile->Get(namex);
      sprintf(namex,"strp_xmulsim_cor_l%i",iki);
      strp_xmulsim_cor[iki] = (TH2D*)collatedRootFile->Get(namex);
      sprintf(namex,"strp_ymulsim_cor_l%i",iki);
      strp_ymulsim_cor[iki] = (TH2D*)collatedRootFile->Get(namex);
    }
    for(int iki=0; iki<numberInLA; iki++) {
      for(int ikj=0; ikj<16; ikj++) {
        for(int ikk=0; ikk<16; ikk++) {
          sprintf(namex,"blk_xmullaysim_m0_xr0_yr0_l%i_%i_%i",iki,ikj,ikk);
          block_xmulsim[iki][ikj][ikk] = (TH2D*)collatedRootFile->Get(namex);
          sprintf(namex,"blk_ymullaysim_m0_xr0_yr0_l%i_%i_%i",iki,ikj,ikk);
          block_ymulsim[iki][ikj][ikk] = (TH2D*)collatedRootFile->Get(namex);
        }
      }
    } // for(int iki=0; iki<numberInLA; iki++) {





    }
  }
}

void MultiSimAnalysisDigi::CloseInputRootFiles() {
  if (inputRootFile) {
    inputRootFile->cd();
    inputRootFile->Close(); //valgrind RSA
    // delete inputEventTree; inputEventTree=0;
    delete inputRootFile; inputRootFile=0;
    cout << "Input root file close !" << endl;
  } else {
    cout << "No input file to close.... !" << endl;
  }

}

void MultiSimAnalysisDigi::CloseOutputRootFiles() {
  if (pRootFile) {
    pRootFile->cd();
    if (pEventTree) { pEventTree->Write(); delete pEventTree; pEventTree=0;}
    pRootFile->Write();

    if (pPosX) {pPosX->Write(); delete pPosX; pPosX=0;}
    if (pPosY) {pPosY->Write(); delete pPosY; pPosY=0;}
    if (pPosZ) {pPosZ->Write(); delete pPosZ; pPosZ=0;}

    if (pPosXX) {pPosXX->Write(); delete pPosXX; pPosXX=0;}
    if (pPosYY) {pPosYY->Write(); delete pPosYY; pPosYY=0;}
    if (pPosZZ) {pPosZZ->Write(); delete pPosZZ; pPosZZ=0;}

    for (int i=0; i<20; i++) {if (pdedz[i]) {pdedz[i]->Write();  delete pdedz[i]; pdedz[i]=0;}};

    if (hitDist) {hitDist->Write(); delete hitDist; hitDist=0;}   //asm
    if (TrkDist) {TrkDist->Write(); delete TrkDist; TrkDist=0;}   //asm
    if (EffDist) {EffDist->Write(); delete EffDist; EffDist=0;}   //asm
    if (InoTrack_listsize) {InoTrack_listsize->Write(); delete InoTrack_listsize;  InoTrack_listsize=0;} //asm

    if (ShwXw) {ShwXw->Write(); delete ShwXw; ShwXw=0;}
    if (ShwYw) {ShwYw->Write(); delete ShwYw; ShwYw=0;}
    if (RC) {RC->Write(); delete RC; RC=0;}
    if (DGap) {DGap->Write(); delete DGap; DGap=0;}

    if (DeadStripX) {DeadStripX->Write(); delete DeadStripX; DeadStripX=0;}
    if (DeadStripY) {DeadStripY->Write(); delete DeadStripY; DeadStripY=0;}
    if (NoisyStripX) {NoisyStripX->Write(); delete NoisyStripX; NoisyStripX=0;}
    if (NoisyStripY) {NoisyStripY->Write(); delete NoisyStripY; NoisyStripY=0;}
    if (DiffTime) {DiffTime->Write(); delete DiffTime; DiffTime=0;}
    if (strpXtime) {strpXtime->Write(); delete strpXtime; strpXtime=0;}
    if (strpYtime) {strpYtime->Write(); delete strpYtime; strpYtime=0;}
    if (strpXtimeCorr) {strpXtimeCorr->Write(); delete strpXtimeCorr; strpXtimeCorr=0;}
    if (strpYtimeCorr) {strpYtimeCorr->Write(); delete strpYtimeCorr; strpYtimeCorr=0;}
    if (hitXtime) {hitXtime->Write(); delete hitXtime; hitXtime=0;}
    if (hitYtime) {hitYtime->Write(); delete hitYtime; hitYtime=0;}
    if (smagFieldX) {smagFieldX->Write(); delete smagFieldX; smagFieldX=0;}
    if (smagFieldY) {smagFieldY->Write(); delete smagFieldY; smagFieldY=0;}
    if (smag2dX) {smag2dX->Write(); delete smag2dX; smag2dX=0;}
    if (smag2dY) {smag2dY->Write(); delete smag2dY; smag2dY=0;}
    if (rmagFieldX) {rmagFieldX->Write(); delete rmagFieldX; rmagFieldX=0;}
    if (rmagFieldY) {rmagFieldY->Write(); delete rmagFieldY; rmagFieldY=0;}
    if (rmag2dX) {rmag2dX->Write(); delete rmag2dX; rmag2dX=0;}
    if (rmag2dY) {rmag2dY->Write(); delete rmag2dY; rmag2dY=0;}
    if (smag2dXYpixel_iron) {smag2dXYpixel_iron->Write(); delete smag2dXYpixel_iron; smag2dXYpixel_iron=0;}
    if (smag2dXYpixel_air) {smag2dXYpixel_air->Write(); delete smag2dXYpixel_air; smag2dXYpixel_air=0;}
    if (rmag2dXYpixel_iron) {rmag2dXYpixel_iron->Write(); delete rmag2dXYpixel_iron; rmag2dXYpixel_iron=0;}
    if (rmag2dXYpixel_air) {rmag2dXYpixel_air->Write(); delete rmag2dXYpixel_air; rmag2dXYpixel_air=0;}
    if (xyvsbxin) {xyvsbxin->Write(); delete xyvsbxin; xyvsbxin=0;}
    if (xyvsbyin) {xyvsbyin->Write(); delete xyvsbyin; xyvsbyin=0;}
    if (xyvsbxdiff) {xyvsbxdiff->Write(); delete xyvsbxdiff; xyvsbxdiff=0;}
    if (xyvsbydiff) {xyvsbydiff->Write(); delete xyvsbydiff; xyvsbydiff=0;}
    // if (xyvsbxindiff) {xyvsbxindiff->Write(); delete xyvsbxindiff; xyvsbxindiff=0;}
    // if (xyvsbyindiff) {xyvsbyindiff->Write(); delete xyvsbyindiff; xyvsbyindiff=0;}
    if (xyvsbxout) {xyvsbxout->Write(); delete xyvsbxout; xyvsbxout=0;}
    if (xyvsbyout) {xyvsbyout->Write(); delete xyvsbyout; xyvsbyout=0;}

    // for(int iai=0; iai<numberInLA; iai++) {
    //   if(hdifftime1[iai]) {hdifftime1[iai]->Write(); delete hdifftime1[iai]; hdifftime1[iai]=0;}
    //   if(hdifftime2[iai]) {hdifftime2[iai]->Write(); delete hdifftime2[iai]; hdifftime2[iai]=0;}
    //   if(hdifftime1[iai]) {hdifftime1[iai]->Write(); delete hdifftime1[iai]; hdifftime1[iai]=0;}
    pRootFile->Close();
    delete pRootFile; pRootFile=0;
    cout << "Output root file  close !" << endl;
  }  else {
    cout << "Output root file cannot close !" << endl;
  } // if (pRootFile) {

  if (pVisFile) {
    pVisFile->cd();
    visTree->Fill(); // fill tree
    H->Clear();
    Hp->Clear();
    H->ClearTracks();
    pVisFile->Write(); //VALGRIND
    cout<< "Write to .inh file  done"<<endl;
    if (visTree) {delete visTree; visTree=0;}
    pVisFile->Close();
    if (Hp) {delete Hp; Hp=0;}
    if (H) {delete H; H=0;}
    delete pVisFile; pVisFile=0;
    cout << "Visualization file  close !" << endl;
  } else {
    cout << "No output Hit Display Tree !" << endl;
  }

  if (collatedRootFile) {
    collatedRootFile->cd();
    for(int iki=0; iki<numberInLA; iki++) {
      delete inefficiency_corx[iki];
      delete inefficiency_uncx[iki];
      delete inefficiency_uncy[iki];
      delete triggereffi_xevt[iki];
      delete triggereffi_yevt[iki];
      delete strp_xmulsim_cor[iki];
      delete strp_ymulsim_cor[iki];
    }
    collatedRootFile->Close();
    delete collatedRootFile; collatedRootFile=0;
    cout<<"Collated root file closed."<<endl;
  } else {
    cout << "No collated histograms !" << endl;
  }

}

void MultiSimAnalysisDigi::SetCorrTimeError(G4double val) {
  cout<<"void MultiSimAnalysisDigi::SetCorrTimeError(G4double "<<val<<")"<<endl;
  CorrTimeError = val;
}

void MultiSimAnalysisDigi::SetUnCorrTimeError(G4double val) {
  cout<<"void MultiSimAnalysisDigi::SetUnCorrTimeError(G4double "<<val<<")"<<endl;
  UnCorrTimeError = val;
}

void MultiSimAnalysisDigi::SetTimeToDigiConvVal(G4double val) {
  cout<<"void MultiSimAnalysisDigi::SetTimeToDigiConvVal(G4double "<<val<<")"<<endl;
  TimeToDigiConv = val;
}

void MultiSimAnalysisDigi::SetSignalSpeedVal(G4double val) {
  cout<<"void MultiSimAnalysisDigi::SetSignalSpeedVal(G4double "<<val<<")"<<endl;
  SignalSpeed = val;
}


void MultiSimAnalysisDigi::SetPhotonSpeedVal(G4double val) {
  cout<<"void MultiSimAnalysis::SetPhotonSpeedVal(G4double "<<val<<")"<<endl;
  PhotonSpeed = val;
}

void MultiSimAnalysisDigi::SetCMVadctons(G4double val) {
  cout<<"void MultiSimAnalysis::SetCMVadctons(G4double "<<val<<")"<<endl;
  CMVadctons = val;
}









void MultiSimAnalysisDigi::SaveGenVisFile() {

  cout<<"void MultiSimAnalysisDigi::SaveGenVisFile() {"<<endl;
  if(isVisOut==1) {
    for(unsigned ij=0;ij<ngent;ij++) {
      H->NParticles++;
      Hp=  H->AddHits(0,0);
      Hp->TrackType=-14;
      Hp->ParCode= pidin[ij];
      Hp->ZZ= (((numberInLA*(parirlay[2]+parlay[2])*cm/m-parlay[2])- poszin[ij]*cm/m))/((parirlay[2]+parlay[2])*2*(1/m));
      Hp->XX= posxin[ij]*cm/m;
      Hp->YY= posyin[ij]*cm/m;
      Hp->pmag= momin[ij];
      Hp->pt= thein[ij];
      Hp->pp= phiin[ij];
    } // for(unsigned ij=0;ij< ngent;ij++) {
  } // if(isVisOut==1) {
}
