/*
  #include <cassert>
  #include <iostream>
  #include <iomanip>
  #include "InoStripCluster.h"
  using namespace std;
  //______________________________________________________________________

  //______________________________________________________________________
  InoStripCluster::InoStripCluster(){
  fUsed = -1;
  fView = -1; 
  fStripCluster =-1000;
  fPlane =-1000; 
  fXYPos =-1000; 
  fZPos = -1000; 
  iTrueTime = -1000; 
  iSmrTime = -1000; 
  fPulse = 1.; //Put it ONE, such that no effect on digital signal

  fMomentum=-1000.;
  fTheta = -1000.;
  fPhi =-1000.;

  fXgen = -1000.;
  fYgen = -1000.;
  fZgen = -1000.;

  fId    = -1;  
  fRPCmod =-1;
  }

  //______________________________________________________________________
  InoStripCluster::InoStripCluster(InoStripCluster* cd) {
  fUsed = cd->fUsed;
  fView = cd->fView;  
  fStripCluster = cd->fStripCluster;
  fPlane = cd->fPlane; 
  fXYPos = cd->fXYPos; 
  fZPos = cd->fZPos;  
  iTrueTime =  cd->iTrueTime; 
  iSmrTime =  cd->iSmrTime; 
  fId   =  cd->fId;
  fRPCmod =  cd->fRPCmod;
  iSmrTime =cd->iSmrTime;
  iTrueTime = cd->iTrueTime;
  

  }

  //______________________________________________________________________
  InoStripCluster::~InoStripCluster() {
  }

  //______________________________________________________________________
  InoStripCluster *InoStripCluster::DupHandle() const {
  return (new InoStripCluster(*this));
  }


  InoStripCluster::InoStripCluster(InoStrip* S1, InoStrip*S2, InoStrip*S3):
  fUsed(-1), fView(-1), fStripCluster(-1), fPlane(-1000), fXYPos(-1000), fZPos(-1000), iTrueTime(-1000), iSmrTime(-1000), fPulse(1), fMomentum(-1000), fTheta(-1000), fPhi(-1000), fXgen(-1000), fYgen(-1000), fZgen(-1000), fId(-1), fRPCmod(-1)
  {

  int nvalidHit=0;
  double pos[2]={0,0};//X/Y,Z
  fRPCmod=-1;	

  if(S1){
    
  if(fRPCmod<0){//to avoid repetition

  fRPCmod = S1->GetRPCmod();
  fView = S1->GetPlaneView();
  fPlane = S1->GetPlane();

  fMomentum =S1->GetMomentum();
  fTheta =S1->GetTheta();
  fPhi =S1->GetPhi();

  fXgen = S1->GetGenPosX();
  fYgen =S1->GetGenPosY();
  fZgen =S1->GetGenPosZ();
      
  }
	  
  fStripCluster = -1000;   
  pos[0] += S0->GetXYPos();
  pos[1] += S0->GetZPos();
  nvalidHit++;
	  
  }


  if(S2){


  if(fRPCmod<0){//to avoid repetition

  fRPCmod = S2->GetRPCmod();
  fView = S2->GetPlaneView();
  fPlane = S2->GetPlane();

  fMomentum =S2->GetMomentum();
  fTheta =S2->GetTheta();
  fPhi =S2->GetPhi();

  fXgen = S2->GetGenPosX();
  fYgen =S2->GetGenPosY();
  fZgen =S2->GetGenPosZ();
      
  }


  pos[0] += S2->GetXYPos();
  pos[1] += S2->GetZPos();
  nvalidHit++;

  }


  if(S3){

  if(fRPCmod<0){//to avoid repetition

  fRPCmod = S3->GetRPCmod();
  fView = S3->GetPlaneView();
  fPlane = S3->GetPlane();

  fMomentum =S3->GetMomentum();
  fTheta =S3->GetTheta();
  fPhi =S3->GetPhi();

  fXgen = S3->GetGenPosX();
  fYgen =S3->GetGenPosY();
  fZgen =S3->GetGenPosZ();
      
  }
    
  pos[0] += S3->GetXYPos();
  pos[1] += S3->GetZPos();
  nvalidHit++;

  }


  if (nvalidHit<1) {nvalidHit=1;}
  fXYPos =pos[0]/nvalidHit;
  fZPos =pos[1]/nvalidHit;

  fclustersize = nvalidHit;


  
  
  if(S1 && S2 && !(S3) ){ //any 2
  int strpid1 = S1->GetId();
  strpid1>>=8;
  int strpno1 = strpid1 & 0x7F;

  int strpid2 = S2->GetId();
  strpid2>>=8;
  int strpno2 = strpid2 & 0x7F;

  if(abs(strpno1-strpno2)==1){ // again check 
  cout<<"Found Two neighbs strip hits"<<endl;
   
  fRPCmod = S1->GetRPCmod();
  fView = S1->GetPlaneView();
  fStripCluster = -1000;
  fXYPos = 0.5*(S1->GetXYPos() + S2->GetXYPos());
  fPlane = S1->GetPlane();

      
  fZPos = 0.5*(S1->GetZPos() + S2->GetZPos());

  pdgid =-1000;

  iTrueTime =-1000;
  iSmrTime =-1000;
  fPulse =1;
      
  fMomentum =-1000;
  fTheta =-1000;
  fPhi =-1000;

  fXgen =-1000;
  fYgen =-1000;
  fZgen =-1000;
      
  fId =-1;
  fUsed = true;
  clustersize = 2;
    
  }
  else if( S2 && S3 && !(S1)){ //any 2
  int strpid1 = S2->GetId();
  strpid1>>=8;
  int strpno1 = strpid1 & 0x7F;

  int strpid2 = S3->GetId();
  strpid2>>=8;
  int strpno2 = strpid2 & 0x7F;
  
  if(abs(strpno1-strpno2)==1){
  cout<<"Found Two neighbs strip hits"<<endl;
  fRPCmod = S2->GetRPCmod();
  fView = S2->GetPlaneView();
    
    
  fStripCluster = -1000;

  fXYPos = 0.5*(S2->GetXYPos() + S3->GetXYPos());
  fPlane = S2->GetPlane();
    
    
  fZPos = 0.5*(S2->GetZPos() + S3->GetZPos());
      
   
    
  pdgid =-1000;
    
  iTrueTime =-1000;
  iSmrTime =-1000;
  fPulse =1;
  fMomentum =-1000;
  fTheta =-1000;
  fPhi =-1000;
    
  fXgen =-1000;
  fYgen =-1000;
  fZgen =-1000;
      
  fId =-1;
  fUsed = -1;
    
  }

  else if(S1 && S3 && !(S2) ){ //any 2
  int strpid1 = S1->GetId();
  strpid1>>=8;
  int strpno1 = strpid1 & 0x7F;

  int strpid2 = S3->GetId();
  strpid2>>=8;
  int strpno2 = strpid2 & 0x7F;

  if(abs(strpno1-strpno2)==1){
  cout<<"Found Two neighbs strip hits"<<endl;
  fRPCmod = S1->GetRPCmod();
    
  fView = S1->GetPlaneView();
 
  fStripCluster = -1000;
  fXYPos =-1000;
    
  fPlane =S1->GetPlane();
    
      
  fZPos = -1000;
    
  pdgid =-1000;
    
  iTrueTime =-1000;
  iSmrTime =-1000;
  fPulse =1;
  fMomentum =-1000;
  fTheta =-1000;
  fPhi =-1000;
    
  fXgen =-1000;
  fYgen =-1000;
  fZgen =-1000;
    
  fId =-1;
  fUsed = -1;
    
  }


  } else if(S1){
   
  fRPCmod = S1->GetRPCmod();
  fView = S1->GetPlaneView();
  fStripCluster = S1->GetStrip();
  fXYPos = S1->GetXYPos();
  fPlane = S1->GetPlane();
  fZPos = S1->GetZPos();
  pdgid = S1->GetpdgId();
      
  iTrueTime = S1->GetTrueTime();
  iSmrTime = S1->GetSmrTime();
  fPulse = S1->GetPulse();
  fMomentum = S1->GetMomentum();
  fTheta = S1->GetTheta();
  fPhi = S1->GetPhi();

  fXgen = S1->GetGenPosX();
  fYgen = S1->GetGenPosY();
  fZgen = S1->GetGenPosZ();
      
  fId = S1-> GetId();
  fUsed = S1->GetUsed();
    
  }// else if(S1){

  else if (S2){
   
  fRPCmod = S2->GetRPCmod();
  fView = S2->GetPlaneView();
  fStripCluster = S2->GetStrip();
  fXYPos = S2->GetXYPos();
  fPlane = S2->GetPlane();
  fZPos = S2->GetZPos();
  pdgid = S2->GetpdgId();
   
  iTrueTime = S2->GetTrueTime();
  iSmrTime = S2->GetSmrTime();
  fPulse = S2->GetPulse();
  fMomentum = S2->GetMomentum();
  fTheta = S2->GetTheta();
  fPhi = S2->GetPhi();
   
  fXgen = S2->GetGenPosX();
  fYgen = S2->GetGenPosY();
  fZgen = S2->GetGenPosZ();
   
  fId = S2-> GetId();
  fUsed = S2->GetUsed();
    

  }

  }

 











  //______________________________________________________________________
  void InoStripCluster::Trace(const char *c) const {
  std::cout<<"**********Begin InoStripCluster::Trace(\"" << c << "\")" << std::endl
  << "Information from InoStripCluster's CandHandle: " << std::endl;
  //  CandHandle::Trace(c);
  std::cout
  << "**********End InoStripCluster::Trace(\"" << c << "\")" << std::endl;
  }


  void InoStripCluster::Print() {
  //	cout<<"----------------------------------------------------------------------"<<endl;
  //	cout<<"Cluster combination "<<endl;
  cout<< "InoStripClusters():" 
  << "Plane "<< std::setw(3)<<  GetPlane()
  << "PlaneView "   << std::setw(3)<<  GetPlaneView()
  << "Strip "<< std::setw(4)<<  GetStripCluster()
  << "RPCmod "<< std::setw(4)<<  GetRPCmod()
  << " XYPos ="   <<std::setw(8)<<  GetXYPos() 
  << " ZPos =" <<std::setw(8)<<  GetZPos()
  << " truetime =" <<std::setw(8)<<  GetTrueTime()
  << " smrtime =" <<std::setw(8)<<  GetSmrTime()
  << " isUsed =" <<std::setw(8)<<  GetUsed()
  //   << " Z_Pos =" <<std::setw(8)<<  GetTruePosZ()
  // // << " chg="   <<std::setw(8)<<  GetPulse()
  //   << " time="  <<std::setw(8)<<  GetLeTime()
  << endl;
  //  cout<<"......................................................................"<<endl;
  }

*/


#include "InoStripCluster.h"
#include "InoStrip.h"
#include "TMath.h"
#include <iostream>

using namespace std;
InoStripCluster::InoStripCluster(InoStrip* strip) :
  fClusterSize(0),fSmrTime(0.0),
fXYPos(-999.9),fTrueTime(0.0),
  fPlane(-1),fRPCmod(-1),fBegXYStrip(-1), fEndXYStrip(-1), 
  fBegTrueTime(0.), fEndTrueTime(0.),
  fBegSmrTime(0.), fEndSmrTime(0.),

    fBegSmrTimeStrip(-1), fEndSmrTimeStrip(-1),
 fBegTrueTimeStrip(-1), fEndTrueTimeStrip(-1),
  fBegXYPos(-999.), fEndXYPos(999.),
  fView(-1), fMomentum(0), fTheta(0), fPhi(0), fXgen(-999.9), fYgen(-999.9), fZgen(-999.9),
   fZPos(0.), fPulse(0.), 
  fTrkFlag(0), fShwFlag(0),
  fTrkPlnFlag(0), fShwPlnFlag(0),
  //  fPlaneView(-1),
  fDigits(0),fNDFlag(1),       //asmQ what is the difference between fDigit and fView
  fXYPosErr(999.),
  //  fXPosErr(0.02/sqrt(12.)),
  //  fYPosErr(0.02/sqrt(12.)),
  InTrack(false), isStraight(true),
  StripXWidth(0.0196),
  StripYWidth(0.0196)
{
  paradef = micalDetectorParameterDef::AnPointer; //AAR: 
  StripXWidth = paradef->GetXStrwd()/1000;
  StripYWidth = paradef->GetYStrwd()/1000;
  //LayerThickness = (1/m)*2*(paradef->GetParlay(2)+paradef->GetParirlay(2)); //(1/m)*2*paradef->GetParlay(2);

  this->AddStrip(strip);
}


InoStripCluster::~InoStripCluster()
{
  cout<<"InoStripCluster::~InoStripCluster()<<endl";
  StripsInStripCluster.clear();
}


void InoStripCluster::AddStrip(InoStrip* strip) {
  if(StripsInStripCluster.size()==0) {
    StripsInStripCluster.push_back(strip);
    fPlane=strip->GetPlane();
    fRPCmod = strip->GetRPCmod();
    fView = strip->GetPlaneView();//X-Z or Y-Z plane
    

    fMomentum =strip->GetMomentum();
    fTheta =strip->GetTheta();
    fPhi =strip->GetPhi();

    fXgen = strip->GetGenPosX();
    fYgen =strip->GetGenPosY();
    fZgen =strip->GetGenPosZ();

    ////  if (strip->GetXYPosErr()<100) {
      fBegXYStrip=strip->GetStrip();
      fEndXYStrip=strip->GetStrip();

      fBegXYPos=strip->GetXYPos();
      fEndXYPos=strip->GetXYPos();
      fXYPos=strip->GetXYPos();
      
      fSmrTime = strip->GetSmrTime();
      fTrueTime = strip->GetTrueTime();

      fBegTrueTime=strip->GetTrueTime();
      fEndTrueTime=strip->GetTrueTime();
      
      fBegSmrTime=strip->GetSmrTime();
      fEndSmrTime=strip->GetSmrTime();
      
    fBegSmrTimeStrip=strip->GetStrip();
    fEndSmrTimeStrip=strip->GetStrip();
    
    fZPos=strip->GetZPos();
    fPulse = strip->GetPulse();
   
  } else {
    if(this->ContainsStrip(strip)==true) {return;}
    StripsInStripCluster.push_back(strip);   
    
    ////  if (strip->GetXYPosErr()<100) {
      if(strip->GetStrip()<fBegXYStrip) fBegXYStrip=strip->GetStrip();
      if(strip->GetStrip()>fEndXYStrip) fEndXYStrip=strip->GetStrip();
      if(strip->GetXYPos()<fBegXYPos) fBegXYPos=strip->GetXYPos();
      if(strip->GetXYPos()>fEndXYPos) fEndXYPos=strip->GetXYPos();
      // if (fView ==1) fView = 2;
      ////  }
    
      if(strip->GetTrueTime()<fBegTrueTime) {fBegTrueTime=strip->GetTrueTime();fBegTrueTimeStrip=strip->GetStrip();}
  if(strip->GetTrueTime()>fEndTrueTime) {fEndTrueTime=strip->GetTrueTime();fBegTrueTimeStrip=strip->GetStrip();}

    if(strip->GetSmrTime()<fBegSmrTime) {fBegSmrTime=strip->GetSmrTime();fBegSmrTimeStrip=strip->GetStrip();}
    if(strip->GetSmrTime()>fEndSmrTime) {fEndSmrTime=strip->GetSmrTime();fEndSmrTimeStrip=strip->GetStrip();}



 
   
    cout<<"fClusterSize:fSmrTime:strip->GetSmrTime() "<<fClusterSize<<" "<<fSmrTime<<" "<<strip->GetSmrTime() <<endl;
    fSmrTime = ( fClusterSize*fSmrTime + strip->GetSmrTime() )/(fClusterSize + 1); //Storing average smr time: raj
    
    fTrueTime = ( fClusterSize*fTrueTime + strip->GetTrueTime() )/(fClusterSize + 1); //Avg.
    fXYPos = ( fClusterSize*fXYPos + strip->GetXYPos() )/(fClusterSize + 1); //Avg.

    cout<<"Avg Smr/True time Pos "<<fSmrTime<<" "<<fTrueTime<<" "<<fXYPos<<endl;

    
  }
  
  //fDigit gives total # of X+Y  strips  //nXstrip -> total # of x strip in this cluster
  ////  if (strip->GetXYPosErr()<100) { 

  // fDigits += 1; // strip->GetCandStripHandle()->GetNDaughters();
  //   unsigned int nxystrip = GetXYEntries();
  //   fXYPos = (fXYPos*(nxystrip-1)+strip->GetXYPos())/(1.0*nxystrip);
  //   if (nxystrip==1) {  
  //     fXYPosErr = strip->GetXYPosErr();
  //   } else {
  //     fXYPosErr = strip->GetXYPosErr()/2; //GMA 09/02/09 Put these separately from data   //asm why so?
  //   }

    
    fPulse += strip->GetPulse();
    fClusterSize++;
  return;
}

bool InoStripCluster::ContainsStrip(InoStrip* strip) {
  for(unsigned int ij=0; ij<StripsInStripCluster.size(); ++ij) {
    if(strip==StripsInStripCluster[ij]) {return true;}
  }
  return false;
}

unsigned int InoStripCluster::GetXYEntries() {
  unsigned int nxystrip = 0;
  for(unsigned int ij=0; ij<StripsInStripCluster.size(); ++ij) {
    //    if (StripsInStripCluster[ij]->GetXYPosErr() < 100)
    nxystrip++;
  }
  return nxystrip;
}


InoStrip* InoStripCluster::GetStrip(unsigned int ij) const {
  cout<<ij<<" "<<StripsInStripCluster.size()<<endl;
  if(ij<StripsInStripCluster.size()) {return StripsInStripCluster[ij];}
  else {return 0;}

}

void InoStripCluster::Print() {

  cout<<" InoStripCluster():"
      <<" Plane = "<<GetPlane()<<" RPCmod "<<GetRPCmod()
      <<" PlaneView = "<<GetPlaneView()
      <<" Beg Strip = "<< std::setw(8)<<GetBegXYStrip() 
      <<" End Strip = "<< std::setw(8)<<GetEndXYStrip() 
      <<" Strip = "<< std::setw(8)<<GetStripCluster()     
      <<" Beg TrueTime = "<< std::setw(8)<<GetBegTrueTime() 
      <<" Beg SmrTime = "<< std::setw(8)<<GetBegSmrTime() 
      <<" Beg SmrTimeStrip = "<< std::setw(8)<<GetBegSmrTimeStrip()
      <<" End TrueTime = "<< std::setw(8)<<GetEndTrueTime() 
      <<" End SmrTime = "<< std::setw(8)<<GetEndSmrTime() 
      <<" End SmrTimeStrip = "<< std::setw(8)<<GetEndSmrTimeStrip() 
      <<" TrueTime = "<< std::setw(8)<<GetTrueTime() 
      <<" SmrTime = "<< std::setw(8)<<GetSmrTime() 
      <<" Beg  Pos = "<< std::setw(8)<<GetBegXYPos() 
      <<" End Pos = "<< std::setw(8)<<GetEndXYPos() 
      <<" XYPos = "<< std::setw(8)<<GetXYPos() 
    
      <<" ClusterSize = "<< std::setw(8)<<GetStripEntries()
      <<endl;
  // for(unsigned int ij=0; ij<StripsInStripCluster.size(); ij++) {
  //   StripsInStripCluster[ij]->Print();
  // }
 
}
