#ifndef INOSTRIPCLUSTER_H
#define INOSTRIPCLUSTER_H


#include <vector>
#include "micalDetectorParameterDef.hh"
#include "InoStrip.h"

using std::vector;
class InoStripCluster
{
 public:
  InoStripCluster(InoStrip* strip);
  virtual ~InoStripCluster();
 
  void Print();
  void AddStrip(InoStrip* strip);
  bool ContainsStrip(InoStrip* strip);

  unsigned int GetStripEntries() const {return StripsInStripCluster.size();}
  unsigned int GetXYEntries();

 InoStrip* GetStrip(unsigned int i) const;

  ////  unsigned int GetXYProjEntries();

  //  //unsigned int GetNXYStripsInClust();

  /* void   SetId(int id) {fId= id;} */
  /* int    GetId() { return fId;}; */
  
  int GetDigits() const { return fDigits; }

  int GetPlane() const {return fPlane;}
  int GetRPCmod() const {return fRPCmod;}
  int GetBegXYStrip() const {return fBegXYStrip;}
  int GetEndXYStrip() const {return fEndXYStrip;}
  int GetStripCluster() const{return 0.5*(fBegXYStrip+fEndXYStrip);}; //change later
 


  double GetZPos() {return fZPos;}



  double GetBegXYPos() const {return fBegXYPos;}
  double GetEndXYPos() const {return fEndXYPos;}

  double GetXYPos() const{return fXYPos;}; 
  
  

  double GetPulse() {return fPulse;};  


    double GetBegTrueTime() const {return fBegTrueTime;}
  double GetEndTrueTime() const {return fEndTrueTime;}

  double GetTrueTime() const {return fTrueTime;}

  double GetBegSmrTime() const {return fBegSmrTime;}
  double GetEndSmrTime() const {return fEndSmrTime;}


  
  int GetBegSmrTimeStrip() const {return fBegSmrTimeStrip;}
  int GetEndSmrTimeStrip() const {return fEndSmrTimeStrip;}

    
  int GetBegTrueTimeStrip() const {return fBegTrueTimeStrip;}
  int GetEndTrueTimeStrip() const {return fEndTrueTimeStrip;}


  double GetSmrTime() const {return fSmrTime;}


  double GetXYPosErr() const    {return fXYPosErr;}
  
  void   SetStrip(int f) { fStrip=f;};
  void SetXYPosErr(double err) {fXYPosErr = err;}
 
  void SetMomentum(double f) {fMomentum = f;};
  void SetTheta(double f) {fTheta = f;};
  void SetPhi(double f) {fPhi = f;};

  double GetMomentum() {return fMomentum ;};
  double GetTheta() {return fTheta ;};
  double GetPhi() {return fPhi;};

  void SetGenPosX(double f) {fXgen = f;}
  void SetGenPosY(double f) {fYgen = f;}
  void SetGenPosZ(double f) {fZgen = f;}


  double GetGenPosX() {return fXgen;}
  double GetGenPosY() {return fYgen;}
  double GetGenPosZ() {return fZgen;}

 int    GetPlaneView() const { return fView;};
    void   SetPlaneView( int f) {fView=f;};

  void   SetpdgId(int id) {pdgid= id;}
  int    GetpdgId() { return pdgid;};
  
 unsigned int fClusterSize;    
    double fMomentum; //Momentum of track which is behind this hit
  double fTheta;  //Theta
  double fPhi;    // Phi
  int    pdgid;  //pdgId of the particle which produce this hit
  double fXgen; //Generated track x-coordinate
  double fYgen; //Generated track y-coordinate
  double fZgen; //Generated track z-coordinate
  int    fId;
  int    fView;   // 0/1 for X/Y-axis
  int fPlane;
  int fRPCmod;
  int fBegXYStrip;
  int fEndXYStrip;

  double fSmrTime;
  double fTrueTime;
  double fBegTrueTime;
  double fEndTrueTime;

  double fBegSmrTime;
  double fEndSmrTime;

  int fBegSmrTimeStrip;
  int fEndSmrTimeStrip;
 

  int fBegTrueTimeStrip;
  int fEndTrueTimeStrip;
 
 
  double fBegXYPos;
  double fEndXYPos;
   double fXYPos;
  double fZPos;
  double fPulse;
  int fTrkFlag;
  int fShwFlag;
  int fTrkPlnFlag;
  int fShwPlnFlag;

  int fDigits;
  int fNDFlag;

  double fXYPosErr; //Error in X-position (m)
  
  
  //GMA put in public place to get strips size in InoTrackFinder.cc file for informations of visualisation2
  vector<InoStrip*> StripsInStripCluster;
  bool isIdentical(InoStripCluster* icls);
 private:
  micalDetectorParameterDef *paradef;
  int fStrip;
  bool   InTrack;
  bool   InShower;
  bool   isStraight; //Is it in straight section of track
  double StripXWidth;
  double StripYWidth;

};

#endif // INOSTRIPCLUSTER_H
