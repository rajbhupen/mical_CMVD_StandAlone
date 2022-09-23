/* #ifndef INOCLUSTER_H */
/* #define INOCLUSTER_H */
/* // ClusterCam */

/* #include <vector> */
/* #include "micalDetectorParameterDef.hh" */
/* using std::vector; */

/* class InoHit; */
/* class InoCluster */
/* { */
/*  public: */
/*   InoCluster(InoHit* hit); */
/*   virtual ~InoCluster(); */

/*   void AddHit(InoHit* hit); */
/*   bool ContainsHit(InoHit* hit); */
/*   void Print(); */
/*   int IsHitAssoc(InoHit* hit) const; */
/*   int IsShwAssoc(InoCluster* clust) const; */
/*   int IsTrkAssoc(InoCluster* clustm, InoCluster* clustp) const; */
/*   int IsDiffuseShwAssoc(InoCluster* clr) const; */

/*   unsigned int GetHitEntries() const {return HitsInCluster.size();} */
/*   unsigned int GetXEntries(); */
/*   unsigned int GetYEntries(); */
/*   unsigned int GetXProjEntries(); */
/*   unsigned int GetYProjEntries(); */
/*   unsigned int GetNXStripsInClust(); */
/*   unsigned int GetNYStripsInClust(); */

/*   InoHit* GetHit(unsigned int i) const; */
/*   int GetDigits() const { return fDigits; } */

/*   int GetZPlane() const {return fZPlane;} */
/*   int GetRPCmod() const {return fRPCmod;} */
/*   int GetBegXStrip() const {return fBegXStrip;} */
/*   int GetEndXStrip() const {return fEndXStrip;} */
/*   int GetBegYStrip() const {return fBegYStrip;} */
/*   int GetEndYStrip() const {return fEndYStrip;} */

/*   int GetView() const {return fView;} */

/*   double GetZPos() {return fZPos;} */
/*   double GetXPos() {return fXPos;} */
/*   double GetYPos() {return fYPos;} */

/*   void SetClusterNum(int val) {fClusterNum = val;} */
/*   int GetClusterNum() {return fClusterNum;} */
/*   double GetPulse() const { */
/*     if (fXPulse >0 && fYPulse >0) { */
/*       return 0.5*(fXPulse+fYPulse); */
/*     } else if (fXPulse >0) { */
/*       return fXPulse; */
/*     } else { */
/*       return fYPulse; */
/*     } */
/*   }; */

/*   double GetXPulse() {return fXPulse;} */
/*   double GetYPulse() {return fYPulse;} */

/*   double GetTime() const {return 0.5*(fBegTime+fEndTime);} */
/*   double GetBegTime() const {return fBegTime;} */
/*   double GetEndTime() const {return fEndTime;} */
/*   // For position correction added Jim */
/*   double GetBegXTime() const {return fBegXTime;} */
/*   double GetEndXTime() const {return fEndXTime;} */
/*   double GetBegYTime() const {return fBegYTime;} */
/*   double GetEndYTime() const {return fEndYTime;} */

/*   double GetBegXPos() const {return fBegXPos;} */
/*   double GetEndXPos() const {return fEndXPos;} */
/*   double GetBegYPos() const {return fBegYPos;} */
/*   double GetEndYPos() const {return fEndYPos;} */

/*   double GetXPosErr() const    {return fXPosErr;} */
/*   double GetYPosErr() const    {return fYPosErr;} */

/*   void SetXPosErr(double err) {fXPosErr = err;} */
/*   void SetYPosErr(double err) {fYPosErr = err;} */

/*   int GetTrkFlag() const {return fTrkFlag;} */
/*   int GetShwFlag() const {return fShwFlag;} */
/*   int GetTrkPlnFlag() const {return fTrkPlnFlag;} */
/*   int GetShwPlnFlag() const {return fShwPlnFlag;} */

/*   void SetTrkFlag(int flag) {fTrkFlag=flag;} */
/*   void SetShwFlag(int flag) {fShwFlag=flag;} */
/*   void SetTrkPlnFlag(int flag) {fTrkPlnFlag=flag;} */
/*   void SetShwPlnFlag(int flag) {fShwPlnFlag=flag;} */

/*   void SetNDFlag(int flag) {fNDFlag=flag;} */
/*   int GetNDFlag() const {return fNDFlag;} */
/*   void SetInTrack (bool tag) { InTrack=tag;} */
/*   bool GetInTrack () { return InTrack;} */
/*   void SetInShower (bool tag) { InShower=tag;} */
/*   bool GetInShower () { return InShower;} */

/*   void SetStraight (bool tag) { isStraight=tag;} */
/*   bool GetStraight () { return isStraight;} */

/*   int fZPlane; */
/*   int fRPCmod; */
/*   int fBegXStrip; */
/*   int fEndXStrip; */
/*   int fBegYStrip; */
/*   int fEndYStrip; */

/*   int kXTogether; */
/*   int kYTogether; */

/*   double fBegTime; */
/*   double fEndTime; */
/*   double fBegXTime; */
/*   double fEndXTime; */
/*   double fBegYTime; */
/*   double fEndYTime; */

/*   int fBegXTimeStrip; */
/*   int fEndXTimeStrip; */
/*   int fBegYTimeStrip; */
/*   int fEndYTimeStrip; */

/*   double fBegXPos; */
/*   double fEndXPos; */
/*   double fBegYPos; */
/*   double fEndYPos; */
/*   double fXPos; */
/*   double fYPos; */
/*   double fZPos; */
/*   double fXPulse; */
/*   double fYPulse; */
/*   int fTrkFlag; */
/*   int fShwFlag; */
/*   int fTrkPlnFlag; */
/*   int fShwPlnFlag; */

/*   int fDigits; */
/*   int fNDFlag; */

/*   double fXPosErr; //Error in X-position (m) */
/*   double fYPosErr; //Error in Y-position (m */
/*   int    fView;    // */
/*                    // 0 : only X(U)-axis hit */
/*                    // 1 : only Y(V)-axis hit */
/*                    // 2 : Both X and Y-axis hit */

/*   //GMA put in public place to get hits size in InoTrackFinder.cc file for informations of visualisation2 */
/*   vector<InoHit*> HitsInCluster; */
/*   bool isIdentical(InoCluster* icls); */
/*  private: */
/*   micalDetectorParameterDef *paradef; */
/*   bool   InTrack; */
/*   bool   InShower; */
/*   bool   isStraight; //Is it in straight section of track */
/*   bool debug_clust; */
/*   double StripXWidth; */
/*   double StripYWidth; */
/*   int fClusterNum; */
/* }; */

/* #endif */

///////////////////////////////////////////////////////
//New Code based on StripCluster

#ifndef INOCLUSTER_H
#define INOCLUSTER_H
// ClusterCam
#include "InoStripCluster.h"
#include <vector>
#include "micalDetectorParameterDef.hh"
#include "MultiSimAnalysisDigi.hh"
using std::vector;

class InoHit;
class InoCluster
{
 public:
  /* InoCluster(InoHit* hit); */
  /* virtual ~InoCluster(); */


  InoCluster();//rsa
  InoCluster(InoStripCluster* fx, InoStripCluster* fy);
  InoCluster(InoStripCluster* fx);
  
  InoCluster(InoCluster* hit);
  ~InoCluster();


  void Print();
  //  void AddHit(InoHit* hit);
  bool ContainsHit(InoHit* hit);

  int IsHitAssoc(InoHit* hit) const;
  int IsShwAssoc(InoCluster* clust) const;
  int IsTrkAssoc(InoCluster* clustm, InoCluster* clustp) const;
  int IsDiffuseShwAssoc(InoCluster* clr) const;


  unsigned int GetXEntries();
  unsigned int GetYEntries();
  unsigned int GetXProjEntries();
  unsigned int GetYProjEntries();
  unsigned int GetNXStripsInClust();
  unsigned int GetNYStripsInClust();
  
  InoStripCluster* GetXStripCluster() const {return fXStripCluster;}
  InoStripCluster* GetYStripCluster() const {return fYStripCluster;}
   
  int GetXStripClusterSize() const {return fXStripClusterSize;}
  int GetYStripClusterSize() const {return fYStripClusterSize;}
  
  int GetDigits() const { return fDigits; }

 
  int GetRPCmod() const {return fRPCmod;}
  int GetBegXStrip() const {return fBegXStrip;}
  int GetEndXStrip() const {return fEndXStrip;}
  int GetBegYStrip() const {return fBegYStrip;}
  int GetEndYStrip() const {return fEndYStrip;}

  int GetXpdgId() const {return fXpdgId;};
  int GetYpdgId() const {return fYpdgId;};

  
  int GetView() const {return fView;}

  double GetZPos() {return fZPos;}
  double GetXPos() {return fXPos-fXpOffset;}
  double GetYPos() {return fYPos-fYpOffset;}

  int GetZPlane() const    {return fZPlane;};
  void SetZPlane(int zpl) {fZPlane = zpl;};
  
  void SetMomentum(double f) {fMomentum = f;};
  void SetTheta(double f) {fTheta = f;};
  void SetPhi(double f) {fPhi = f;};

  double GetMomentum() {return fMomentum ;};
  double GetTheta() {return fTheta ;};
  double GetPhi() {return fPhi;};


  
  //jim_pos_corr
  void SetClusterNum(int val) {fClusterNum = val;}
  int GetClusterNum() {return fClusterNum;}
  //
  double GetPulse() const { 
    if (fXPulse >0 && fYPulse >0) {
      return 0.5*(fXPulse+fYPulse);
    } else if (fXPulse >0) {
      return fXPulse;
    } else {
      return fYPulse;
    }
  };  

  void SetXtOffset(double q) {fXtOffset=q;};
  void SetYtOffset(double q) {fYtOffset=q;};

  void SetXpOffset(double q) {fXpOffset=q;};
  void SetYpOffset(double q) {fYpOffset=q;};

  
  double GetXPulse() {return fXPulse;}
  double GetYPulse() {return fYPulse;}

  double GetTime() const;

  double GetBegXTime() const {return fBegXTime;} 
  double GetEndXTime() const {return fEndXTime;}
 
  double GetBegXTimeStrip() const {return fBegXTimeStrip;} 
  double GetEndXTimeStrip() const {return fEndXTimeStrip;}

 
  double GetBegYTime() const {return fBegYTime;}
  double GetEndYTime() const {return fEndYTime;}

  double GetBegYTimeStrip() const {return fBegYTimeStrip;} 
  double GetEndYTimeStrip() const {return fEndYTimeStrip;}
  
  double GetBegTime() const {return 0.5*(fBegXTime+fBegYTime);} 
  double GetEndTime() const {return 0.5*(fEndXTime+fEndYTime);}

  
  //not using anywhere
  double GetXTime() const {return DigiToTimeConv*fXTime;}; // Returns Smeared Time of X strip (in ns).
  double GetYTime() const {return DigiToTimeConv*fYTime;};

  
  double GetBegXPos() const {return fBegXPos-fXpOffset;}
  double GetEndXPos() const {return fEndXPos-fXpOffset;}
  double GetBegYPos() const {return fBegYPos-fYpOffset;}
  double GetEndYPos() const {return fEndYPos-fYpOffset;}

  double GetXPosErr() const    {return fXPosErr;}
  double GetYPosErr() const    {return fYPosErr;}

  void SetXPosErr(double err) {fXPosErr = err;}
  void SetYPosErr(double err) {fYPosErr = err;}

  int GetTrkFlag() const {return fTrkFlag;}
  int GetShwFlag() const {return fShwFlag;}
  int GetTrkPlnFlag() const {return fTrkPlnFlag;}
  int GetShwPlnFlag() const {return fShwPlnFlag;}

  void SetTrkFlag(int flag) {fTrkFlag=flag;}
  void SetShwFlag(int flag) {fShwFlag=flag;}
  void SetTrkPlnFlag(int flag) {fTrkPlnFlag=flag;}
  void SetShwPlnFlag(int flag) {fShwPlnFlag=flag;}
  
  void SetNDFlag(int flag) {fNDFlag=flag;}
  int GetNDFlag() const {return fNDFlag;}
  void SetInTrack (bool tag) { InTrack=tag;}
  bool GetInTrack () { return InTrack;}
  void SetInShower (bool tag) { InShower=tag;}
  bool GetInShower () { return InShower;}

  void SetStraight (bool tag) { isStraight=tag;}
  bool GetStraight () { return isStraight;}


  //GMA put in public place to get hits size in InoTrackFinder.cc file for informations of visualisation2

  
  bool isIdentical(InoCluster* icls);
 private:
  micalDetectorParameterDef *paradef;

  MultiSimAnalysisDigi* pAnalysis;

  int fDigits;
  int fRPCmod;
  int fBegXStrip;
  int fEndXStrip;
  int fBegYStrip;
  int fEndYStrip;
  
  int fXpdgId; // Particle ID which produces hit in X-strip
  int fYpdgId; // Particle ID which produces hit in Y-strip

  int    fView;    // 
                   // 0 : only X(U)-axis hit
                   // 1 : only Y(V)-axis hit
                   // 2 : Both X and Y-axis hit
  double fZPos;
  double fXPos;
  double fYPos;
  int    fZPlane;  //Z-plane
  double fMomentum; //Momentum of track which is behind this hit
  double fTheta;  //Theta
  double fPhi;    // Phi
  double fXPulse;
  double fYPulse;  
  double fBegXTime;
  double fEndXTime;
  double fBegYTime;
  double fEndYTime;

  
  //jim_pos_corr
  int kXTogether;
  int kYTogether;
  
  double fXTime;
  double fXTrueTime;
  double fYTime;
  double fYTrueTime;


  int fBegXTimeStrip;
  int fEndXTimeStrip;
  int fBegYTimeStrip;
  int fEndYTimeStrip;

  
  double fBegTime;
  double fEndTime;
  double fBegXPos;
  double fEndXPos;
  double fBegYPos;
  double fEndYPos;

  double fXPosErr; //Error in X-position (m)
  double fYPosErr; //Error in Y-position (m
  int fTrkFlag;
  int fShwFlag;
  int fTrkPlnFlag;
  int fShwPlnFlag;

  int fNDFlag;
  
  int fXPlane;     // =0 indicating strip measured X-axis planeview
  int fYPlane;     // =1 indicating strip measured y-axis
  
  bool   InTrack;
  bool   InShower;
  bool   isStraight; //Is it in straight section of track

  int fClusterNum;

  double DigiToTimeConv;
  double SignalSpeed;

  /* double fXtOffset=0.0; */
  /* double fYtOffset=0.0; */
   
  int fXStripClusterSize;
  int fYStripClusterSize;
 
  InoStripCluster* fXStripCluster;   //StripCluster for X-axis 
  InoStripCluster* fYStripCluster;   //StripCluster for Y-axis

  
  
  double StripXWidth;
  double StripYWidth;

  double fXtOffset;
  double fYtOffset;
  double fXpOffset;
  double fYpOffset;
  
};

#endif
