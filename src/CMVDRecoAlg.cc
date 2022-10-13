#include "CMVDRecoAlg.hh"

CMVDRecoAlg::CMVDRecoAlg(int isInOut) {
  grecoi = GeneralRecoInfo::GnPointer;
  pAnalysis = MultiSimAnalysisDigi::AnPointer;
  MessengerFile = ParameterMessenger::AnPointer;
  paradef = micalDetectorParameterDef::AnPointer;
   SipmHit_pointer = SipmHit_Manager::APointer; 


  
 CmvStrip_pointer = CmvStrip_Manager::APointer;
  
  //  DetectorType = paradef->GetDetectorType();
  isData = (isInOut==2) ? 1 : 0;

 

 for(int ii=0; ii<3; ii++) {
  partopscint[ii]=paradef->GetPartopscint(ii);
}

	
  for(int op=0; op<3;op++) {
    partopscint[op] = paradef->partopscint[op];
  }
	
  AirGapScintTop= paradef->AirGapScintTop;
	
       
  int jmax;
  for(int i =0;i<7;i++){//<4
   
    jmax =  (i==0)? 4:3;
    // if(i==0){
    //   jmax=4;
    // }
    // else{
    //   jmax=3;
    // }
    for(int j=0;j<jmax;j++){
      for(int k=0;k<3;k++){
	PhyVolGlPos[i][j][k] = paradef->ScintLayGPos[i][j][k];
	// cout<< PhyVolGlPos[i][j][k]<<" ";
      }//k
      // cout<<endl;
    }//j
    // cout<<endl<<endl;
  }//i
	
  NoScntStrpTop = paradef->GetNoScntStrpTop();//88
  NoScntStrpSide = paradef->GetNoScntStrpSide();//40
  NoScntStrpSideSmallay= paradef->GetNoScntStrpSideSmallay();	 //8
  SidePlaneHalfLength = paradef->GetSidePlaneHalfLength();
  SideSmallPlaneHalfLength = paradef->GetSideSmallPlaneHalfLength();
  ScntLayShifSide = paradef->GetScntLayShifSide();

  
}
CMVDRecoAlg::~CMVDRecoAlg() {
  SipmHit_pointer->SipmHit_list.clear();
  if (SipmHit_pointer) {delete SipmHit_pointer; SipmHit_pointer=0;}

  for (unsigned ij=0; ij<CmvHit_pointer->CmvHit_list.size(); ij++) {
    if (CmvHit_pointer->CmvHit_list[ij]) {
      cout <<"ij "<< ij<<" "<<CmvHit_pointer->CmvHit_list.size()<<endl;
      delete CmvHit_pointer->CmvHit_list[ij]; CmvHit_pointer->CmvHit_list[ij]=0;
    }
  }

  CmvHit_pointer->CmvHit_list.clear();
  if (CmvHit_pointer) {delete CmvHit_pointer; CmvHit_pointer=0;}




  for (unsigned ij=0; ij<CmvCluster_pointer->CmvCluster_list.size(); ij++) {
    if (CmvCluster_pointer->CmvCluster_list[ij]) {
      cout <<"ij "<< ij<<" "<<CmvCluster_pointer->CmvCluster_list.size()<<endl;
      delete CmvCluster_pointer->CmvCluster_list[ij]; CmvCluster_pointer->CmvCluster_list[ij]=0;
    }
  }

  CmvCluster_pointer->CmvCluster_list.clear();
  if (CmvCluster_pointer) {delete CmvCluster_pointer; CmvCluster_pointer=0;}
  
  cout<<"check cluster pointer"<<endl;

  
  // for (unsigned ij=0; ij<CmvLayExtra_pointer->CmvLayExtra_list.size(); ij++) {
  //   if (CmvLayExtra_pointer->CmvLayExtra_list[ij]) {
  //       cout <<"ij "<< ij<<" "<<CmvLayExtra_pointer->CmvLayExtra_list.size()<<endl;
  //     delete CmvLayExtra_pointer->CmvLayExtra_list[ij]; CmvLayExtra_pointer->CmvLayExtra_list[ij]=0;
  //   }
  // }

  // CmvLayExtra_pointer->CmvLayExtra_list.clear();
  // if (CmvLayExtra_pointer) {delete CmvLayExtra_pointer; CmvLayExtra_pointer=0;}


   cout << "Deleting CMVDRecoAlg ..." << endl;


   
}



void CMVDRecoAlg::ReadCMVDdata(int evt){
  cout<<"void CMVDRecoAlg::ReadCMVDdata(int dataMC){"<<endl;
  SipmHit_pointer = new SipmHit_Manager();

    pAnalysis->inputRootFile->cd();
   pAnalysis->inputEventTree->GetEntry(evt);// Already did while reading mical data


  if(isData){}
  else{
  for(unsigned ij=0;ij<pAnalysis->cmv_ndigihit;ij++) {
      //      unsigned istrp = pAnalysis->stripid[ij];
      SipmHit*  sipmht = new SipmHit(); //VALGRIND
      sipmht->SetId(pAnalysis->cmv_sipmid[ij]);
      sipmht->SetpdgId(pAnalysis->cmv_digipdgid[ij]);
      sipmht->SetTimePulse(pAnalysis->cmv_digitimpul[ij]);
      sipmht->SetXPos(pAnalysis->cmv_digiposx[ij]);
      sipmht->SetYPos(pAnalysis->cmv_digiposy[ij]);
      sipmht->SetZPos(pAnalysis->cmv_digiposz[ij]);
      sipmht->SetSimMom(pAnalysis->cmv_digimom[ij]);
      sipmht->SetSimThe(pAnalysis->cmv_digithe[ij]);
      sipmht->SetSimPhi(pAnalysis->cmv_digiphi[ij]);
      sipmht->SetXLocPos(pAnalysis->cmv_digilocx[ij]);
      sipmht->SetYLocPos(pAnalysis->cmv_digilocy[ij]);
    sipmht->SetZLocPos(pAnalysis->cmv_digilocz[ij]);
      int isipmid = sipmht->GetId();
      int oldid=-1;
      for (unsigned int jk=0; jk<SipmHit_pointer->SipmHit_list.size(); jk++) {
	if (isipmid==SipmHit_pointer->SipmHit_list[jk]->GetId()) {
	  SipmHit_pointer->SipmHit_list[jk]->Update(sipmht->GetPulse(), sipmht->GetTime());
	  oldid = jk;
	  break;
	}
      }
      if (oldid<0 && sipmht->GetPulse() > 0.0001) {
	//	sipmht->Print();
	SipmHit_pointer->SipmHit_list.push_back(sipmht);
      }
      //
  }// for(unsigned ij=0;ij<pAnalysis->cmv_ndigihit;ij++) {


  
  }

    cout<<SipmHit_pointer<<endl;

}



bool CMVDRecoAlg::LinePlaneInt(double* Line, double* Plane, double* Point){
  //  cout<<"bool micalEventAction::LinePlaneInt(double* Line, double* Plane, double* Point){"<<endl;
  //	G4double Dist;
  //	G4double a, b;
  bool ok;
	
  double b = Line[3]*Plane[3] + Line[4]*Plane[4] + Line[5]*Plane[5];
  //  cout<<Line[3]<<" "<<Plane[3]<<" "<<Line[4]<<" "<<Plane[4]<<" "<<Line[5]<<" "<<Plane[5]<<" "<<b<<endl;
  // ok= (fabs(b) > 1e-10) ? 1 : 0;//
  
  ok= ((1.0*b)> 1e-10) ? 1 : 0;//#we only do intersection with those planes which muon passed first..
  if(ok==1) {
    double a=(Plane[0]-Line[0])*Plane[3] +
      (Plane[1]-Line[1])*Plane[4] +
      (Plane[2]-Line[2])*Plane[5];
    //    cout<<"a "<<a<<endl;
    G4double Dist = a/b;
    //    cout<<"dist: "<<Dist<<endl;
    
    Point[0] = Line[0] + Line[3]*Dist;
    Point[1] = Line[1] + Line[4]*Dist;
    Point[2] = Line[2] + Line[5]*Dist;
    //    cout<<"Point "<<Point[0]<<" "<<Point[1]<<" "<<Point[2]<<endl;
  } else {
    //    cout<<"Setting Point =-100000 "<< endl;
    Point[0]=-1000000; Point[1]=-1000000; Point[2]=-1000000;
  }
  return ok;
}
//


void CMVDRecoAlg::CMVD_Extrapolation(){   
 
  
  cout<<"......................CMVD straight line extrapolation........................."<<endl;

  cout<<"Store in CmvClusterBank"<<endl;
  cout<<"CmvCluster_pointer->CmvCluster_list.size() "<<CmvCluster_pointer->CmvCluster_list.size()<<endl;
  for(int pl=0;pl<7;pl++){
    for(int la=0;la<4;la++){

      CmvClusterBank[pl][la].clear();


    }



  }
  for (unsigned int ix=0; ix<CmvCluster_pointer->CmvCluster_list.size(); ix++) {
    CmvCluster* tmpcluster = CmvCluster_pointer->CmvCluster_list[ix];
    tmpcluster->Print();
    CmvClusterBank[tmpcluster->GetPlane()-1][tmpcluster->GetLayer()].push_back(tmpcluster);
  }
  cout<<" Fill # of hits and clusters in each layer "<<endl;
  unsigned int arcnt=0;
  for(int ij=0;ij<7;ij++){
    if(ij==4)continue; //Reserved for Front Wall
    for(int jk=0;jk<4;jk++){
      if(ij>0 && jk>2) continue; //only 3 layers in side walls
      arcnt =ij+1;
      arcnt<<=2;
      arcnt+=jk;
      cout<<"arcnt "<<arcnt<<endl;
      pAnalysis->CMVDHitsInLay[arcnt] = CmvHitBank[ij][jk].size();
      pAnalysis->CMVDClustsInLay[arcnt] = CmvClusterBank[ij][jk].size();
      
      //      cout<<ij<<" "<<jk<<" "<<pAnalysis->CMVDHitsInLay[arcnt]<<" "<<pAnalysis->CMVDClustsInLay[arcnt]<<" "<<CmvHitBank[ij][jk].size()<<" "<<CmvClusterBank[ij][jk].size()<<endl;
  }
  }
  
  
  
  int ijmax=0;
  cout<<"check ab "<<pAnalysis->ntrkt<<endl;
  for (unsigned jk=0; jk<pAnalysis->ntrkt ; jk++) {

    pAnalysis->CMVDnLayer =pAnalysis->cmvdnlayermx;
    cout<<"ntrkt CMVDnLayer "<<pAnalysis->ntrkt<<" "<<pAnalysis->CMVDnLayer<< endl;
  
    //		double    momvx = pAnalysis->momvx[jk];
    //This theta phi represents a track going downward. The dirvector using this theta phi matched with generated direction vector (in PGA after...)
    double theta =pAnalysis->thevx[jk];
    double phi = pAnalysis->phivx[jk];
    double posx = pAnalysis->posxvx[jk]*m; //these were stored in metre and convert in mm unit
    double posy =pAnalysis->posyvx[jk]*m;
    double posz= pAnalysis->poszvx[jk]*mm;//already stored in mm only


    double atimslope = pAnalysis->atimslope[jk];
    double atiminter = pAnalysis->atiminter[jk];


    cout<<"atimslope "<<atimslope<<"atiminter "<<atiminter<<endl;
    cout<<"...jk..... "<<jk<<" "<<posx<<" "<<posy<<" "<<posz<<" "<<theta<<" "<<phi<<endl;
    //err calculation:
    
    double therr = pAnalysis->therr[jk];
    double pherr = pAnalysis->pherr[jk];
    // these are errors i.e. just sigma and are stored in meters..
    double xxerr = pAnalysis->xxerr[jk]*m;
    double yyerr = pAnalysis->yyerr[jk]*m;
    double txerr = pAnalysis->txerr[jk];// x = a + b z , err in b (no units)
    double tyerr = pAnalysis->tyerr[jk];// y = a + b z , err in b (no units)
    cout<<"error in positions and angles  "<<xxerr<<" "<<yyerr<<" "<<txerr<<" "<<tyerr<<" "<<therr<<" "<<pherr<<endl;
    
    double xxtxerr = pAnalysis->xxtxerr[jk]*m;
    double xxtyerr = pAnalysis->xxtyerr[jk]*m;
    double yytyerr = pAnalysis->yytyerr[jk]*m;
    double yytxerr = pAnalysis->yytxerr[jk]*m;
    double txtyerr = pAnalysis->txtyerr[jk];
    cout<<"xxtxerr "<<xxtxerr<<"yytyerr "<<yytyerr<<endl;
    //err cal.

    //    pAnalysis->chisq[jk]=-1.0; pAnalysis->chisq2[jk]=-1.0;
    //    pAnalysis->posxvx[jk]=0.0;  pAnalysis->posyvx[jk]=0.0;  pAnalysis->poszvx[jk]=0.0;
   
   for(int pq=0;pq<32;pq++){
      pAnalysis->CMVDisExpIn[pq] = 0;
      pAnalysis->CMVDisHit[pq] = 0;
    }



    //convert theta and phi to dxdz and dydz
    
    double dxdz = tan(theta)*cos(phi);   
    double dydz = tan(theta)*sin(phi);
    
    cout<< "Downward theta "<< theta<<" phi  "<< phi<< " dxdz "<< dxdz<<" dydz "<<dydz<<endl;


    double PI = acos(-1.0);
	
    //This gives a downward going track
    // G4ThreeVector dirVector1(0,0,1);       
    // dirVector1.setTheta(theta);
    // dirVector1.setPhi(phi);
    //We need to do extrapolation for those cmvd plane which a muon have passed before entering RPC stack. Thus we need to reverse sign of dirvector(parity) thus we add Pi to phi and do theta-pi..


    double dxdz_upward = tan(PI-theta)*cos(phi+PI);   //modified on 10.01.2022
    double dydz_upward = tan(PI-theta)*sin(phi+PI);
    //slope doesnt change on inverting all xyz to negative
    cout<< "Upward  theta "<< PI-theta<<" phi  "<< phi+PI<< " dxdz "<< dxdz_upward<<" dydz "<<dydz_upward<<endl;



    
    G4ThreeVector dirVector(0,0,1);  //     
    dirVector.setTheta(PI-theta);// this works atleast for SL
    dirVector.setPhi(phi+PI);

    //   dirVector.setTheta(theta);
    //  dirVector.setPhi(phi);
    

    cout<<"dirvector "<<dirVector<<" "<<endl;
    //thus this dirvector represents a track going backward towards the vertex/gen point.
    G4double G_Point[3]={0};
    G4double Ptxerr, Ptyerr, Ptzerr =0;
    G4double Line[6]={posx,posy,posz,dirVector.x(),dirVector.y(),dirVector.z()};
    
    cout<<"Line: "<<Line[0]<<" "<<Line[1]<<" "<<Line[2]<<" "<<Line[3]<<" "<<Line[4]<<" "<<Line[5]<<endl;



    unsigned int layid = 0;
    double  xhat=0,yhat=0,zhat=0; //area unit vector of planes

    double layhalflength;
    for(int ijk=0;ijk<7;ijk++){//loc_no loop 0:Top, 1:Left 2:Right 3:Back 4:Front 5: miniLeft 6:miniRight //4
      if(ijk==4){continue;} //Reserved for Front Wall
      cout<<"locno: "<<ijk<<endl;
      if(ijk==0){ijmax=4;xhat=0;yhat=0;zhat=1;layhalflength = paradef->GetTopPlaneHalfLength();}//top
      else if(ijk==1 ){ijmax=3;xhat=-1;yhat=0;zhat=0; layhalflength = paradef->GetSidePlaneHalfLength(); }//left
      else if(ijk==2){ijmax=3;xhat=1;yhat=0;zhat=0;layhalflength = paradef->GetSidePlaneHalfLength(); }//right
      else if(ijk==3){ijmax=3;xhat=0;yhat=1;zhat=0;layhalflength = paradef->GetSidePlaneHalfLength();}//back
      else if(ijk==4){ijmax=3;xhat=0;yhat=-1;zhat=0;layhalflength = paradef->GetSidePlaneHalfLength();}//Front
      else if(ijk==5){ijmax=3;xhat=-1;yhat=0;zhat=0;layhalflength = paradef->GetSideSmallPlaneHalfLength();}//miniLeft
      else if(ijk==6){ijmax=3;xhat=1;yhat=0;zhat=0;layhalflength = paradef->GetSideSmallPlaneHalfLength();}//miniRight

      for(int ij=0;ij<ijmax;ij++){//layer loop
	double diffmx=1000000;
	cout<<" Layer No. "<<ij<<endl;
	//	cout<<ijk<<" "<<ij<<endl;
	G4double Plane[6]={PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2],xhat,yhat,zhat};

        cout<<"Plane: "<<Plane[0]<<" "<<Plane[1]<<" "<<Plane[2]<<" "<<Plane[3]<<" "<<Plane[4]<<" "<<Plane[5]<<endl;
	

        bool pl2 = LinePlaneInt (Line, Plane, G_Point);

	double layhalfbreadth = paradef->partopscint[1];     //2300 for left right
	if(ijk==3){
	  layhalfbreadth = paradef->partopscint[1]+50;//4.7 m    

	}
	if(ijk==0){
	  layhalfbreadth = paradef->partopscint[1]-50;//4.5 m    

	}

	if(ijk==5 || ijk==6){
	  layhalfbreadth = 1000;

	}



		
		cout<<"Loc_no: "<<ijk<<" layhalflength "<<layhalflength<<"layhalfbreadth  "<<layhalfbreadth<<endl;

	// cout<<"Point: "<<G_Point[0]<<" "<<G_Point[1]<<" "<<G_Point[2]<<endl;
	vector <double> edge[4];
    
	if(ijk==0){

	  edge[0] =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]+layhalfbreadth, PhyVolGlPos[ijk][ij][2],-1,0,0}; //backside edge
	  edge[1] =	      {PhyVolGlPos[ijk][ij][0]-layhalflength, PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2],0,-1,0}; //leftside edge
	  edge[2] =	      {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]-layhalfbreadth, PhyVolGlPos[ijk][ij][2],1,0,0}; //frontside edge
	  edge[3] =	    {PhyVolGlPos[ijk][ij][0]+layhalflength, PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2],0,1,0}; //rightside edge


	}
	      
	if(ijk==1){

	  edge[0]    =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]+layhalflength,0,1,0}; //topside edge
	  edge[1] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]+layhalfbreadth, PhyVolGlPos[ijk][ij][2],0,0,-1}; //backside edge
	  edge[2] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]-layhalflength,0,-1,0}; //bottomside edge
	  edge[3] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]-layhalfbreadth, PhyVolGlPos[ijk][ij][2],0,0,1} ; //frontside edge
		
		
	}
	      
	if(ijk==2){
		
	  edge[0]  =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]+layhalflength,0,-1,0}; //topside edge
	  edge[1] =     {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]+layhalfbreadth, PhyVolGlPos[ijk][ij][2],0,0,-1}; //backside edge
	  edge[2] =    {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]-layhalflength,0,1,0}; //bottomside edge
	  edge[3] =   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]-layhalfbreadth, PhyVolGlPos[ijk][ij][2],0,0,1}; //frontside edge
		
		  
	}
	      
	if(ijk==3){
		
	  edge[0]   =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]+layhalflength,1,0,0}; //topside edge
	  edge[1] =	{PhyVolGlPos[ijk][ij][0]+layhalfbreadth, PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2],0,0,-1}; //rightside edge
	  edge[2] =	  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]-layhalflength,-1,0,0}; //bottomside edge
	  edge[3] =    {PhyVolGlPos[ijk][ij][0]-layhalfbreadth, PhyVolGlPos[ijk][ij][1],PhyVolGlPos[ijk][ij][2],0,0,1}; //leftside edge
		
	       
	}
	//miniLeft
	if(ijk==5){

	  edge[0]    =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]+layhalfbreadth,0,1,0}; //topside edge
	  edge[1] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]+layhalflength, PhyVolGlPos[ijk][ij][2],0,0,-1}; //backside edge
	  edge[2] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]-layhalfbreadth,0,-1,0}; //bottomside edge
	  edge[3] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]-layhalflength, PhyVolGlPos[ijk][ij][2],0,0,1} ; //frontside edge
		
		
	}
	      
	if(ijk==6){

	  edge[0]  =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]+layhalfbreadth,0,-1,0}; //topside edge
	  edge[1] =     {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]+layhalflength, PhyVolGlPos[ijk][ij][2],0,0,-1}; //backside edge
	  edge[2] =    {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]-layhalfbreadth,0,1,0}; //bottomside edge
	  edge[3] =   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]-layhalflength, PhyVolGlPos[ijk][ij][2],0,0,1}; //frontside edge
		
		  
	}

        
	     layid = ijk+1;
	     layid<<=2;
	     layid+=ij;

	     cout<<"layid: "<<layid<<endl;
	
	     G4ThreeVector extposvec;
	     G4ThreeVector tmphtvec;
	     G4ThreeVector truehtvec;   
        if(pl2){

	  pAnalysis->CMVDpl2[layid]=1;
	  
	  cout<<"---Intersection with Plane found---"<<endl;
	 
	  bool isInside= false;
	  int delta = 0;// 100cm
	  //  cout<<"checkkk: "<<G_Point[2]<<" "<<PhyVolGlPos[ijk][ij][2]<<endl;
	    
	  switch (ijk){
	    //  cout<<"checkkk: "<<G_Point[2]<<PhyVolGlPos[ijk][ij][2]<<endl;
          case 0: isInside = ( (G_Point[0])<  (PhyVolGlPos[ijk][ij][0]+layhalflength+delta) &&  (G_Point[0])>  (PhyVolGlPos[ijk][ij][0]-layhalflength-delta) &&  (G_Point[1])<  (PhyVolGlPos[ijk][ij][1]+layhalfbreadth+delta) &&  (G_Point[1])>  (PhyVolGlPos[ijk][ij][1]-layhalfbreadth-delta)    && abs(G_Point[2]-PhyVolGlPos[ijk][ij][2])<1.e-10 );
	    // cout<<"case 0"<<endl;
	    //	    cout<<fixed<<G_Point[2]<<" "<<fixed<<PhyVolGlPos[ijk][ij][2]<<endl;
	    // cout<<"0"<<PhyVolGlPos[ijk][ij][0]<<" "<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<endl;
	    
	    // cout<<"   (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    "<<    (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   "<<     (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][0])+layhalflength   "<<    (PhyVolGlPos[ijk][ij][0])+layhalflength   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][0])-layhalflength     "<<     (PhyVolGlPos[ijk][ij][0])-layhalflength   <<endl;

   
	    
            break;
	    
          case 1: isInside = (  (G_Point[1])<  (PhyVolGlPos[ijk][ij][1]+layhalfbreadth+delta) &&  (G_Point[1])>  (PhyVolGlPos[ijk][ij][1]-layhalfbreadth-delta) &&  (G_Point[2])<  (PhyVolGlPos[ijk][ij][2]+layhalflength+delta) &&  (G_Point[2])>  (PhyVolGlPos[ijk][ij][2]-layhalflength-delta)       &&     abs(G_Point[0]-PhyVolGlPos[ijk][ij][0])<1.e-10  );
	    //	    cout<<"case 1"<<endl;
	    
	    //	    cout<<"1"<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<endl;
	    
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    "<<    (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   "<<     (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][2])+layhalflength   "<<    (PhyVolGlPos[ijk][ij][2])+layhalflength   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][2])-layhalflength     "<<     (PhyVolGlPos[ijk][ij][2])-layhalflength   <<endl;

   
	    
	    break;
	    
	  case 2: isInside = ( (G_Point[1])<  (PhyVolGlPos[ijk][ij][1]+layhalfbreadth+delta) &&  (G_Point[1])>  (PhyVolGlPos[ijk][ij][1]-layhalfbreadth-delta) &&  (G_Point[2])<  (PhyVolGlPos[ijk][ij][2]+layhalflength+delta) &&  (G_Point[2])>  (PhyVolGlPos[ijk][ij][2]-layhalflength-delta)   &&   abs(G_Point[0] -  PhyVolGlPos[ijk][ij][0])<1.e-10  );
	    //	    cout<<"case 2"<<endl;

	    //	    cout<<"2 "<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<endl;
	    
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    "<<    (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   "<<     (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][2])+layhalflength   "<<    (PhyVolGlPos[ijk][ij][2])+layhalflength   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][2])-layhalflength     "<<     (PhyVolGlPos[ijk][ij][2])-layhalflength   <<endl;

	    break;
	    
	  case 3: isInside = (   (G_Point[0])<  (PhyVolGlPos[ijk][ij][0]+layhalfbreadth+delta) &&  (G_Point[0])>  (PhyVolGlPos[ijk][ij][0]-layhalfbreadth-delta) &&  (G_Point[2])<  (PhyVolGlPos[ijk][ij][2]+layhalflength+delta) &&  (G_Point[2])> (PhyVolGlPos[ijk][ij][2]-layhalflength-delta)     &&  abs(G_Point[1] -PhyVolGlPos[ijk][ij][1])<1.e-10 );
	    //	    cout<<"case 3:: "<<endl;


	    //	    cout<<"3"<<PhyVolGlPos[ijk][ij][0]<<" "<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<" "<<G_Point[0]<<endl;
	    
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][0])+layhalfbreadth    "<<    (PhyVolGlPos[ijk][ij][0])+layhalfbreadth    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][0])-layhalfbreadth   "<<     (PhyVolGlPos[ijk][ij][0])-layhalfbreadth   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][2])+layhalflength   "<<    (PhyVolGlPos[ijk][ij][2])+layhalflength   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][2])-layhalflength     "<<     (PhyVolGlPos[ijk][ij][2])-layhalflength   <<endl;





	    
	    break;



	  case 5: isInside = (  (G_Point[1])<  (PhyVolGlPos[ijk][ij][1]+layhalflength+delta) &&  (G_Point[1])>  (PhyVolGlPos[ijk][ij][1]-layhalflength-delta) &&  (G_Point[2])<  (PhyVolGlPos[ijk][ij][2]+layhalfbreadth+delta) &&  (G_Point[2])>  (PhyVolGlPos[ijk][ij][2]-layhalfbreadth-delta)  &&     abs(G_Point[0]-  PhyVolGlPos[ijk][ij][0])<1.e-10  );
	    //	    cout<<"case 1"<<endl;
	    
	    //	    cout<<"1"<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<endl;
	    
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][1])+layhalflength    "<<    (PhyVolGlPos[ijk][ij][1])+layhalflength    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][1])-layhalflength   "<<     (PhyVolGlPos[ijk][ij][1])-layhalflength   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][2])+layhalfbreadth   "<<    (PhyVolGlPos[ijk][ij][2])+layhalfbreadth   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][2])-layhalfbreadth     "<<     (PhyVolGlPos[ijk][ij][2])-layhalfbreadth   <<endl;

   
	    
	    break;



	  case 6: isInside = (  (G_Point[1])<  (PhyVolGlPos[ijk][ij][1]+layhalflength+delta) &&  (G_Point[1])>  (PhyVolGlPos[ijk][ij][1]-layhalflength-delta) &&  (G_Point[2])<  (PhyVolGlPos[ijk][ij][2]+layhalfbreadth+delta) &&  (G_Point[2])>  (PhyVolGlPos[ijk][ij][2]-layhalfbreadth-delta)       &&     abs(G_Point[0]-  PhyVolGlPos[ijk][ij][0])<1.e-10  );
	    //	    cout<<"case 1"<<endl;
	    
	    //	    cout<<"1"<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<endl;
	    
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][1])+layhalflength    "<<    (PhyVolGlPos[ijk][ij][1])+layhalflength    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][1])-layhalflength   "<<     (PhyVolGlPos[ijk][ij][1])-layhalflength   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][2])+layhalfbreadth   "<<    (PhyVolGlPos[ijk][ij][2])+layhalfbreadth   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][2])-layhalfbreadth     "<<     (PhyVolGlPos[ijk][ij][2])-layhalfbreadth   <<endl;

   
	    
	    break;

	    
          default :isInside=false;
	    //	    cout<<"default case"<<endl;
	    break;
          }

	
	  cout<<"isinside: "<<isInside<<endl;

	
	  // bool hitpresent = false;

	  // layexp->SetExtXPos(G_Point[0]);
	  // layexp->SetExtYPos(G_Point[1]);
	  // layexp->SetExtZPos(G_Point[2]);
	  // layexp->SetUsed(true);
	  // layexp->Print();

 

          if(isInside) {
	    
	    cout<<"Extrapolation inside the detector boundary"<<endl;


	    
	    //error calculation
	    if(ijk==0){ //top wall
	      Ptxerr = sqrt((xxerr*xxerr)+pow((G_Point[2]-posz),2)*txerr*txerr+2*(G_Point[2]-posz)*xxtxerr);
	      Ptyerr = sqrt((yyerr*yyerr)+pow((G_Point[2]-posz),2)*tyerr*tyerr+ 2*(G_Point[2]-posz)*yytyerr);
	      Ptzerr = 0;
		    
	    }
	    else if(ijk==3) {//for back wall
	      Ptxerr =sqrt((xxerr*xxerr)+pow((G_Point[1]-posy)*txerr/dydz,2)+ pow(((G_Point[1]-posy)*dxdz*tyerr)/(dydz*dydz),2)+       pow((dxdz*yyerr)/dydz,2)            + 2*((G_Point[1]-posy)*xxtxerr)/dydz-  2*pow((G_Point[1]-posy),2)*dxdz*txtyerr/pow(dydz,3) - 2*(G_Point[1]-posy)*dxdz*xxtyerr/pow(dydz,2)    - 2*(dxdz*xxtyerr)/dydz - 2*   (G_Point[1]-posy)*dxdz*yytxerr/pow(dydz,2) +2*dxdz*dxdz* (G_Point[1]-posy)  *yytyerr/pow(dydz,3)  );
	      Ptyerr =0 ;
	      Ptzerr =sqrt(  pow((G_Point[1]-posy)*tyerr/(dydz*dydz),2) + pow((yyerr/dydz),2) + 2*(G_Point[1]-posy)*yytyerr/pow(dydz,3)   ) ;
	    }
		  
	    else { //for left right walls
	      Ptxerr =0 ;
	      Ptyerr =sqrt((yyerr*yyerr)+pow((G_Point[0]-posx)*tyerr/dxdz,2)+ pow(((G_Point[0]-posx)*dydz*txerr)/(dxdz*dxdz),2)+pow((dydz*xxerr)/dxdz,2)          + 2*((G_Point[0]-posx)*yytyerr)/dxdz- 2*(dydz*xxtyerr)/dxdz -  2*pow((G_Point[0]-posx),2)*dydz*txtyerr/pow(dxdz,3)-2*((G_Point[0]-posx)*dydz*xxtyerr)/(dxdz*dxdz)+2*(G_Point[0]-posx)*dydz*dydz*xxtxerr/pow(dxdz,3) - 2*(G_Point[0]-posx)*dydz*yytxerr/pow(dxdz,2));
	      Ptzerr =sqrt(pow((G_Point[0]-posx)*txerr/(dxdz*dxdz),2) +pow((xxerr/dxdz),2) +2* ((G_Point[0]-posx)*xxtxerr)/pow(dxdz,3)) ;
	    }
	    // cout<<"Ptxerr "<<Ptxerr<<" Ptyerr "<<Ptyerr<<"Ptzerr "<<Ptzerr<<endl;

		  

	   
	     cout<<"Point inside the detector boundary"<<endl;

	   


	    extposvec.setX(G_Point[0]);
	    extposvec.setY(G_Point[1]);
	    extposvec.setZ(G_Point[2]);
 
	    
	     cout<<"Point:"<<extposvec<<endl;
	     cout<<" CmvClusterBank[ijk][ij].size() "<<CmvClusterBank[ijk][ij].size()<<endl;
	    if(CmvClusterBank[ijk][ij].size()>0){
	      cout<<"Hit Present"<<endl;
            for (unsigned int ix=0; ix<CmvClusterBank[ijk][ij].size(); ix++) {

	      cout<<"## "<<CmvClusterBank[ijk][ij][ix]->GetPlane()-1<<" "<<CmvClusterBank[ijk][ij][ix]->GetLayer()<<" "<<ijk<<" "<<ij<<" "<<ix<<endl;
	      
                
        
		CmvClusterBank[ijk][ij][ix]->Print();


	     	
		//		cout<<"CmvClusterBank[ijk][ij][ix] inside loop "<<CmvClusterBank[ijk][ij][ix]->GetPlane()<<" "<<CmvClusterBank[ijk][ij][ix]->GetLayer()<<endl;
		tmphtvec.setX(CmvClusterBank[ijk][ij][ix]->GetRecoPosX());
		tmphtvec.setY(CmvClusterBank[ijk][ij][ix]->GetRecoPosY());
		tmphtvec.setZ(CmvClusterBank[ijk][ij][ix]->GetRecoPosZ());

		//		cout<<"recohtvec: "<<tmphtvec<<endl;
        	double difx = (extposvec-tmphtvec).mag();

		truehtvec.setX(CmvClusterBank[ijk][ij][ix]->GetTruePosX());
		truehtvec.setY(CmvClusterBank[ijk][ij][ix]->GetTruePosY());
		truehtvec.setZ(CmvClusterBank[ijk][ij][ix]->GetTruePosZ());
		//		cout<<" truehtvvec: "<<truehtvec<<endl;
		//Store the differences in three


		//For multiple clusters in a layer consider one with minimum difx
		if (difx < diffmx) {

		  diffmx = difx;
		  //		  cout<<" Diff "<<difx<<endl;
                
	
     		

		  cout<<" Extrapolated position "<<G_Point[0]<<" "<<G_Point[1]<<" "<<G_Point[2]<<endl;
		  cout<<" Reconstructed position "<<tmphtvec.x()<<" "<<tmphtvec.y()<<" "<<tmphtvec.z()<<endl;
		  cout<<" True Position "<<truehtvec.x()<<" "<<truehtvec.y()<<" "<<truehtvec.z()<<endl;

		  cout<<" Reco-Extrap: "<<  tmphtvec.x()-G_Point[0]<<" "<<  tmphtvec.y()-G_Point[1]<<" "    <<tmphtvec.z()-G_Point[2]<<endl;
		  cout<<" True-Extra: "<<truehtvec.x()-G_Point[0]<<" " << truehtvec.y()-G_Point[1]<<" " << truehtvec.z()-G_Point[2] <<endl;
    


		  
		  pAnalysis->CMVDTruePosX[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosX();
		  pAnalysis->CMVDTruePosY[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosY();
		  pAnalysis->CMVDTruePosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosZ();

		  
		  pAnalysis->CMVDRecoPosX[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosX();
		  pAnalysis->CMVDRecoPosY[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosY();
		  pAnalysis->CMVDRecoPosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosZ();


		  pAnalysis->CMVDWRecoPosX[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosX();
		  pAnalysis->CMVDWRecoPosY[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosY();
		  pAnalysis->CMVDWRecoPosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosZ();

		  pAnalysis->CMVDClusterSize[layid] = CmvClusterBank[ijk][ij][ix]->GetClusterSize();
		   
		  pAnalysis->CMVDExpnHit[layid] = 3;		  
		  pAnalysis->CMVDisExpIn[layid] =1;
                  pAnalysis->CMVDisHit[layid] = 1;
	  
	  
		}//difx
		
	      

		
            }//for
	    }//if
	   
	    else{
	      cout<<"Extrapolation is there but no hit"<<endl;
	      pAnalysis->CMVDExpnHit[layid] = 2;

	      pAnalysis->CMVDisExpIn[layid] =1;
              pAnalysis->CMVDisHit[layid] =0;

	      pAnalysis->CMVDTruePosX[layid] = -10000;
	      pAnalysis->CMVDTruePosY[layid] = -10000;
	      pAnalysis->CMVDTruePosZ[layid] = -10000;
	      
	      
	      pAnalysis->CMVDRecoPosX[layid] =  -10000;
	      pAnalysis->CMVDRecoPosY[layid] =  -10000;
	      pAnalysis->CMVDRecoPosZ[layid] =  -10000;
	      
	      
	      pAnalysis->CMVDWRecoPosX[layid] =  -10000;
	      pAnalysis->CMVDWRecoPosY[layid] =  -10000;
	      pAnalysis->CMVDWRecoPosZ[layid] =  -10000;
	      
	      pAnalysis->CMVDClusterSize[layid] =0;
	      




	      
	    }
          }//isinside


	  else{
	    cout<<"..Extrapolated point outside.."<<endl;

	    if(CmvClusterBank[ijk][ij].size()>0){

	      cout<<"But Hit Present"<<endl;
	      for (unsigned int ix=0; ix<CmvClusterBank[ijk][ij].size(); ix++) {
		
		//		cout<<"## "<<CmvClusterBank[ijk][ij][ix]->GetPlane()-1<<" "<<CmvClusterBank[ijk][ij][ix]->GetLayer()<<endl;
        
		CmvClusterBank[ijk][ij][ix]->Print();
	     	
		//		cout<<"CmvClusterBank[ijk][ij] inside loop "<<CmvClusterBank[ijk][ij][ix]->GetPlane()<<" "<<CmvClusterBank[ijk][ij][ix]->GetLayer()<<endl;
		tmphtvec.setX(CmvClusterBank[ijk][ij][ix]->GetRecoPosX());
		tmphtvec.setY(CmvClusterBank[ijk][ij][ix]->GetRecoPosY());
		tmphtvec.setZ(CmvClusterBank[ijk][ij][ix]->GetRecoPosZ());

		//		cout<<"recohtvec: "<<tmphtvec<<endl;
        	 double difx = (extposvec-tmphtvec).mag();

		truehtvec.setX(CmvClusterBank[ijk][ij][ix]->GetTruePosX());
		truehtvec.setY(CmvClusterBank[ijk][ij][ix]->GetTruePosY());
		truehtvec.setZ(CmvClusterBank[ijk][ij][ix]->GetTruePosZ());
		cout<<" truehtvvec: "<<truehtvec<<endl;
		//Store the differences in three


		//For multiple clusters in a layer consider one with minimum difx
		if (difx < diffmx) {

		  diffmx = difx;
		  //cout<<" Diff "<<difx<<endl;
                
	
     		

		  cout<<" Extrapolated position "<<G_Point[0]<<" "<<G_Point[1]<<" "<<G_Point[2]<<endl;
		  cout<<" Reconstructed position "<<tmphtvec.x()<<" "<<tmphtvec.y()<<" "<<tmphtvec.z()<<endl;
		  cout<<" True Position "<<truehtvec.x()<<" "<<truehtvec.x()<<" "<<truehtvec.y()<<" "<<truehtvec.z()<<endl;

		  cout<<"Reco-Extrap: "<<  tmphtvec.x()-G_Point[0]<<" "<<  tmphtvec.y()-G_Point[1]<<" "    <<tmphtvec.z()-G_Point[2]<<endl;
		  cout<<"True-Extra: "<<truehtvec.x()-G_Point[0]<<" " << truehtvec.y()-G_Point[1]<<" " << truehtvec.z()-G_Point[2] <<endl;
    


		  
		  pAnalysis->CMVDTruePosX[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosX();
		  pAnalysis->CMVDTruePosY[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosY();
		  pAnalysis->CMVDTruePosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosZ();

		  
		  pAnalysis->CMVDRecoPosX[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosX();
		  pAnalysis->CMVDRecoPosY[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosY();
		  pAnalysis->CMVDRecoPosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosZ();


		  pAnalysis->CMVDWRecoPosX[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosX();
		  pAnalysis->CMVDWRecoPosY[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosY();
		  pAnalysis->CMVDWRecoPosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosZ();

		  pAnalysis->CMVDClusterSize[layid] = CmvClusterBank[ijk][ij][ix]->GetClusterSize();
		   
		  pAnalysis->CMVDExpnHit[layid] = 1;		  
		  pAnalysis->CMVDisExpIn[layid] =0;
		  pAnalysis->CMVDisHit[layid] =1;

	  
	  
		}//difx
		
	      

		
	      }//  for (unsigned int ix=0; ix<CmvClusterBank[ijk][ij].size(); ix++) {
	    }//    if(CmvClusterBank[ijk][ij].size()>0){
	    else{
	      cout<<"Extrapolation is outside and no hits"<<endl;
	    pAnalysis->CMVDExpnHit[layid] = 0;// exp point is outside and no hit is there
              pAnalysis->CMVDisExpIn[layid] =0;
              pAnalysis->CMVDisHit[layid] =0;


	      pAnalysis->CMVDTruePosX[layid] = -10000;
	      pAnalysis->CMVDTruePosY[layid] = -10000;
	      pAnalysis->CMVDTruePosZ[layid] = -10000;
	      
	      
	      pAnalysis->CMVDRecoPosX[layid] =  -10000;
	      pAnalysis->CMVDRecoPosY[layid] =  -10000;
	      pAnalysis->CMVDRecoPosZ[layid] =  -10000;
	      
	      
	      pAnalysis->CMVDWRecoPosX[layid] =  -10000;
	      pAnalysis->CMVDWRecoPosY[layid] =  -10000;
	      pAnalysis->CMVDWRecoPosZ[layid] =  -10000;
	      
	      pAnalysis->CMVDClusterSize[layid] =0;

	    
	    }

	  }


		  pAnalysis->CMVDExpPosX[layid] = G_Point[0];
		  pAnalysis->CMVDExpPosY[layid] = G_Point[1];
		  pAnalysis->CMVDExpPosZ[layid] = G_Point[2];

        

	}	else{//pl2
	  cout<<"..No Line PLane Intersection.."<<endl;
	  
   if(CmvClusterBank[ijk][ij].size()>0){

	      cout<<"But Hit Present"<<endl;
	      for (unsigned int ix=0; ix<CmvClusterBank[ijk][ij].size(); ix++) {
		
		//		cout<<"## "<<CmvClusterBank[ijk][ij][ix]->GetPlane()-1<<" "<<CmvClusterBank[ijk][ij][ix]->GetLayer()<<endl;
	      
                
        
		CmvClusterBank[ijk][ij][ix]->Print();


	     	
		//		cout<<"CmvClusterBank[ijk][ij] inside loop "<<CmvClusterBank[ijk][ij][ix]->GetPlane()<<" "<<CmvClusterBank[ijk][ij][ix]->GetLayer()<<endl;
		tmphtvec.setX(CmvClusterBank[ijk][ij][ix]->GetRecoPosX());
		tmphtvec.setY(CmvClusterBank[ijk][ij][ix]->GetRecoPosY());
		tmphtvec.setZ(CmvClusterBank[ijk][ij][ix]->GetRecoPosZ());

		//		cout<<"recohtvec: "<<tmphtvec<<endl;
        	 double difx = (extposvec-tmphtvec).mag();

		truehtvec.setX(CmvClusterBank[ijk][ij][ix]->GetTruePosX());
		truehtvec.setY(CmvClusterBank[ijk][ij][ix]->GetTruePosY());
		truehtvec.setZ(CmvClusterBank[ijk][ij][ix]->GetTruePosZ());
		cout<<" truehtvvec: "<<truehtvec<<endl;
		//Store the differences in three


		//For multiple clusters in a layer consider one with minimum difx
		if (difx < diffmx) {

		  diffmx = difx;
		  //cout<<" Diff "<<difx<<endl;
                
	
     		

		  cout<<" Extrapolated position "<<G_Point[0]<<" "<<G_Point[1]<<" "<<G_Point[2]<<endl;
		  cout<<" Reconstructed position "<<tmphtvec.x()<<" "<<tmphtvec.y()<<" "<<tmphtvec.z()<<endl;
		  cout<<" True Position "<<truehtvec.x()<<" "<<truehtvec.x()<<" "<<truehtvec.y()<<" "<<truehtvec.z()<<endl;

		  cout<<"Reco-Extrap: "<<  tmphtvec.x()-G_Point[0]<<" "<<  tmphtvec.y()-G_Point[1]<<" "    <<tmphtvec.z()-G_Point[2]<<endl;
		  cout<<"True-Extra: "<<truehtvec.x()-G_Point[0]<<" " << truehtvec.y()-G_Point[1]<<" " << truehtvec.z()-G_Point[2] <<endl;
    


		  
		  pAnalysis->CMVDTruePosX[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosX();
		  pAnalysis->CMVDTruePosY[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosY();
		  pAnalysis->CMVDTruePosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosZ();

		  
		  pAnalysis->CMVDRecoPosX[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosX();
		  pAnalysis->CMVDRecoPosY[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosY();
		  pAnalysis->CMVDRecoPosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosZ();


		  pAnalysis->CMVDWRecoPosX[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosX();
		  pAnalysis->CMVDWRecoPosY[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosY();
		  pAnalysis->CMVDWRecoPosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosZ();

		  pAnalysis->CMVDClusterSize[layid] = CmvClusterBank[ijk][ij][ix]->GetClusterSize();
		   
		  pAnalysis->CMVDExpnHit[layid] = 1;		  
		 		  
	  
		}//difx
		
	      

		
	      }
   }  else{
		    cout<<"no lineplane intersection and no hits"<<endl;


	      pAnalysis->CMVDTruePosX[layid] = -10000;
	      pAnalysis->CMVDTruePosY[layid] = -10000;
	      pAnalysis->CMVDTruePosZ[layid] = -10000;
	      
	      
	      pAnalysis->CMVDRecoPosX[layid] =  -10000;
	      pAnalysis->CMVDRecoPosY[layid] =  -10000;
	      pAnalysis->CMVDRecoPosZ[layid] =  -10000;
	      
	      
	      pAnalysis->CMVDWRecoPosX[layid] =  -10000;
	      pAnalysis->CMVDWRecoPosY[layid] =  -10000;
	      pAnalysis->CMVDWRecoPosZ[layid] =  -10000;
	      
	      pAnalysis->CMVDClusterSize[layid] =0;
	      
        

		  pAnalysis->CMVDClusterSize[layid] = 0;
		   
		  pAnalysis->CMVDExpnHit[layid] = 0;


		  }
		  
        	  pAnalysis->CMVDpl2[layid]=0;
		  pAnalysis->CMVDExpPosX[layid] = -10000;
		  pAnalysis->CMVDExpPosY[layid] = -10000;
		  pAnalysis->CMVDExpPosZ[layid] = -10000;		  
		 		  

	}

	  
      }//layer loop
  
  
    }//loc_no loop

    cout<<"Printing CMVDExpnHit "<<endl;
    for(int iij=0;iij<32;iij++){
      cout<<iij<<" "<<pAnalysis->CMVDExpnHit[iij]<<endl;
  }



    
  }


	
}







void CMVDRecoAlg::CreateCmvHit() {
  //   SipmHit_pointer = SipmHit_Manager::APointer; 
  CmvHit_pointer = new CmvHit_Manager();

  CmvHit_pointer->CmvHit_list.clear(); 

  // Convert SiPM hit to Cmv Hit
  // SipmHit* foursipm[4]={0}; //GMA memory leakage ?
  cout<<SipmHit_pointer<<endl;
  cout <<"micalEventAction::CreateCmvHit() "<<pAnalysis->isInOut<<" "<<SipmHit_pointer->SipmHit_list.size()<<endl;


    cout<<" SipmHit_pointer->SipmHit_list.size() "<<SipmHit_pointer->SipmHit_list.size()<<endl;
    for (unsigned int ijj=0; ijj<SipmHit_pointer->SipmHit_list.size(); ijj++) {
      SipmHit_pointer->SipmHit_list[ijj]->Print();  
    }
   
    
    for (unsigned int ij=0; ij<SipmHit_pointer->SipmHit_list.size(); ij++) {
      //  if(SipmHit_pointer->SipmHit_list[ij]->GetPulse()<160) continue;
      SipmHit* foursipm[4]={0}; //GMA memory leakage ?

      int isfoursipm[4]={0};
      int tmpstripid = -1;
      int tmpside = -1; //Used this to find global position of the layer
      int tmplayer=-1;
      if (!(SipmHit_pointer->SipmHit_list[ij]->GetUsed())) {
	tmpstripid = SipmHit_pointer->SipmHit_list[ij]->GetStripId();
	//	cout <<"tmp "<< tmpstripid<<" "<<int(SipmHit_pointer->SipmHit_list[ij]->GetUsed())<<" "<<SipmHit_pointer->SipmHit_list[ij]->GetSiPM()<<endl;
        int isipm = SipmHit_pointer->SipmHit_list[ij]->GetSiPM();

	isfoursipm[isipm]=1;
	foursipm[isipm] = SipmHit_pointer->SipmHit_list[ij];
	SipmHit_pointer->SipmHit_list[ij]->SetUsed(true);
				
	tmpside = SipmHit_pointer->SipmHit_list[ij]->GetPlane()-1; //We had added 1 while storing it.
	tmplayer = SipmHit_pointer->SipmHit_list[ij]->GetLayer();

	//	cout<<tmpside+1<<" "<<tmplayer<<" "<<isipm<<endl;
	
	//Look for all other SiPM of same strip
	for (unsigned int jk=ij+1; jk<SipmHit_pointer->SipmHit_list.size(); jk++) {
	  // if(SipmHit_pointer->SipmHit_list[jk]->GetPulse()<160) continue;
	  if (!(SipmHit_pointer->SipmHit_list[jk]->GetUsed())) {	
	    int tmpstripid2 = SipmHit_pointer->SipmHit_list[jk]->GetStripId();
	    //  cout <<"tmp2 "<< tmpstripid2 <<endl;
	    int isipm2 = SipmHit_pointer->SipmHit_list[jk]->GetSiPM();
	    int tmpside2 = SipmHit_pointer->SipmHit_list[jk]->GetPlane()-1;
	    int tmplayer2 = SipmHit_pointer->SipmHit_list[jk]->GetLayer();

	    //	    cout<<tmpside2+1<<" "<<tmplayer2<<" "<<isipm2<<endl;
	    if (tmpstripid !=tmpstripid2) continue; //this continues to next iteration
	    //	cout<<tmpside2+1<<" "<<tmplayer2<<" "<<isipm2<<endl;
	    ////	if(tmpside !=tmpside2 || tmplayer !=tmplayer2) continue;
        
            isfoursipm[isipm2]=1;
	    foursipm[isipm2] = SipmHit_pointer->SipmHit_list[jk];
	    SipmHit_pointer->SipmHit_list[jk]->SetUsed(true);
	  }
	}
      }
      //      cout<<foursipm[0]<<" "<< foursipm[1]<< " "<< foursipm[2]<<" "<< foursipm[3]<<endl;
      //      cout<<"isfoursipm: "<< isfoursipm[0]<<" "<<isfoursipm[1]<<" "<<isfoursipm[2]<<" "<<isfoursipm[3]<<endl;
      if (tmpside>=0 && isfoursipm[0]+isfoursipm[1]+isfoursipm[2]+isfoursipm[3]>1) { //# atleast 2 sipms must have hit
	double pos[3];
	//	cout<<"posvec: ";
	for (int ix=0; ix<3; ix++) {pos[ix] = PhyVolGlPos[tmpside][tmplayer][ix];cout<<pos[ix]<<" ";}
	//	cout<<endl;

	//	cout<<"creating cmv hit "<<endl;
	if(foursipm[0]){	cout<<foursipm[0]->GetStripId()<<" "<<foursipm[0]->GetSiPM()<< endl;}
	if(foursipm[1]){	cout<<foursipm[1]->GetStripId()<<" "<<foursipm[1]->GetSiPM()<< endl;}
	if(foursipm[2]){	cout<<foursipm[2]->GetStripId()<<" "<<foursipm[2]->GetSiPM()<< endl;}
	if(foursipm[3]){	cout<<foursipm[3]->GetStripId()<<" "<<foursipm[3]->GetSiPM()<< endl;}
	CmvHit* tmpcmvHit = new CmvHit(foursipm[0], foursipm[1], foursipm[2], foursipm[3], pos);
	//	tmpcmvHit->Print();		     
	CmvHit_pointer->CmvHit_list.push_back(tmpcmvHit);
      }
    }

      

  //Calling Save cmv hit data here
  //  if (pAnalysis->cmv_nhit >pAnalysis->cmv_nhtmx) pAnalysis->cmv_nhit =pAnalysis->cmv_nhtmx;

    CmvHit_Manager *pCmvHit = CmvHit_Manager::APointer;
  if(pCmvHit){
    
    if(pCmvHit->CmvHit_list.size()) {
      
      pAnalysis->pRootFile->cd();
      pAnalysis->cmv_nhit = pCmvHit->CmvHit_list.size();
      for (unsigned int ij=0; ij<pCmvHit->CmvHit_list.size() && ij<pAnalysis->cmv_nhit; ij++) {
	
	SaveCMVHitData(ij,pCmvHit->CmvHit_list[ij]);
	
      }
    }
  }

  
}  

// cmv hit



void CMVDRecoAlg::FormCmvCluster() {

  cout<<"Form Cmv_Cluster"<<endl;
 CmvCluster_pointer = new CmvCluster_Manager();
  CmvCluster_pointer->CmvCluster_list.clear(); 

  // Form cmv cluster

  cout <<"micalEventAction::CreateCmvHit() "<<pAnalysis->isInOut<<" "<<CmvHit_pointer->CmvHit_list.size()<<endl<<endl;
  // if (pAnalysis->isInOut==0 || pAnalysis->isInOut ==3 || pAnalysis->isInOut==5) {//0:GEN->RECO, 1:GEN->DIGI, 2:GEN->SIM, 3: SIM -> RECO, 4: SIM -> DIFI, 5 : DIGI -> RECO
    cout<<"CmvHit_list size: "<<CmvHit_pointer->CmvHit_list.size()<<endl;


    cout<<"CmvHit_pointer->CmvHit_list.size() "<<CmvHit_pointer->CmvHit_list.size()<<endl;
    for (unsigned int ijj=0; ijj<CmvHit_pointer->CmvHit_list.size(); ijj++) {
      CmvHit_pointer->CmvHit_list[ijj]->Print();
    }

  cout<<"...........First Sorting hits acc to stripno in each layer..................."<<endl;
  //Note: Here we compare strip numbers which ranges from 0-88/40 in each layers rather than stripid.. since for two nearby strips strip no. diff is 1 but stripId diff is not one.
  
    for (unsigned int ix=0; ix<CmvHit_pointer->CmvHit_list.size(); ix++) {
    
      int side1 = CmvHit_pointer->CmvHit_list[ix]->GetPlane()-1;
      int lay1 = CmvHit_pointer->CmvHit_list[ix]->GetLayer();
      int strp1id  = CmvHit_pointer->CmvHit_list[ix]->GetStripId();
      cout<<ix<<" "<<strp1id<<endl;
      
      for (unsigned int ixi=ix+1; ixi<CmvHit_pointer->CmvHit_list.size(); ixi++) {
      int side2 = CmvHit_pointer->CmvHit_list[ixi]->GetPlane()-1;
      int lay2 = CmvHit_pointer->CmvHit_list[ixi]->GetLayer();
      int strp2id  = CmvHit_pointer->CmvHit_list[ixi]->GetStripId(); // this id is full id plus right shifting 2 bits i.e. excluding sipm 2-bits
  cout<<ixi<<" "<<strp2id<<endl;
        
      //	if(side1==side2 && lay1==lay2){

	  
  //  if(strp1id<strp2id){
  if(CmvHit_pointer->CmvHit_list[ixi]->GetStripId()<CmvHit_pointer->CmvHit_list[ix]->GetStripId()){
    swap(CmvHit_pointer->CmvHit_list[ix],CmvHit_pointer->CmvHit_list[ixi]);
    cout<<"swap"<<endl;
	  }//  if(strp1<strp2){
	
	  //	}// 	if(side1==side2 && lay1==lay2){

    	
      }// for (unsigned ixi=ix+1; ixi<CmvHit_pointer->CmvHit_list.size(); ixi++) {
    
    }//  for (unsigned ix=0; ix<CmvHit_pointer->CmvHit_list.size(); ix++) {


    cout<<"................sorted array: ..............."<<endl;
    for (unsigned int ix=0; ix<CmvHit_pointer->CmvHit_list.size(); ix++) {
      CmvHit* cmvhit = CmvHit_pointer->CmvHit_list[ix];
      cmvhit->SetUsed(false);
      cmvhit->Print();
      CmvHitBank[cmvhit->GetPlane()-1][cmvhit->GetLayer()].push_back(cmvhit);
    }//  for (unsigned int ix=0; ix<CmvHit_pointer->CmvHit_list.size(); ix++) {
  
    //Fill number of hits in a layer:

    
    for(unsigned int tmpside=0;tmpside<7;tmpside++){
      
      for(unsigned int tmplay=0;tmplay<4;tmplay++){
	if(tmpside>0 && tmplay==3) continue;
	//	cout<<"tmpside tmplay hits "<<tmpside<<" "<<tmplay <<" "<<CmvHitBank[tmpside][tmplay].size()<< endl;   
    for(unsigned int jk=0; jk<CmvHitBank[tmpside][tmplay].size();jk++) {
      // cout<<"jk "<<jk<<endl;

      	if(!(CmvHitBank[tmpside][tmplay][jk]->GetUsed())) {

	  CmvCluster* tmpclust = new CmvCluster(CmvHitBank[tmpside][tmplay][jk]);

	  //Change Used from 0 to 1 to show that hit has been added
	  CmvHitBank[tmpside][tmplay][jk]->SetUsed(true);

	  //loop over other entries in same side and same lay
	  for(unsigned int kl=jk+1; kl<CmvHitBank[tmpside][tmplay].size(); kl++) {
	    //	    cout<<"kl "<<kl<<endl;

	    if( abs(CmvHitBank[tmpside][tmplay][kl-1]->GetStrip() - CmvHitBank[tmpside][tmplay][kl]->GetStrip()) >=2 ) break;//Geometrically nearby  xx x xx 

	    
	    if(  !( CmvHitBank[tmpside][tmplay][kl]->GetUsed())  ){

	      tmpclust->AddHits(CmvHitBank[tmpside][tmplay][kl]);

	      
	      CmvHitBank[tmpside][tmplay][kl]->SetUsed(true);

	      
	    }//  if(  !( CmvHitBank[tmpside][tmplay][kl]->GetUsed())  ){
	  
	    
	    
	  }//	  for(unsigned int kl=jk+1; kl<CmvHitBank[tmpside][tmplay].size(); kl++) {


	  CmvCluster_pointer->CmvCluster_list.push_back(tmpclust);
	  
	}//  	if(!(CmvHitBank[tmpside][tmplay][jk]->GetUsed())) {




    }//  for(unsigned int jk=0; jk<CmvHitBank[tmpside][tmplay].size();jk++) {
      }// for(unsigned int tmplay=0;tmplay<4;tmplay++){
      
    }// for(unsigned int tmpside=0;tmpside<7;tmpside++){
    

    cout<<" CmvCluster_list.size() "<<CmvCluster_pointer->CmvCluster_list.size()<<endl;
    for (unsigned int ix=0; ix<CmvCluster_pointer->CmvCluster_list.size(); ix++) {
      CmvCluster* cmvcluster = CmvCluster_pointer->CmvCluster_list[ix];  
      cmvcluster->Print();
    }//  for (unsigned int ix=0; ix<CmvCluster_pointer->CmvCluster_list.size(); ix++) {
    










    

  //Calling Save cluster  data here
   CmvCluster_Manager *pCmvClust = CmvCluster_Manager::APointer;
  if(pCmvClust){
    
    if(pCmvClust->CmvCluster_list.size()) {
      
      pAnalysis->pRootFile->cd();
      pAnalysis->cmv_nclust = pCmvClust->CmvCluster_list.size();
      for (unsigned int ij=0; ij<pCmvClust->CmvCluster_list.size(); ij++) {
	
	SaveCMVClusterData(ij,pCmvClust->CmvCluster_list[ij]);
	
      }
    }
  }



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CMVDRecoAlg::SaveCMVHitData(unsigned int hits,CmvHit* pCmvHit ) {

 
  //  cout<<"  pCmvHit->CmvHit_list.size():  "<<  pCmvHit->CmvHit_list.size()<<endl;
  //  cout <<"  pAnalysis->cmv_nclusthit "<<  pAnalysis->cmv_nhit<<endl;
  
  unsigned int hitId =0;
  
    hitId = pCmvHit->GetPlane();
    hitId<<=2;
    hitId+=pCmvHit->GetLayer();
    hitId<<=7;
    hitId+=pCmvHit->GetStrip();
    hitId<<=2; //just shifted by 2 as we did this while forming strip	
    //pAnalysis->cmv_clustpdgid[hits] =pCmvHit->CmvHit_list[hits]->GetpdgId();
    
    cout<<"clustid "<<hitId<<" "<<endl;
    
    pAnalysis->cmv_hitid[hits] =hitId;
    
    cout<< pAnalysis->cmv_hitid[hits]<<endl;
    
    pAnalysis->cmv_hitLeTim[hits] =pCmvHit->GetLeTime();
    pAnalysis->cmv_hitRiTim[hits] =pCmvHit->GetRiTime();
    
    pAnalysis->cmv_hitLePul[hits] =pCmvHit->GetLePulse();
    pAnalysis->cmv_hitRiPul[hits] =pCmvHit->GetRiPulse();
    
    
    pAnalysis->cmv_hitTrueposx[hits] =pCmvHit->GetTruePosX();
    pAnalysis->cmv_hitTrueposy[hits] =pCmvHit->GetTruePosY();
    pAnalysis->cmv_hitTrueposz[hits] =pCmvHit->GetTruePosZ();
    
    pAnalysis->cmv_hitRecoposx[hits] =pCmvHit->GetRecoPosX();
    pAnalysis->cmv_hitRecoposy[hits] =pCmvHit->GetRecoPosY();
    pAnalysis->cmv_hitRecoposz[hits] =pCmvHit->GetRecoPosZ();
    
    //	pAnalysis->cmv_clustsiz[hits] =pCmvHit->GetHitsize();
    // pAnalysis->cmv_clustmom[hits] =pCmvHit->GetSimMom();
    // pAnalysis->cmv_clustthe[hits] =pCmvHit->GetSimThe();
    // pAnalysis->cmv_clustphi[hits] =pCmvHit->GetSimPhi();
      
  //			pAnalysis->pEventTree->Fill();
    
}


//..................

void CMVDRecoAlg::SaveCMVClusterData(unsigned int clusthits,CmvCluster* pCmvCluster ) {
 
  //  cout<<"  pCmvCluster->CmvCluster_list.size():  "<<  pCmvCluster->CmvCluster_list.size()<<endl;
  //  cout <<"  pAnalysis->cmv_nclusthit "<<  pAnalysis->cmv_nclust<<endl;
  
  unsigned int clustId =0;
  
    clustId = pCmvCluster->GetPlane();

    clustId<<=2;
    clustId+=pCmvCluster->GetLayer();
    clustId<<=7;
    clustId+=pCmvCluster->GetStrip();
    clustId<<=2; //just shifted by 2 as we did this while forming strip	
    //pAnalysis->cmv_clustpdgid[clusthits] =pCmvCluster->GetpdgId();
    
    cout<<"clustid "<<clustId<<" "<<endl;
    
    pAnalysis->cmv_clustid[clusthits] =clustId;
    
    cout<< pAnalysis->cmv_clustid[clusthits]<<endl;
    
    // pAnalysis->cmv_clustLeTim[clusthits] =pCmvCluster->GetLeTime();
    // pAnalysis->cmv_clustRiTim[clusthits] =pCmvCluster->GetRiTime();
    
    // pAnalysis->cmv_clustLePul[clusthits] =pCmvCluster->GetLePulse();
    // pAnalysis->cmv_clustRiPul[clusthits] =pCmvCluster->GetRiPulse();
    
    
    pAnalysis->cmv_clustTrueposx[clusthits] =pCmvCluster->GetTruePosX();
    pAnalysis->cmv_clustTrueposy[clusthits] =pCmvCluster->GetTruePosY();
    pAnalysis->cmv_clustTrueposz[clusthits] =pCmvCluster->GetTruePosZ();
    
    pAnalysis->cmv_clustRecoposx[clusthits] =pCmvCluster->GetRecoPosX();
    pAnalysis->cmv_clustRecoposy[clusthits] =pCmvCluster->GetRecoPosY();
    pAnalysis->cmv_clustRecoposz[clusthits] =pCmvCluster->GetRecoPosZ();


       
    	pAnalysis->cmv_clustsiz[clusthits] =pCmvCluster->GetClusterSize();
    // pAnalysis->cmv_clustmom[clusthits] =pCmvCluster->GetSimMom();
    // pAnalysis->cmv_clustthe[clusthits] =pCmvCluster->GetSimThe();
    // pAnalysis->cmv_clustphi[clusthits] =pCmvCluster->GetSimPhi();
    
    
   
  
  //			pAnalysis->pEventTree->Fill();
  


  
}



//




