TString outName = "efficiency2.root";

void anaEmbedClusters();
void drawEmbedClusters();
void anaResiduals();

#include "myAnalysis.h"
#include "ImageProcessing.h"
using namespace ImageProcessing;

std::vector<Int_t> types = {0};
std::vector<Int_t> parts = {0,1,2,3,4};
std::vector<TString> typeNames = {"embed132","embed108"};
std::vector<TString> partNames = {"proton","deuteron","triton","he3","alpha"};

std::vector<TString> recoTypes = {"va","reco","vapri"};

const int effBins[5] = { 150, 150, 150, 150, 150 };
const double effXmin[5] = { 100, 200, 100, 100, 100 };
const double effXmax[5] = { 1500, 2200, 3200, 3200, 4200 };


TString embedFileName(int type, int part)
{
	TString path = "/data/Q23393/common/efficiency/embed108/";
	//TString path = "/data/Q18393/isobe/20190822SingleIonMC/20200115/tree_proton_embed132.left.vaprip.root";
	//TString path = "/data/Q18393/common/efficiency/";
	//TString path = "./test/";

	//return path+"file_"+partNames[part]+".root";
	return path+"tree_"+partNames[part]+"_embed108.root";
}
//TString outName = "~/data/EmbedAna.20200115.left.vaprip.root";

TString nameId(int type, int pid=-1, int rtype=-1)
{
	TString name = Form("_Embed%dSn",beamA[type]);
	if(pid!=-1) name += "_"+pidName[pid];
	if(rtype!=-1) name += "_"+recoTypes[rtype];
	return name;
}

std::vector<int> multBins={0};

const Double_t invProjA[3] = {44.1/1000., 55.2/1000.,45./1000. };
const Double_t beta[3]  = {-0.364, -0.349, -0.355 };
const Double_t ycoll[3] = {0.382, 0.365, 0.371 };
const Double_t ucoll[] = {0.392, 0.373, 0.371 };
const Double_t momCut[]  = { 2000., 4000., 4200. };

const double yAA[]   = {0.3822, 0.3647, 0.3538, 0.3902};
const double yNN[]   = {0.3696, 0.3697, 0.3705, 0.3706};
const double yBeam[] = {0.7421, 0.7423, 0.7439, 0.7441};
const double uNN[]   = {0.3709, 0.3709, 0.3718, 0.3719};

const Int_t nType = 1;
const Int_t nPart = 5;
const Int_t nReco = 3;
const Int_t nCent = 5; // same as FOPI

TFile *effFile[nType][nPart];
TTree *mctrkTree[nType][nPart], *trkTree[nType][nPart], *eveTree[nType][nPart];
TEntryList *mcelist[nType][nPart], *trkelist[nType][nPart], *eveelist[nType][nPart];
Bool_t   isEffLoad=kFALSE;
Double_t zet;
//,aoq;
Double_t vtxx,vtxy,vtxz;
Double_t dist,recodist;
Double_t theta;
Double_t px,py,pz;
Double_t vapripx,vapripy,vapripz;
Double_t recopx,recopy,recopz;
Double_t recophi,recotheta;
Double_t mcpx,mcpy,mcpz;
Double_t mctheta;
Int_t    nclus, nlclus, nrclus, neclus;
Int_t    ntrk,ngrtrk;
Int_t    ggclose;

void fillEmbedding();
//void drawEmbedding();

void anaEmbedding()
{
	fillEmbedding();
	//drawEmbedding();
}

void fillEmbedding()
{
	for(auto i: types)for(auto j: parts){
		effFile[i][j] = TFile::Open(embedFileName(i,j));
        cout << effFile[i][j]->GetName() << endl;
		effFile[i][j]->GetObject("mctrktree",mctrkTree[i][j]);
		effFile[i][j]->GetObject("trktree",trkTree[i][j]);
		effFile[i][j]->GetObject("evetree",eveTree[i][j]);
	}

	auto setupTrkTree = [&](TTree *t)
	{
		t->SetBranchAddress("zet",&zet);
		t->SetBranchAddress("aoq",&aoq);
		 
		t->SetBranchAddress("vtxx",&vtxx);
		t->SetBranchAddress("vtxy",&vtxy);
		t->SetBranchAddress("vtxz",&vtxz);
		t->SetBranchAddress("dist",&dist);
		t->SetBranchAddress("theta",&theta);
		
		t->SetBranchAddress("px",&px);
		t->SetBranchAddress("py",&py);
		t->SetBranchAddress("pz",&pz);
		t->SetBranchAddress("vapripx",&vapripx);
		t->SetBranchAddress("vapripy",&vapripy);
		t->SetBranchAddress("vapripz",&vapripz);
		t->SetBranchAddress("recodist",&recodist);
		t->SetBranchAddress("recopx",&recopx);
		t->SetBranchAddress("recopy",&recopy);
		t->SetBranchAddress("recopz",&recopz);
		t->SetBranchAddress("recophi",&recophi);
		t->SetBranchAddress("recotheta",&recotheta);
		t->SetBranchAddress("nclus",&nclus);
		t->SetBranchAddress("nlclus",&nlclus);
		t->SetBranchAddress("nrclus",&nrclus);
		t->SetBranchAddress("neclus",&neclus);
		t->SetBranchAddress("mcpx",&mcpx);
		t->SetBranchAddress("mcpy",&mcpy);
		t->SetBranchAddress("mcpz",&mcpz);
		t->SetBranchAddress("ntrk",&ntrk); 
		t->SetBranchAddress("ngrtrk",&ngrtrk); 
	};
	auto setupMCTrkTree = [&](TTree *t)
	{
		t->SetBranchAddress("zet",&zet);
		t->SetBranchAddress("aoq",&aoq);
		t->SetBranchAddress("vtxx",&vtxx);
		t->SetBranchAddress("vtxy",&vtxy);
		t->SetBranchAddress("vtxz",&vtxz);
		t->SetBranchAddress("mcpx",&mcpx);
		t->SetBranchAddress("mcpy",&mcpy);
		t->SetBranchAddress("mcpz",&mcpz);
		t->SetBranchAddress("mctheta",&mctheta);
		t->SetBranchAddress("ggclose",&ggclose);
		t->SetBranchAddress("ntrk",&ntrk); 
	};

	for(auto i: types)for(auto j: parts){
		setupTrkTree(trkTree[i][j]);
		setupMCTrkTree(mctrkTree[i][j]);
	}

	auto vcutFile = TFile::Open("./Vertex.root");
	TF1* fitVz[2]; TF2* fitVxy[2];
	Double_t vzPar[2][4]={}, vxyPar[2][6]={};
	auto beamCutFile = TFile::Open("./beamGate.root");
	TCutG* beamPICut[2];
	for(auto i: types){
		auto beamMass = i==0?132:108;
		beamCutFile->GetObject(Form("sigma30_%dSn",beamMass),beamPICut[i]);
		beamPICut[i]->SetVarX("aoq");
		beamPICut[i]->SetVarY("zet");
		vcutFile->GetObject(Form("f1Vz_%dSn",beamMass),fitVz[i]);
		fitVz[i]->GetParameters(vzPar[i]);
	}

    std::cout << "vzPar[0][0] : " << vzPar[0][0] << std::endl;
    std::cout << "vzPar[0][1] : " << vzPar[0][1] << std::endl;
    std::cout << "vzPar[0][2] : " << vzPar[0][2] << std::endl;
    std::cout << "vzPar[0][3] : " << vzPar[0][3] << std::endl;
    std::cout << "vzPar[1][0] : " << vzPar[1][0] << std::endl;
    std::cout << "vzPar[1][1] : " << vzPar[1][1] << std::endl;
    std::cout << "vzPar[1][2] : " << vzPar[1][2] << std::endl;
    std::cout << "vzPar[1][3] : " << vzPar[1][3] << std::endl;

    std::cout << vxyPar[0][0] << std::endl;
    std::cout << vxyPar[0][1] << std::endl;
    std::cout << vxyPar[0][2] << std::endl;
    std::cout << vxyPar[0][3] << std::endl;
    std::cout << vxyPar[0][4] << std::endl;
    std::cout << vxyPar[0][5] << std::endl;

    std::cout << vxyPar[1][0] << std::endl;
    std::cout << vxyPar[1][1] << std::endl;
    std::cout << vxyPar[1][2] << std::endl;
    std::cout << vxyPar[1][3] << std::endl;
    std::cout << vxyPar[1][4] << std::endl;
    std::cout << vxyPar[1][5] << std::endl;

	auto eventCut = [&](int i, double z, double aoq, double vz, double vx, double vy, bool ggclose)
	{ 
        if (/*beamPICut[i]->IsInside(aoq,z) &&*/
            true
            //&& TMath::Abs(vz-vzPar[i][1])<=3.*vzPar[i][2]
            && TMath::Abs(vz-vzPar[i][1])<=10.*vzPar[i][2]
			&& TMath::Abs(vx)<=15.
            //&& TMath::Abs(vy+205)<=20.
            //&& ggclose==0
            )
            return true;

		return /*beamPICut[i]->IsInside(aoq,z) &&*/ TMath::Abs(vz-vzPar[i][1])<=3.*vzPar[i][2] && 
			TMath::Abs(vx)<=15. && TMath::Abs(vy+205)<=20. && ggclose==0;
	};

	auto outFile = new TFile(outName,"recreate");
	TH1 *h1NECRatio[nType][nPart];
	TH2 *h2PtY[nType][nPart];
	TH2 *h2ThetaPhi[nType][nPart];
	TH2 *h2MomTheta[nType][nPart];
	TH3 *h3MomThetaPhi[nType][nPart];

	TH2 *h2ThetaCor[nType][nPart][nReco][2], *h2PhiCor[nType][nPart][nReco][2], *h2MomCor[nType][nPart][nReco][2];
	TH2 *h2MomRes[nType][nPart][nReco], *h2RapRes[nType][nPart][nReco], *h2PtRes[nType][nPart][nReco];
	TH2 *h2PxRes[nType][nPart][nReco], *h2PyRes[nType][nPart][nReco], *h2PzRes[nType][nPart][nReco];
	TH2 *h2ERes[nType][nPart][nReco];
	TH2 *h2ThetaRes[nType][nPart][nReco], *h2PhiRes[nType][nPart][nReco];
	
	TH2 *h2RapDiff[nType][nPart][nReco], *h2PtDiff[nType][nPart][nReco], *h2UtDiff[nType][nPart][nReco];

	TEfficiency *e1Mom[nType][nPart];
	//TEfficiency *e2MomTheta[nType][nPart];
	//TEfficiency *e2MomTheta_ForCut[nType][nPart];
	TEfficiency *e2PtY[nType][nPart];
	TEfficiency *e2PtYSevere[nType][nPart][nReco];
	TEfficiency *e3MomThetaPhi[nType][nPart];

	TH2D *e2MomTheta[nType][nPart];
	TH2D *e2MomTheta_ForCut[nType][nPart];
	//TH3D *e3MomThetaPhi[nType][nPart];
	TH2D h2Bin[nPart];
	for(auto j: parts)h2Bin[j] = TH2D(Form("h2Bin%d",j),"",40,-2,2,50,0,2);

	double residRange=0.8;
	double phiRange=0.5*TMath::Pi();
	for(auto i: types)
        for(auto j: parts)
        {
            TString title = partNames[j]+" "+typeNames[i];
            h1NECRatio[i][j] = new TH1D("h1NEmbedClRatio"+nameId(i,j),title+";#it{n}_{EmbedCluster}/#it{n}_{Cluster};",110,0,1.1);
            h2PtY[i][j] = new TH2D("h2PtY"+nameId(i,j),title+";#it{y}^{(0)};#it{p_{T}} (GeV/#it{c})",40,-2,2,50,0,2);
            h2ThetaPhi[i][j] = new TH2D("h2ThetaPhi"+nameId(i,j),title+";#Theta (radian);#varphi (radian)",100,0,TMath::PiOver2(),100,-phiRange,phiRange);
            h2MomTheta[i][j] = new TH2D("h2MomTheta"+nameId(i,j),title+";p_{MC};#theta_{lab}",150,effXmin[j],effXmax[j],150,0,90);
            h3MomThetaPhi[i][j] = new TH3D("h3MomThetaPhi"+nameId(i,j),title+";p_{MC};#theta_{lab}",150,effXmin[j],effXmax[j],150,0,90,300,-180,180);
            for(auto k: TSeqI(3)){
                for(auto l: TSeqI(2)){
                    h2ThetaCor[i][j][k][l] = new TH2D(Form("h2ThetaCor"+nameId(i,j,k)+"_%d",l),title+";#Theta^{MC};#Theta^{"+recoTypes[k]+"}",100,0,1.6,100,0,1.6);
                    h2PhiCor[i][j][k][l] = new TH2D(Form("h2PhiCor"+nameId(i,j,k)+"_%d",l),title+";#varphi^{MC};#varphi^{"+recoTypes[k]+"}",100,-phiRange,phiRange,100,-phiRange,phiRange);
                    h2MomCor[i][j][k][l] = new TH2D(Form("h2MomCor"+nameId(i,j,k)+"_%d",l),title+";#it{p}^{MC} (MeV/c);#it{p}^{"+recoTypes[k]+"} (MeV/c)",100,0,2400+600*j,100,0,2400+600*j);
                }
                h2MomRes[i][j][k] = new TH2D("h2MomResid"+nameId(i,j,k),title+";#it{p}^{MC} (MeV/c);Residual: (#it{p}^{MC}-#it{p}^{"+recoTypes[k]+"})/#it{p}^{MC}",100,0,2400+1000*j,100,-residRange,residRange);
                h2RapRes[i][j][k] = new TH2D("h2RapResid"+nameId(i,j,k),title+";#it{y}^{MC} ;Residual: (#it{y}^{MC}-#it{y}^{"+recoTypes[k]+"})/#it{y}^{MC}",100,-2.5,2.5,100,-residRange,residRange);
                h2ERes[i][j][k]  = new TH2D("h2EResid"+nameId(i,j,k),title+";#it{E}^{MC} ;Residual: (#it{E}^{MC}-#it{E}^{"+recoTypes[k]+"})/#it{E}^{MC}",100,0.9*(j+1),2.+j*1.2,100,-residRange,residRange);
                h2PtRes[i][j][k]  = new TH2D("h2PtResid"+nameId(i,j,k),title+";#it{p}^{MC}_{#it{T}} (|#it{y}^{(0)}| #leq 0.2) (GeV/c);Residual: (#it{p}^{MC}_{#it{T}}-#it{p}^{"+recoTypes[k]+"}_{#it{T}})/#it{p}^{MC}_{#it{T}}",100,0,1.3+0.2*j,100,-residRange,residRange);
                h2PxRes[i][j][k]  = new TH2D("h2PxResid"+nameId(i,j,k),title+";#it{p}^{MC}_{#it{x}} (MeV/c);Residual: (#it{p}^{MC}_{#it{x}}-#it{p}^{"+recoTypes[k]+"}_{#it{x}})/#it{p}^{MC}_{#it{x}}",100,0,1300+200*j,100,-residRange,residRange);
                h2PyRes[i][j][k]  = new TH2D("h2PyResid"+nameId(i,j,k),title+";#it{p}^{MC}_{#it{y}} (MeV/c);Residual: (#it{p}^{MC}_{#it{y}}-#it{p}^{"+recoTypes[k]+"}_{#it{y}})/#it{p}^{MC}_{#it{y}}",100,-600-200*j,600+200*j,100,-residRange,residRange);
                h2PzRes[i][j][k]  = new TH2D("h2PzResid"+nameId(i,j,k),title+";#it{p}^{MC}_{#it{z}} (MeV/c);Residual: (#it{p}^{MC}_{#it{z}}-#it{p}^{"+recoTypes[k]+"}_{#it{z}})/#it{p}^{MC}_{#it{z}}",100,0,2000+1200*j,100,-residRange,residRange);
                h2ThetaRes[i][j][k] = new TH2D("h2ThetaResid"+nameId(i,j,k),title+";#Theta^{MC} (radian);Residual: (#Theta^{MC}-#Theta^{"+recoTypes[k]+"})/#Theta^{MC}",100,0,1.6,100,-residRange,residRange);
                h2PhiRes[i][j][k] = new TH2D("h2PhiResid"+nameId(i,j,k),title+";#varphi^{MC} (radian);Residual: (#varphi^{MC}-#varphi^{"+recoTypes[k]+"})/#varphi^{MC}",100,-phiRange,phiRange,100,-residRange,residRange);

                h2RapDiff[i][j][k] = new TH2D("h2RapDiff"+nameId(i,j,k),title+";#it{y}^{MC} ;difference: (#it{y}^{MC}-#it{y}^{"+recoTypes[k]+"})",100,-2.5,2.5,100,-1,1);
                h2PtDiff[i][j][k]  = new TH2D("h2PtDiff"+nameId(i,j,k),title+";#it{p}^{MC}_{#it{T}} (|#it{y}^{(0)}| #leq 0.2) (GeV/c);difference: (#it{p}^{MC}_{#it{T}}-#it{p}^{"+recoTypes[k]+"}_{#it{T}})",100,0,1.3+0.2*j,100,-0.3,0.3);
                h2UtDiff[i][j][k]  = new TH2D("h2UtDiff"+nameId(i,j,k),title+";#it{u}^{MC}_{#it{T}} (|#it{y}^{(0)}| #leq 0.2) (GeV/c);difference: (#it{u}^{MC}_{#it{T}}-#it{u}^{"+recoTypes[k]+"}_{#it{T}})",100,0,1.5,100,-0.3,0.3);
            }
            e1Mom[i][j] = new TEfficiency("e1Mom"+nameId(i,j),title+";#it{p}^{MC};Efficiency",100,0,2000+600*j);
            //e2MomTheta[i][j] = new TEfficiency("e2MomTheta"+nameId(i,j),title+";#it{p}^{MC};Theta",150,effXmin[j],effXmax[j],150,0,90);
            //e2MomTheta_ForCut[i][j] = new TEfficiency("e2MomTheta_ForCut"+nameId(i,j),title+";#it{p}^{MC};Theta",150,effXmin[j],effXmax[j],150,0,90);
            e2PtY[i][j] = new TEfficiency("e2PtY"+nameId(i,j),title+";#it{y}^{(0)};#it{p_{T}} (GeV/#it{c})",40,-2,2,50,0,2);
            for(auto k: TSeqI(3))
                e2PtYSevere[i][j][k] = new TEfficiency("e2PtYSevere"+nameId(i,j,k),title+" "+recoTypes[k]+"Momentum;#it{y}^{(0)};#it{p_{T}} (GeV/#it{c})",40,-2,2,50,0,2);
            e2MomTheta[i][j] = new TH2D("e2MomTheta"+nameId(i,j),title+";#it{p}^{va};Theta_{lab}",150,effXmin[j],effXmax[j],150,0,90);
            e2MomTheta_ForCut[i][j] = new TH2D("e2MomTheta_ForCut"+nameId(i,j),title+";#it{p}^{va};Theta_{lab}",150,effXmin[j],effXmax[j],150,0,90);
            //e3MomThetaPhi[i][j] = new TH3D("e3MomThetaPhi"+nameId(i,j),title+";#it{p}^{va};Theta_{lab}",150,effXmin[j],effXmax[j],150,0,90,300,-180,180);
            e3MomThetaPhi[i][j] = new TEfficiency("e3MomThetaPhi"+nameId(i,j),title+";#it{p}^{va};Theta_{lab}",150,effXmin[j],effXmax[j],150,0,90,300,-180,180);


            std::cout<<" Analyzing "<<effFile[i][j]->GetName()<<std::endl;

            auto nEntries = mctrkTree[i][j]->GetEntries();
            for(auto entry: MakeSeq(nEntries)){
                if(entry==0) std::cout<<" MC track Tree ana";
                if(entry%100000==0) std::cout<<" ."<<std::flush;
                if(entry==nEntries-1) std::cout<<" -> finish analysis!!"<<std::endl;
                mctrkTree[i][j]->GetEntry(entry);

                if(!eventCut(i,zet,aoq,vtxz,vtxx,vtxy,ggclose)) continue;

                TVector3 labMom(mcpx,mcpy,mcpz);
                labMom.RotateY(invProjA[i]);
                double mcphi = TMath::RadToDeg()*labMom.Phi();
                //std::cout << mcphi << std::endl;
                //if(!(abs(mcphi)<=30.||abs(mcphi)>=150.)) continue;

                Double_t E = TMath::Sqrt(labMom.Mag2()+pidMass[j]*pidMass[j]);
                TLorentzVector lab4Mom(labMom,E);
                double mcy  = lab4Mom.Rapidity()/yNN[i]-1.;
                double mcpt = lab4Mom.Pt()/1000.;

                h2PtY[i][j]->Fill(mcy,mcpt);
                h2ThetaPhi[i][j]->Fill(labMom.Theta(),labMom.Phi());
                h2MomTheta[i][j]->Fill(labMom.Mag(),labMom.Theta()*TMath::RadToDeg());
                h3MomThetaPhi[i][j]->Fill(labMom.Mag(),labMom.Theta()*TMath::RadToDeg(),labMom.Phi()*TMath::RadToDeg());

            }
            mctrkTree[i][j]->Delete();


            TVector3 oldMom(0.,0.,0.);
            Bool_t   isTrackFound=kFALSE;
            nEntries = trkTree[i][j]->GetEntries();
            for(auto entry: MakeSeq(nEntries)){
                if(entry==0) std::cout<<" Reco. track Tree ana";
                if(entry%100000==0) std::cout<<" ."<<std::flush;
                if(entry==nEntries-1) std::cout<<" -> finish analysis!!"<<std::endl;
                trkTree[i][j]->GetEntry(entry);

                if(!eventCut(i,zet,aoq,vtxz,vtxx,vtxy,ggclose)) continue;


                TVector3 initMom(mcpx,mcpy,mcpz);
                initMom.RotateY(invProjA[i]);
                Double_t initE = TMath::Sqrt(initMom.Mag2()+pidMass[j]*pidMass[j]);
                TLorentzVector init4Mom(initMom,initE);

                double mcP  = initMom.Mag();
                double mcy  = init4Mom.Rapidity()/yNN[i]-1;
                double mcpt = init4Mom.Pt()/1000.;
                double mcut = init4Mom.Pt()/init4Mom.M()/uNN[i];
                double mcPx = initMom.X();
                double mcPy = initMom.Y();
                double mcPz = initMom.Z();
                double mcE  = init4Mom.Mt()*TMath::CosH(init4Mom.Rapidity()-yNN[i])/1000.;


                TVector3 recoMom, vaMom, vaPriMom;
                recoMom.SetXYZ(recopx,recopy,recopz);
                vaMom.SetXYZ(px,py,pz);
                vaPriMom.SetXYZ(vapripx,vapripy,vapripz);
                recoMom.RotateY(invProjA[i]);
                vaMom.RotateY(invProjA[i]);
                vaPriMom.RotateY(invProjA[i]);

                double recoP  = recoMom.Mag();
                double vaP    = vaMom.Mag();
                double vaPriP = vaPriMom.Mag();
                bool isRecoNotLowMom  = (j==0 && recoP>=100)  || (j==1&&recoP>=200)  || (j==2&&recoP>=300);
                bool isVANotLowMom    = (j==0 && vaP>=100)    || (j==1&&vaP>=200)    || (j==2&&vaP>=300);
                bool isVAPriNotLowMom = (j==0 && vaPriP>=100) || (j==1&&vaPriP>=200) || (j==2&&vaPriP>=300);

                Double_t recoE  = TMath::Sqrt(recoMom.Mag2()+pidMass[j]*pidMass[j]);
                Double_t vaE    = TMath::Sqrt(vaMom.Mag2()+pidMass[j]*pidMass[j]);
                Double_t vaPriE = TMath::Sqrt(vaPriMom.Mag2()+pidMass[j]*pidMass[j]);
                TLorentzVector reco4Mom(recoMom,recoE);
                TLorentzVector va4Mom(vaMom,vaE);
                TLorentzVector vapri4Mom(vaPriMom,vaPriE);

                double recoy  = reco4Mom.Rapidity()/yNN[i]-1;
                double vay    = va4Mom.Rapidity()/yNN[i]-1;
                double vaPriy = vapri4Mom.Rapidity()/yNN[i]-1;
                double recopt = reco4Mom.Pt()/1000.;
                double recout = reco4Mom.Pt()/reco4Mom.M()/uNN[i];
                double vapt    = va4Mom.Pt()/1000.;
                double vaut    = va4Mom.Pt()/va4Mom.M()/uNN[i];
                double vaPript = vapri4Mom.Pt()/1000.;
                double vaPriut = vapri4Mom.Pt()/vapri4Mom.M()/uNN[i];
                double recophi = TMath::RadToDeg()*recoMom.Phi();
                double vaphi   = TMath::RadToDeg()*vaMom.Phi();
                double vaPriphi= TMath::RadToDeg()*vaPriMom.Phi();


                int ctype = nlclus==0;
                double RecoP[nReco]={vaP,recoP,vaPriP};
                double RecoTheta[nReco]={vaMom.Theta(),recoMom.Theta(),vaPriMom.Theta()};
                double RecoPhi[nReco]={vaMom.Phi(),recoMom.Phi(),vaPriMom.Phi()};
                double RecoY[nReco]={vay,recoy,vaPriy};
                double RecoPt[nReco]={vapt,recopt,vaPript};
                double RecoUt[nReco]={vaut,recout,vaPriut};
                double RecoPx[nReco]={vaMom.X(),recoMom.X(),vaPriMom.X()};
                double RecoPy[nReco]={vaMom.Y(),recoMom.Y(),vaPriMom.Y()};
                double RecoPz[nReco]={vaMom.Z(),recoMom.Z(),vaPriMom.Z()};
                double RecoE[nReco]={va4Mom.Mt()*TMath::CosH(vay*yNN[i])/1000.,reco4Mom.Mt()*TMath::CosH(recoy*yNN[i])/1000.,vapri4Mom.Mt()*TMath::CosH(vaPriy*yNN[i])/1000.};
                bool   IsNotLow[nReco]={isVANotLowMom,isRecoNotLowMom,isVAPriNotLowMom};

                if( recodist<=20&&nclus>=15 ){
                    h1NECRatio[i][j]->Fill((double)neclus/nclus);
                    for(auto recotype: TSeqI(3)){
                        if(!(RecoPhi[recotype]>=-30&&RecoPhi[recotype]<=20))continue;
                        int isembed = neclus>=0.5*nclus;
                        h2ThetaCor[i][j][recotype][!isembed]->Fill(initMom.Theta(),RecoTheta[recotype]);
                        h2PhiCor[i][j][recotype][!isembed]->Fill(initMom.Phi(),RecoPhi[recotype]);
                        h2MomCor[i][j][recotype][!isembed]->Fill(mcP,RecoP[recotype]);
                        if(isembed){
                            h2MomRes[i][j][recotype]->Fill(mcP,(mcP-RecoP[recotype])/mcP);
                            h2RapRes[i][j][recotype]->Fill(mcy,(mcy-RecoY[recotype])/mcy);
                            if(abs(mcy)<=0.3)h2PtRes[i][j][recotype]->Fill(mcpt,(mcpt-RecoPt[recotype])/mcpt);
                            h2ERes[i][j][recotype]->Fill(mcE,(mcE-RecoE[recotype])/mcE);
                            h2PxRes[i][j][recotype]->Fill(mcPx,(mcPx-RecoPx[recotype])/mcPx);
                            h2PyRes[i][j][recotype]->Fill(mcPy,(mcPy-RecoPy[recotype])/mcPy);
                            h2PzRes[i][j][recotype]->Fill(mcPz,(mcPz-RecoPz[recotype])/mcPz);
                            h2ThetaRes[i][j][recotype]->Fill(initMom.Theta(),(initMom.Theta()-RecoTheta[recotype])/initMom.Theta());
                            h2PhiRes[i][j][recotype]->Fill(initMom.Phi(),(initMom.Phi()-RecoPhi[recotype])/initMom.Phi());

                            h2RapDiff[i][j][recotype]->Fill(mcy,mcy-RecoY[recotype]);
                            if(abs(mcy)<=0.3){
                                h2PtDiff[i][j][recotype]->Fill(mcpt,mcpt-RecoPt[recotype]);
                                h2UtDiff[i][j][recotype]->Fill(mcut,mcut-RecoUt[recotype]);
                            }

                            isTrackFound=kTRUE;
                        }

                    }
                }

                if( oldMom != initMom ){
                    oldMom = initMom;
                    if(recodist<=20&&nclus>=15) {
                        e2PtY[i][j]->            Fill(isTrackFound,mcy,mcpt);
                        e3MomThetaPhi[i][j]->    Fill(isTrackFound,vaMom.Mag(),vaMom.Theta()*TMath::RadToDeg(),vaMom.Phi()*TMath::RadToDeg());
                        e2MomTheta[i][j]->       Fill(vaMom.Mag(),vaMom.Theta()*TMath::RadToDeg());
                        e2MomTheta_ForCut[i][j]->Fill(vaMom.Mag(),vaMom.Theta()*TMath::RadToDeg());
                        e2MomTheta[i][j]->       Fill(isTrackFound,vaMom.Mag(),vaMom.Theta()*TMath::RadToDeg());
                        e2MomTheta_ForCut[i][j]->Fill(isTrackFound,vaMom.Mag(),vaMom.Theta()*TMath::RadToDeg());
                        e2MomTheta[i][j]->       Fill(isTrackFound,mcP,initMom.Theta()*TMath::RadToDeg());
                        e2MomTheta_ForCut[i][j]->Fill(isTrackFound,mcP,initMom.Theta()*TMath::RadToDeg());
                        e1Mom[i][j]->Fill(isTrackFound,mcP);
                        for(auto recotype: TSeqI(3)){
                            bool isSameBin = h2Bin[j].FindBin(mcy,mcpt)==h2Bin[j].FindBin(RecoY[recotype],RecoPt[recotype]);
                            e2PtYSevere[i][j][recotype]->Fill(isTrackFound&&isSameBin,mcy,mcpt);
                        }
                    }
                    isTrackFound=kFALSE;
                }
            }
            trkTree[i][j]->Delete();

            //---------------------------------------------------------------------------------

            e2PtY[i][j]->Write();
            e3MomThetaPhi[i][j]->Write();
            e2MomTheta[i][j]->Write();
            e2MomTheta_ForCut[i][j]->Write();
            e1Mom[i][j]->Write();
            for(auto recotype: TSeqI(3)){
                e2PtYSevere[i][j][recotype]->Write();
            }

            //---------------------------------------------------------------------------------

        }
	outFile->cd();
	outFile->Write();
}
