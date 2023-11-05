void makeEfficiencyHist()
{
    const int iSys = 0;

    // general ---------------------------------------------------------------------------------------------------

    const char* pathToEmbeddingTrees = "/home/ejungwoo/data/spirit/efficiency/"; // tree_proton_embed108 ...
    bool writeOutputTree = false;
    bool writeEfficiency2 = true;
    bool drawEfficiency2 = true;

    // general parameters ---------------------------------------------------------------------------------------------------

    const double invProjA[]  = { 44.1/1000., 55.2/1000.,45./1000. };
    const double beta[]      = {-0.364,-0.349,-0.355 };
    const double ycoll[]     = { 0.382, 0.365, 0.371 };
    const double ucoll[]     = { 0.392, 0.373, 0.371 };
    const double momCut[]    = { 2000., 4000., 4200. };
    const double yAA[]       = { 0.3822, 0.3647, 0.3538, 0.3902 };
    const double yNN[]       = { 0.3696, 0.3697, 0.3705, 0.3706 };
    const double yBeam[]     = { 0.7421, 0.7423, 0.7439, 0.7441 };
    const double uNN[]       = { 0.3709, 0.3709, 0.3718, 0.3719 };
    const double pidMass[]   = { 938.2720813, 1875.612762, 2808.921112, 2808.39132, 3727.379378};
    const double p1[5]       = { 100, 200, 100, 100, 100 };
    const double p2[5]       = { 1500, 2200, 3200, 3200, 4200 };
    const TString pidNames[] = { "proton","deuteron","triton","he3","alpha" };

    int numTP = 2;
    int numResolution = 2;
    int numTestMom = 8;
    double testMomBinSize = 400;

    const int kTotal = 0;
    const int kPassed = 1;
    const int kHigh = 0;
    const int kLow = 1;
    const int numP0 = 150;
    const int numTheta0 = 150;
    const int numPhi0 = 300;

    // vertex and beam ---------------------------------------------------------------------------------------------------

    auto fileVT = new TFile("Vertex.root");
    auto fitVz = (TF1*) fileVT -> Get("f1Vz_108Sn");
    double vzPar[4] = {0};
    fitVz -> GetParameters(vzPar);

    // output ---------------------------------------------------------------------------------------------------

    auto fileOut = new TFile("efficiency.108Sn.root","recreate");

    //for (auto iParticle : {0,1,2,3,4})
    for (auto iParticle : {0})
    {
        const char *namePID = pidNames[iParticle].Data();
        cout << "== " << namePID << endl;
        // output ---------------------------------------------------------------------------------------------------

        fileOut -> cd();
        TH3D* h3[numTP][numResolution];
        TH2D* h2[numTP][numResolution][numTestMom];
        for (auto iTP : {kTotal,kPassed}) {
            const char* nameTP = (iTP==kTotal?"Total":"Passed");
            for (auto iHL : {kHigh,kLow}) {
                int numP = numP0;
                int numTheta = numTheta0;
                int numPhi = numPhi0;
                if (iHL==kLow) {
                    numP = numP0/3;
                    numTheta = numTheta0/3;
                    numPhi = numPhi0/3;
                }
                const char* nameHL = (iHL==kHigh?"Normal":"LowRes");
                TString name3 = Form("h3_%s_%s_%s",namePID,nameTP,nameHL);
                TString title3 = Form("%s, %s, %s;#it{p}^{va};#theta_{lab};#phi_{lab}",namePID,nameTP,nameHL);
                h3[iTP][iHL] = new TH3D(name3,title3,numP,p1[iParticle],p2[iParticle],numTheta,0,90,numPhi,-180,180);
                for (auto iMom=0; iMom<numTestMom; ++iMom) {
                    TString name2 = Form("h2_%s_%s_%s_%d",namePID,nameTP,nameHL,iMom);
                    TString title2 = Form("%s, %s, %s, p=%d-%d;#theta_{lab};#phi_{lab}",namePID,nameTP,nameHL,int(testMomBinSize*iMom),int(testMomBinSize*(iMom+1)));
                    h2[iTP][iHL][iMom] = new TH2D(name2,title2,numTheta,0,90,numPhi,-180,180);
                }
            }
        }

        int oDist, oNCluster, oNClusterE;
        bool oCutCluster, oIsTrackFound;
        double oMomMag, oTheta, oPhi;
        TTree* treeOut = nullptr;
        if (writeOutputTree) {
            treeOut = new TTree("track","");
            treeOut -> Branch("dist",&oDist);
            treeOut -> Branch("nCluster",&oNCluster);
            treeOut -> Branch("nClusterE",&oNClusterE);
            treeOut -> Branch("cutCluster",&oCutCluster);
            treeOut -> Branch("isTrackFound",&oIsTrackFound);
            treeOut -> Branch("momMag",&oMomMag);
            treeOut -> Branch("theta",&oTheta);
            treeOut -> Branch("phi",&oPhi);
        }

        // input file ---------------------------------------------------------------------------------------------------

        const char* nameFile = Form("%s/tree_%s_embed108.root",pathToEmbeddingTrees,namePID);
        auto fileIn = new TFile(nameFile);
        auto treeTrack = (TTree*) fileIn -> Get("trktree");

        double zet,aoq;
        double vtxx,vtxy,vtxz;
        double dist,recodist;
        double px,py,pz;
        double vapripx,vapripy,vapripz;
        double recopx,recopy,recopz;
        double mcpx,mcpy,mcpz;
        int    nclus, nlclus, nrclus, neclus;

        treeTrack -> SetBranchAddress("zet",&zet);
        treeTrack -> SetBranchAddress("aoq",&aoq);
        treeTrack -> SetBranchAddress("vtxx",&vtxx);
        treeTrack -> SetBranchAddress("vtxy",&vtxy);
        treeTrack -> SetBranchAddress("vtxz",&vtxz);
        treeTrack -> SetBranchAddress("dist",&dist);
        treeTrack -> SetBranchAddress("px",&px);
        treeTrack -> SetBranchAddress("py",&py);
        treeTrack -> SetBranchAddress("pz",&pz);
        treeTrack -> SetBranchAddress("vapripx",&vapripx);
        treeTrack -> SetBranchAddress("vapripy",&vapripy);
        treeTrack -> SetBranchAddress("vapripz",&vapripz);
        treeTrack -> SetBranchAddress("recodist",&recodist);
        treeTrack -> SetBranchAddress("recopx",&recopx);
        treeTrack -> SetBranchAddress("recopy",&recopy);
        treeTrack -> SetBranchAddress("recopz",&recopz);
        treeTrack -> SetBranchAddress("nclus",&nclus);
        treeTrack -> SetBranchAddress("nlclus",&nlclus);
        treeTrack -> SetBranchAddress("nrclus",&nrclus);
        treeTrack -> SetBranchAddress("neclus",&neclus);
        treeTrack -> SetBranchAddress("mcpx",&mcpx);
        treeTrack -> SetBranchAddress("mcpy",&mcpy);
        treeTrack -> SetBranchAddress("mcpz",&mcpz);

        // track loop ---------------------------------------------------------------------------------------------------

        TVector3 oldMom(0.,0.,0.);
        bool isTrackFound = false;
        auto numTracks = treeTrack -> GetEntries();
        for (auto iTrack=0; iTrack<numTracks; ++iTrack)
        {
            treeTrack -> GetEntry(iTrack);

            // vertex cut ---------------------------------------------------------------------------------------------------

            if (!(true
                        //&& beamPICut[i]->IsInside(aoq,z)
                        && TMath::Abs(vtxz-vzPar[1])<=3.*vzPar[2]
                        && TMath::Abs(vtxx)<=15.
                        && TMath::Abs(vtxy+205)<=20.
                 )) continue;

            TVector3 initMom(mcpx,mcpy,mcpz);
            initMom.RotateY(invProjA[iSys]);
            TVector3 recoMom, vaMom, vaPriMom;
            recoMom.SetXYZ(recopx,recopy,recopz);
            vaMom.SetXYZ(px,py,pz);
            vaPriMom.SetXYZ(vapripx,vapripy,vapripz);
            recoMom.RotateY(invProjA[iSys]);
            vaMom.RotateY(invProjA[iSys]);
            vaPriMom.RotateY(invProjA[iSys]);
            auto momMag = vaMom.Mag();
            auto theta_deg = vaMom.Theta()*TMath::RadToDeg();
            auto phi_deg = vaMom.Phi()*TMath::RadToDeg();

            // cluster cut ---------------------------------------------------------------------------------------------------

            isTrackFound = false;
            if (recodist<=20 && nclus>=15 && neclus>=0.5*nclus)
                isTrackFound = true;

            // fill ---------------------------------------------------------------------------------------------------

            if (oldMom!=initMom) {
                oldMom = initMom;

                for (auto iTP : {kTotal,kPassed}) {
                    if (iTP==kPassed && !isTrackFound)
                        continue;

                    for (auto iHL : {kHigh,kLow}) {
                        h3[iTP][iHL] -> Fill(momMag,theta_deg,phi_deg);
                        for (auto iMom=0; iMom<numTestMom; ++iMom) {
                            for (auto iMom=0; iMom<numTestMom; ++iMom) {
                                if (momMag>testMomBinSize*iMom && momMag<testMomBinSize*(iMom+1))
                                    h2[iTP][iHL][iMom] -> Fill(theta_deg,phi_deg);
                            }
                        }
                    }
                }
            }

            if (writeOutputTree)
            {
                oDist = recodist;
                oNCluster = neclus;
                oNClusterE = nclus;
                oCutCluster = (neclus>=0.5*nclus);
                oIsTrackFound = isTrackFound;
                oMomMag = momMag;
                oTheta = theta_deg;
                oPhi =  phi_deg;
                treeOut -> Fill();
            }
        }
        fileIn -> Close();

        // write and draw ---------------------------------------------------------------------------------------------------

        fileOut -> cd();

        for (auto iHL : {kHigh,kLow}) {
            const char* nameHL = (iHL==kHigh?"HighRes":"LowRes");
            TString name3 = Form("e3_%s_%s",namePID,nameHL);
            auto e3 = new TEfficiency(*h3[kPassed][iHL], *h3[kTotal][iHL]);
            e3 -> SetNameTitle(name3,h3[kPassed][iHL]->GetTitle());
            e3 -> Write();

            TCanvas* cvs = nullptr;
            if (drawEfficiency2) {
                cvs = new TCanvas(Form("cvs_e2_%s_%s",namePID,nameHL),"",3000,2000);
                cvs -> Divide(4,2);
            }

            if (writeEfficiency2) {
                for (auto iMom=0; iMom<numTestMom; ++iMom) {
                    const char* nameHL = (iHL==kHigh?"HighRes":"LowRes");
                    TString name2 = Form("e2_%s_%s_%d",namePID,nameHL,iMom);
                    auto e2 = new TEfficiency(*h2[kPassed][iHL][iMom], *h2[kTotal][iHL][iMom]);
                    e2 -> SetNameTitle(name2,h2[kPassed][iHL][iMom]->GetTitle());
                    e2 -> Write();
                    if (drawEfficiency2) {
                        cvs -> cd(iMom+1);
                        e2 -> Draw("colz");
                    }
                }
            }
        }

        if (writeOutputTree)
            treeOut -> Write();
    }

    fileOut -> ls();
}
