void makeEfficiencyHist()
{
    // option 1) pt phi rapidity efficiency ----------------------------------------------------
    const int histType = 1;
    const bool makeH3FromH3 = true;
    int smoothingValue = 4;

    // option 2) mom theta phi efficiency ----------------------------------------------------
    //const int histType = 0;
    //const bool makeH3FromH3 = false;
    //int smoothingValue = 4;

    // write & draw options ----------------------------------------------------
    const bool writeOutputTree = false;
    const bool writeEfficiency2 = false;
    const bool drawEfficiency2 = false;
    const bool drawHist2 = false;
    int numTestVar = 20;

    // names ---------------------------------------------------------------------------------------------------

    TString nameTest; if ((drawEfficiency2 || drawHist2) && numTestVar>0) nameTest = "test.";
    TString nameFileOut = Form("efficiency.MomThetaPhi.sm%d.%s.%sroot",smoothingValue,(makeH3FromH3?"proj":"fill"),nameTest.Data());
    if (histType==1) nameFileOut = Form("efficiency.PtPhiRapidity.sm%d.%s.%sroot",smoothingValue,(makeH3FromH3?"proj":"fill"),nameTest.Data());

    //const char* pathToEmbeddingTrees = "/home/ejungwoo/data/spirit/efficiency/"; // tree_proton_embed108 ...
    const char* pathToEmbeddingTrees = "/Users/ejungwoo/data/spirit/efficiency/"; // tree_proton_embed108 ...

    gStyle -> SetOptStat(0);

    cout << "== hist-type :              " << histType << endl;
    cout << "== smoothing value :        " << smoothingValue << endl;
    cout << "== make 3d from smooth 2d : " << makeH3FromH3 << endl;
    cout << "== input path :             " << pathToEmbeddingTrees << endl;
    cout << "== output file :            " << nameFileOut << endl;

    // general ---------------------------------------------------------------------------------------------------

    if (!makeH3FromH3)
        smoothingValue = 0;

    const int numTotalPass = 2;
    const int kTotal = 0;
    const int kPassed = 1;

    const int numResolution = 2;
    const int kHigh = 0;
    const int kLow = 1;
    int arrayHL[] = {kHigh};
    //int arrayHL[] = {kHigh,kLow};

    //int arrayParticles[] = {0};
    int arrayParticles[] = {0,1,2,3,4};

    // general parameters ---------------------------------------------------------------------------------------------------

    const int iSys = 0;
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
    //const double pt2[5]      = { 1200, 1400, 1600, 1600, 1800 };
    const double pt2[5]      = { 1400, 1600, 1800, 1800, 2000 };
    const double y2          = 1.6;
    const TString pidNames[] = { "proton","deuteron","triton","he3","alpha" };

    int numP0 = 100;
    const int numTheta0 = 100;
    const int numRapidity0 = 100;
    const int numPt0 = 100;
    int numPhi0 = 100;
    if ((drawEfficiency2 || drawHist2) && numTestVar>0)
        numP0 = numTestVar;
    numTestVar = numP0;
    double testVarBinSize = p2[4]/numTestVar;

    if (histType==1) {
        numPhi0 = 100;
        if ((drawEfficiency2 || drawHist2) && numTestVar>0)
            numPhi0 = numTestVar;
        numTestVar = numPhi0;
        testVarBinSize = 360./numTestVar;
    }
    cout << "== # of test variables :    " << numTestVar << endl;

    // vertex and beam ---------------------------------------------------------------------------------------------------

    auto fileVT = new TFile("Vertex.root");
    auto fitVz = (TF1*) fileVT -> Get("f1Vz_108Sn");
    double vzPar[4] = {0};
    fitVz -> GetParameters(vzPar);

    // output ---------------------------------------------------------------------------------------------------

    auto fileOut = new TFile(nameFileOut,"recreate");

    for (auto iParticle : arrayParticles)
    {
        const char *namePID = pidNames[iParticle].Data();
        cout << "== " << namePID << endl;
        // output ---------------------------------------------------------------------------------------------------

        fileOut -> cd();
        TH3D* h3[numTotalPass][numResolution];
        TH2D* h2[numTotalPass][numResolution][numTestVar];
        for (auto iTP : {kTotal,kPassed}) {
            const char* nameTP = (iTP==kTotal?"Total":"Passed");
            for (auto iHL : arrayHL) {
                int numP = numP0;
                int numPt = numPt0;
                int numTheta = numTheta0;
                int numPhi = numPhi0;
                int numRapidity = numRapidity0;
                if (iHL==kLow) {
                    numP = numP0/3;
                    numPt = numPt0/3;
                    numTheta = numTheta0/3;
                    numPhi = numPhi0/3;
                    numRapidity = numRapidity0/3;
                }
                const char* nameHL = (iHL==kHigh?"Normal":"LowRes");
                TString name3 = Form("h3_%s_%s_%s",namePID,nameTP,nameHL);
                if (histType==0) {
                    TString title3 = Form("%s, %s, %s;#it{p}^{va};#theta_{lab};#phi_{lab}",namePID,nameTP,nameHL);
                    h3[iTP][iHL] = new TH3D(name3,title3,numP,p1[iParticle],p2[iParticle],numTheta,0,90,numPhi,-180,180);
                    cout << "== histogram :              " << name3 << ", " << title3 << ", " << numP << ", " << numTheta << ", " << numPhi << endl;
                    for (auto iVar=0; iVar<numTestVar; ++iVar) {
                        TString name2 = Form("h2_%s_%s_%s_%d",namePID,nameTP,nameHL,iVar);
                        TString title2 = Form("%s, %s, %s, p=%d-%d;#theta_{lab};#phi_{lab}",namePID,nameTP,nameHL,int(testVarBinSize*iVar),int(testVarBinSize*(iVar+1)));
                        h2[iTP][iHL][iVar] = new TH2D(name2,title2,numTheta,0,90,numPhi,-180,180);
                    }
                }
                else if (histType==1) {
                    TString title3 = Form("%s, %s, %s;#it{p}_{T};#phi_{CM};y_{CM}",namePID,nameTP,nameHL);
                    h3[iTP][iHL] = new TH3D(name3,title3,numPt,0,pt2[iParticle],numPhi,0,360,numRapidity,0,y2);
                    cout << "== histogram :              " << name3 << ", " << title3 << ", " << numPt << ", " << numPhi << ", " << numRapidity << endl;
                    for (auto iVar=0; iVar<numTestVar; ++iVar) {
                        TString name2 = Form("h2_%s_%s_%s_%d",namePID,nameTP,nameHL,iVar);
                        TString title2 = Form("%s, %s, %s, phi=%d-%d;y_{CM};p_{T}",namePID,nameTP,nameHL,int(testVarBinSize*iVar),int(testVarBinSize*(iVar+1)));
                        h2[iTP][iHL][iVar] = new TH2D(name2,title2,numRapidity,0,y2,numPt,0,pt2[iParticle]);
                    }
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

        TString nameFile = Form("%s/tree_%s_embed108.root",pathToEmbeddingTrees,namePID);
        auto fileIn = new TFile(nameFile);
        auto treeTrack = (TTree*) fileIn -> Get("trktree");

        double zet,aoq;
        double vtxx,vtxy,vtxz;
        double dist,recodist;
        double px,py,pz;
        double vapripx,vapripy,vapripz;
        double recopx,recopy,recopz;
        double mcpx,mcpy,mcpz;
        double mccmpx ,mccmpy ,mccmy ,mcphi;
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
        treeTrack -> SetBranchAddress("mccmpx",&mccmpx);
        treeTrack -> SetBranchAddress("mccmpy",&mccmpy);
        treeTrack -> SetBranchAddress("mccmy",&mccmy );
        treeTrack -> SetBranchAddress("mccmphi",&mcphi);

        // track loop ---------------------------------------------------------------------------------------------------

        TVector3 oldMom(0.,0.,0.);
        bool isTrackFound = false;
        auto numTracks = treeTrack -> GetEntries();
        for (auto iTrack=0; iTrack<numTracks; ++iTrack)
        {
            treeTrack -> GetEntry(iTrack);
            mcphi = mcphi*TMath::RadToDeg();
            mcphi = mcphi + 180;

            // vertex cut ---------------------------------------------------------------------------------------------------

            if (!(true
                        //&& beamPICut[i]->IsInside(aoq,z)
                        && TMath::Abs(vtxz-vzPar[1])<=3.*vzPar[2]
                        && TMath::Abs(vtxx)<=15.
                        && TMath::Abs(vtxy+205)<=20.
                 )) continue;

            TVector3 initMom(mcpx,mcpy,mcpz);
            initMom.RotateY(invProjA[iSys]);
            TVector3 recoMom, vaMom, vaPriVar;
            recoMom.SetXYZ(recopx,recopy,recopz);
            vaMom.SetXYZ(px,py,pz);
            vaPriVar.SetXYZ(vapripx,vapripy,vapripz);
            recoMom.RotateY(invProjA[iSys]);
            vaMom.RotateY(invProjA[iSys]);
            vaPriVar.RotateY(invProjA[iSys]);
            auto momMag = vaMom.Mag();
            auto theta_deg = vaMom.Theta()*TMath::RadToDeg();
            auto phi_deg = vaMom.Phi()*TMath::RadToDeg();
            double mcpt = sqrt(mccmpx*mccmpx+mccmpy*mccmpy);

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

                    for (auto iHL : arrayHL) {
                        if (histType==0) {
                            if (!makeH3FromH3)
                                h3[iTP][iHL] -> Fill(momMag,theta_deg,phi_deg);
                            for (auto iVar=0; iVar<numTestVar; ++iVar) {
                                if (momMag>testVarBinSize*iVar && momMag<testVarBinSize*(iVar+1))
                                    h2[iTP][iHL][iVar] -> Fill(theta_deg,phi_deg);
                            }
                        }
                        else if (histType==1) {
                            if (!makeH3FromH3)
                                h3[iTP][iHL] -> Fill(mcpt, mcphi, mccmy);
                            for (auto iVar=0; iVar<numTestVar; ++iVar) {
                                if (mcphi>testVarBinSize*iVar && mcphi<testVarBinSize*(iVar+1))
                                    h2[iTP][iHL][iVar] -> Fill(mccmy, mcpt);
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

        for (auto iHL : arrayHL)
        {
            const char* nameHL = (iHL==kHigh?"HighRes":"LowRes");

            TCanvas* cvs = nullptr;
            if (drawEfficiency2) {
                cvs = new TCanvas(Form("cvs_e2_%s_%s",namePID,nameHL),Form("cvs_e2_%s_%s",namePID,nameHL),1150,700);
                     if (numTestVar>50) cvs -> Divide(10,8,0,0);
                else if (numTestVar>40) cvs -> Divide(10,5,0,0);
                else if (numTestVar>35) cvs -> Divide(8,5,0,0);
                else if (numTestVar>20) cvs -> Divide(7,5,0,0);
                else if (numTestVar>12) cvs -> Divide(5,4,0,0);
                else if (numTestVar>8)  cvs -> Divide(4,3,0,0);
                else if (numTestVar>4)  cvs -> Divide(4,2,0,0);
                else if (numTestVar>1)  cvs -> Divide(2,2,0,0);
            }

            for (auto iVar=0; iVar<numTestVar; ++iVar) {
                const char* nameHL = (iHL==kHigh?"HighRes":"LowRes");
                TString name2 = Form("e2_%s_%s_%d",namePID,nameHL,iVar);
                if (smoothingValue>0) {
                    for (auto iSmooth=0; iSmooth<smoothingValue; ++iSmooth) {
                        h2[kTotal][iHL][iVar] -> Smooth();
                        h2[kPassed][iHL][iVar] -> Smooth();
                    }
                }

                auto max = h2[kTotal][iHL][iVar] -> GetMaximum();
                h2[kPassed][iHL][iVar] -> SetMaximum(max);
                auto e2 = new TEfficiency(*h2[kPassed][iHL][iVar], *h2[kTotal][iHL][iVar]);
                e2 -> SetNameTitle(name2,h2[kPassed][iHL][iVar]->GetTitle());
                if (drawEfficiency2) {
                    cvs -> cd(iVar+1);
                    e2 -> Draw("colz");
                }

                if (writeEfficiency2)
                    e2 -> Write();
            }

            if (makeH3FromH3)
            {
                int numP = numP0;
                int numPt = numPt0;
                int numTheta = numTheta0;
                int numPhi = numPhi0;
                int numRapidity = numRapidity0;
                if (iHL==kLow) {
                    numP = numP0/3;
                    numPt = numPt0/3;
                    numTheta = numTheta0/3;
                    numPhi = numPhi0/3;
                    numRapidity = numRapidity0/3;
                }
                for (auto iTP : {kTotal,kPassed}) {
                    if (histType==0) {
                        for (auto iP=0; iP<numP; ++iP) {
                            for (auto iTheta=0; iTheta<numTheta; ++iTheta) {
                                for (auto iPhi=0; iPhi<numPhi; ++iPhi) {
                                    auto content = h2[iTP][iHL][iP] -> GetBinContent(iTheta+1,iPhi+1);
                                    h3[iTP][iHL] -> SetBinContent(iP+1,iTheta+1,iPhi+1,content);
                                }
                            }
                        }
                    }
                    else if (histType==1) {
                        for (auto iPt=0; iPt<numPt; ++iPt) {
                            for (auto iPhi=0; iPhi<numPhi; ++iPhi) {
                                for (auto iRapidity=0; iRapidity<numRapidity; ++iRapidity) {
                                    auto content = h2[iTP][iHL][iPhi] -> GetBinContent(iRapidity+1,iPt+1);
                                    h3[iTP][iHL] -> SetBinContent(iPt+1,iPhi+1,iRapidity+1,content);
                                }
                            }
                        }
                    }
                }
            }

            TString name3 = Form("e3_%s_%s",namePID,nameHL);
            auto e3 = new TEfficiency(*h3[kPassed][iHL], *h3[kTotal][iHL]);
            e3 -> SetNameTitle(name3,h3[kPassed][iHL]->GetTitle());
            e3 -> Write();
        }

        if (drawHist2)
        {
            for (auto iTP : {kTotal,kPassed}) {
                const char* nameTP = (iTP==kTotal?"Total":"Passed");
                auto cvs = new TCanvas(Form("cvs_%s_%s",nameTP,namePID),Form("cvs_%s_%s",nameTP,namePID),1150,700);
                     if (numTestVar>50) cvs -> Divide(10,8,0,0);
                else if (numTestVar>40) cvs -> Divide(10,5,0,0);
                else if (numTestVar>35) cvs -> Divide(8,5,0,0);
                else if (numTestVar>20) cvs -> Divide(7,5,0,0);
                else if (numTestVar>12) cvs -> Divide(5,4,0,0);
                else if (numTestVar>8)  cvs -> Divide(4,3,0,0);
                else if (numTestVar>4)  cvs -> Divide(4,2,0,0);
                else if (numTestVar>1)  cvs -> Divide(2,2,0,0);
                if (drawHist2) {
                    for (auto iVar=0; iVar<numTestVar; ++iVar) {
                        cvs -> cd(iVar+1);
                        h2[iTP][kHigh][iVar] -> Draw("colz");
                    }
                }
            }
        }


        if (writeOutputTree)
            treeOut -> Write();
    }
}
