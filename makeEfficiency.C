void makeEfficiency()
{
    const int iSys = 0;

    // general flags  ---------------------------------------------------------------------------------------------------

    bool writeOutputTree = false;
    bool writeEfficiency2 = true;
    bool drawEfficiency2 = true;

    // general parameters ---------------------------------------------------------------------------------------------------

    const double invProjA[] = {44.1/1000., 55.2/1000.,45./1000. };
    const double beta[]     = {-0.364, -0.349, -0.355 };
    const double ycoll[]    = { 0.382, 0.365, 0.371 };
    const double ucoll[]    = { 0.392, 0.373, 0.371 };
    const double momCut[]   = { 2000., 4000., 4200. };
    const double yAA[]      = {0.3822, 0.3647, 0.3538, 0.3902};
    const double yNN[]      = {0.3696, 0.3697, 0.3705, 0.3706};
    const double yBeam[]    = {0.7421, 0.7423, 0.7439, 0.7441};
    const double uNN[]      = {0.3709, 0.3709, 0.3718, 0.3719};
    const double amu        = 931.478;
    const double pidMass[]  = { 938.2720813, 1875.612762, 2808.921112, 2808.39132, 3727.379378};
    const TString fPIDName[] = {"proton","deuteron","triton","he3","alpha"};
    const double p1[5]      = { 100, 200, 100, 100, 100 };
    const double p2[5]      = { 1500, 2200, 3200, 3200, 4200 };

    // vertex and beam ---------------------------------------------------------------------------------------------------

    auto fileVT = new TFile("Vertex.root");
    auto fitVz = (TF1*) fileVT -> Get("f1Vz_108Sn");
    double vzPar[4] = {0};
    fitVz -> GetParameters(vzPar);

    // output ---------------------------------------------------------------------------------------------------

    auto fileOut = new TFile("efficiency.108Sn.root","recreate");

    for (auto iParticle : {0,1,2,3,4})
    {
        // output ---------------------------------------------------------------------------------------------------

        fileOut -> cd();
        const char *pidName = fPIDName[iParticle].Data();
        TEfficiency *e3MomThetaPhi = new TEfficiency(Form("e3MomThetaPhi_%s",pidName), ";#it{p}^{va};#theta_{lab};#phi_{lab}",
                150,p1[iParticle],p2[iParticle],150,0,90,300,-180,180);
        TEfficiency *e2ThetaPhi[4] = {0};
        if (writeEfficiency2) {
            e2ThetaPhi[0] = new TEfficiency(Form("e2_%s_0",pidName),    "p=0-500;#theta_{lab};#phi_{lab}",150,0,90,300,-180,180);
            e2ThetaPhi[1] = new TEfficiency(Form("e2_%s_1",pidName), "p=500-1000;#theta_{lab};#phi_{lab}",150,0,90,300,-180,180);
            e2ThetaPhi[2] = new TEfficiency(Form("e2_%s_2",pidName),"p=1000-1500;#theta_{lab};#phi_{lab}",150,0,90,300,-180,180);
            e2ThetaPhi[3] = new TEfficiency(Form("e2_%s_3",pidName),"p=1500-2000;#theta_{lab};#phi_{lab}",150,0,90,300,-180,180);
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

        auto fileIn = new TFile(Form("/Users/ejungwoo/data/spirit/efficiency/tree_%s_embed108.root",pidName));
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
                e3MomThetaPhi -> Fill(isTrackFound,momMag,theta_deg,phi_deg);
                if (writeEfficiency2) {
                    if (momMag>   0 && momMag< 500) e2ThetaPhi[0] -> Fill(isTrackFound,theta_deg,phi_deg);
                    if (momMag> 500 && momMag<1000) e2ThetaPhi[1] -> Fill(isTrackFound,theta_deg,phi_deg);
                    if (momMag>1000 && momMag<1500) e2ThetaPhi[2] -> Fill(isTrackFound,theta_deg,phi_deg);
                    if (momMag>1500)                e2ThetaPhi[3] -> Fill(isTrackFound,theta_deg,phi_deg);
                }
            }

            oDist = recodist;
            oNCluster = neclus;
            oNClusterE = nclus;
            oCutCluster = (neclus>=0.5*nclus);
            oIsTrackFound = isTrackFound;
            oMomMag = momMag;
            oTheta = theta_deg;
            oPhi =  phi_deg;
            if (writeOutputTree)
                treeOut -> Fill();
        }
        fileIn -> Close();

        // write and draw ---------------------------------------------------------------------------------------------------

        fileOut -> cd();
        e3MomThetaPhi -> Write();
        if (writeEfficiency2) {
            e2ThetaPhi[0] -> Write();
            e2ThetaPhi[1] -> Write();
            e2ThetaPhi[2] -> Write();
            e2ThetaPhi[3] -> Write();
        }
        if (writeOutputTree)
            treeOut -> Write();

        if (writeEfficiency2 && drawEfficiency2) {
            auto cvs = new TCanvas(Form("cvs_e2_%s",pidName),"",1000,720);
            cvs -> Divide(2,2);
            cvs -> cd(1); e2ThetaPhi[0] -> Draw("colz");
            cvs -> cd(2); e2ThetaPhi[1] -> Draw("colz");
            cvs -> cd(3); e2ThetaPhi[2] -> Draw("colz");
            cvs -> cd(4); e2ThetaPhi[3] -> Draw("colz");
            cvs -> Write();
        }
    }

    fileOut -> ls();
}
