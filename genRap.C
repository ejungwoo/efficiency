void genRap(Int_t system=108, Int_t pid=1)
{
    if (system!=108 && system!=132) {
        cout << "system is " << system << "! exit" << endl;
        return;
    }
    if (pid<0 || pid>4) {
        cout << "pid is " << pid << "! exit" << endl;
        return;
    }
    Int_t nfile=1000000;
    Int_t neve=10000;
    Double_t maxp = 1; 
    Double_t proton_mass = 0.94; // GeV/cc
    Double_t beamangle = -44.06/1000;// rad, Sn-132
    if (system==108)
        beamangle = -55.19/1000;// rad, Sn-108

    Double_t pt,px,py,pz,e,rap,mass;
    Double_t A,B;
    Int_t m, q;
    if (pid==0) { m = 1; q = 1; }
    if (pid==1) { m = 2; q = 1; }
    if (pid==2) { m = 3; q = 1; }
    if (pid==3) { m = 3; q = 2; }
    if (pid==4) { m = 4; q = 2; }
    mass = m*proton_mass;
    TRandom3 *ran = new TRandom3;
    TFile *fout;
    TNtuple *ntp;

    for(int i=0;i<nfile;i++)
    {
        std::cout << system << " " << pid << " " << i << std::endl;
        gSystem -> mkdir(Form("genRap%d/%d/%02d/",system,pid,i/10000),kTRUE);
        fout = new TFile(Form("genRap%d/%d/%02d/mc%06d_%dA.root",system,pid,i/10000,i,m),"RECREATE");
        ntp = new TNtuple("mc","eta","px:py:pz:y");
        Int_t ceve = 0;

        while(ceve<neve)
        {
            pt = ran->Uniform(0,maxp);
            ran->Circle(px,py,pt);
            rap = ran->Uniform(-0.3,1.0);
            A = TMath::Exp(2*rap);
            B = (A-1)/(A+1);
            e = sqrt((px*px+py*py+mass*mass)/(1-B*B));
            pz = B*e;

            TLorentzVector vec;
            vec.SetPxPyPzE(px,py,pz,e);
            Double_t phibefore = vec.Phi();
            vec.RotateY(beamangle);
            Double_t phiafter = vec.Phi();

            //if( abs(phibefore) <35./180.*TMath::Pi() || abs(phiafter ) <35./180.*TMath::Pi() )
            if( abs(phibefore)  < 1 || abs(phiafter )  < 1 )
            {
                if(vec.Pz()>0)
                {
                    ntp->Fill((Float_t)(vec.Px()),
                            (Float_t)(vec.Py()),
                            (Float_t)(vec.Pz()),
                            (Float_t)(vec.Rapidity()));
                    ceve ++;
                }
            }

        }
        fout->Write();
        fout->Close();
        fout = NULL;
        ntp = NULL;
    }
}
