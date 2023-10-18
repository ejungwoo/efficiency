#include "ejungwooA.h"

void drawEfficiency()
{
    auto file = new TFile("efficiency2.root");
    auto eff = (TEfficiency*) file -> Get("e3MomThetaPhi_Embed132Sn_Proton");
    //auto eff = (TEfficiency*) file -> Get("e3MomThetaPhi_Embed132Sn_Alpha");
    auto hist = (TH3D*) eff -> GetTotalHistogram();
    auto bx = ejungwoo::Binning(hist -> GetXaxis());
    auto by = ejungwoo::Binning(hist -> GetYaxis());
    auto bz = ejungwoo::Binning(hist -> GetZaxis());
    auto histxy = new TH2D("histxy",";x;y",bx.nx(),bx.x1(),bx.x2(),by.nx(),by.x1(),by.x2());
    auto histyz = new TH2D("histyz",";y;z",by.nx(),by.x1(),by.x2(),bz.nx(),bz.x1(),bz.x2());
    auto histzx = new TH2D("histzx",";z;x",bz.nx(),bz.x1(),bz.x2(),bx.nx(),bx.x1(),bx.x2());

    for (auto x=bx.x1(); x<bx.nx(); x+=bx.dx()) {
        for (auto y=by.x1(); y<by.nx(); y+=by.dx()) {
            for (auto z=bz.x1(); z<bz.nx(); z+=bz.dx()) {
                auto bin = eff -> FindFixBin(x,y,z);
                auto val = eff -> GetEfficiency(bin);
                if (val>0) {
                    histxy -> Fill(x,y,val);
                    histyz -> Fill(y,z,val);
                    histzx -> Fill(z,x,val);
                    cout << x << " " << y << " " << z << " " << val << endl;
                }
            }
        }
    }

    auto cvs = ejungwoo::Canvas("cvs",100,50,3,1);
    cvs -> cd(1); histxy -> Draw("colz");
    cvs -> cd(2); histyz -> Draw("colz");
    cvs -> cd(3); histzx -> Draw("colz");
}
