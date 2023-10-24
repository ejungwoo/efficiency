#include "ejungwooA.h"

void drawEfficiency()
{
    gStyle -> SetOptStat(0);
    //auto file = new TFile("efficiency3.root");
    //auto file = new TFile("efficiency_normal.root");
    auto file = new TFile("efficiency_allfound.root");
    auto eff = (TEfficiency*) file -> Get("e3MomThetaPhi_Embed132Sn_Proton");
    auto hist = (TH3D*) eff -> GetTotalHistogram();
    auto bx = ejungwoo::Binning(hist->GetXaxis());
    auto by = ejungwoo::Binning(hist->GetYaxis());
    auto bz = ejungwoo::Binning(hist->GetZaxis());
    auto hist1 = new TH2D("hist1",";p;theta",bx.nx(),bx.x1(),bx.x2(),by.nx(),by.x1(),by.x2());
    auto hist2 = new TH2D("hist2",";p;phi",bx.nx(),bx.x1(),bx.x2(),bz.nx(),bz.x1(),bz.x2());
    auto hist3 = new TH2D("hist3",";theta;phi",by.nx(),by.x1(),by.x2(),bz.nx(),bz.x1(),bz.x2());
    auto histEfficiency = new TH1D("hist",";efficiency;",201,0,1.05);

    for (auto x=bx.x1()+0.5*bx.dx(); x<bx.x2(); x+=bx.dx())
    for (auto y=by.x1()+0.5*by.dx(); y<by.x2(); y+=by.dx())
    for (auto z=bz.x1()+0.5*bz.dx(); z<bz.x2(); z+=bz.dx())
    {
        auto bin = eff -> FindFixBin(x,y,z);
        auto val = eff -> GetEfficiency(bin);
        hist1 -> Fill(x,y,val);
        hist2 -> Fill(x,z,val);
        hist3 -> Fill(y,z,val);
        histEfficiency -> Fill(val);
    }

    auto cvs = ejungwoo::Canvas("cvs",100,50,3,1);
    cvs -> cd(1); hist1 -> Draw("colz");
    cvs -> cd(3); hist2 -> Draw("colz");
    cvs -> cd(2); hist3 -> Draw("colz");

    auto cvs2 = ejungwoo::Canvas("cvs2",2,2);
    histEfficiency -> Draw();
}
