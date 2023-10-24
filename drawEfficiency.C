//#include "ejungwooA.h"

void drawEfficiency()
{
    gStyle -> SetOptStat(0);
    //auto file = new TFile("efficiency3.root");
    //auto file = new TFile("efficiency_normal.root");
    //auto file = new TFile("efficiency_allfound.root");
    //auto file = new TFile("efficiency_allfound.root");
    auto file = new TFile("/home/ejungwoo/efficiency_hokusai_ejungwoo/efficiency2.root");
    auto eff = (TEfficiency*) file -> Get("e3MomThetaPhi_Embed132Sn_Proton");
    auto hist = (TH3D*) eff -> GetTotalHistogram();
    auto nx = hist -> GetXaxis() -> GetNbins();
    auto x1 = hist -> GetXaxis() -> GetXmin();
    auto x2 = hist -> GetXaxis() -> GetXmax();
    auto dx = (x2-x1)/nx;
    auto ny = hist -> GetYaxis() -> GetNbins();
    auto y1 = hist -> GetYaxis() -> GetXmin();
    auto y2 = hist -> GetYaxis() -> GetXmax();
    auto dy = (y2-y1)/ny;
    auto nz = hist -> GetZaxis() -> GetNbins();
    auto z1 = hist -> GetZaxis() -> GetXmin();
    auto z2 = hist -> GetZaxis() -> GetXmax();
    auto dz = (z2-z1)/nz;
    //auto hist1 = new TH2D("hist1",";p;theta",nx,x1,x2,ny,y1,y2);
    //auto hist2 = new TH2D("hist2",";p;phi",nx,x1,x2,ny,y1,y2);
    //auto hist3 = new TH2D("hist3",";theta;phi",nx,x1,x2,ny,y1,y2);
    auto histEfficiency = new TH1D("hist",";efficiency;",201,0,1.05);

    for (auto x=x1+0.5*dx; x<x2; x+=dx)
    for (auto y=y1+0.5*dy; y<y2; y+=dy)
    for (auto z=z1+0.5*dz; z<z2; z+=dz)
    {
        auto bin = eff -> FindFixBin(x,y,z);
        auto val = eff -> GetEfficiency(bin);
        //hist1 -> Fill(x,y,val);
        //hist2 -> Fill(x,z,val);
        //hist3 -> Fill(y,z,val);
        histEfficiency -> Fill(val);
    }

    //auto cvs = ejungwoo::Canvas("cvs",100,50,3,1);
    //cvs -> cd(1); hist1 -> Draw("colz");
    //cvs -> cd(3); hist2 -> Draw("colz");
    //cvs -> cd(2); hist3 -> Draw("colz");

    auto cvs2 = new TCanvas("cvs2","cvs2",1000,800);
    histEfficiency -> Draw();
}
