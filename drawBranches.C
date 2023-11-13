void drawBranches()
{
    // input file ---------------------------------------------------------------------------------------------------

    const char* pidNames[] = { "proton","deuteron","triton","he3","alpha" };

    for (auto iParticle : {0,1,2,3,4})
    {
        auto file = new TFile(Form("/Users/ejungwoo/data/spirit/efficiency/tree_%s_embed108.root",pidNames[iParticle]));
        auto tree = (TTree*) file -> Get("trktree");

        auto cvs = new TCanvas(Form("cvs%s",pidNames[iParticle]),pidNames[iParticle],1200,700);
        cvs -> Divide(3,2);

        if (1) {
            cvs -> cd(1); tree -> Draw("mccmp");
            cvs -> cd(2); tree -> Draw("mccmtheta");
            cvs -> cd(3); tree -> Draw("mccmphi");
            cvs -> cd(4); tree -> Draw("sqrt(mccmpx*mccmpx+mccmpy*mccmpy)");
            cvs -> cd(5); tree -> Draw("mccmy");
            cvs -> cd(6); tree -> Draw("mccmphi");
        }
        else {
            cvs -> cd(1); tree -> Draw("mcp");
            cvs -> cd(2); tree -> Draw("mctheta");
            cvs -> cd(3); tree -> Draw("mcphi");
            cvs -> cd(4); tree -> Draw("sqrt(mcpx*mcpx+mcpy*mcpy)");
            cvs -> cd(5); tree -> Draw("mcy");
            cvs -> cd(6); tree -> Draw("mcphi");
        }
    }
}
