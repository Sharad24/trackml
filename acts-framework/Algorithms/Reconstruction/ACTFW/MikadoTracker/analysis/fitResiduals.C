/*
  root -l analysis/fitResiduals.C'("bV0&&pt>=0.5")'
  root -l analysis/fitResiduals.C'("prim&&bV0&&phiV0<0.025")'
 */



int fitResiduals( const char *inputCuts=0 )
{
  TFile *f = new TFile("statFit.root","READ");
  
  TNtuple *fit = (TNtuple*) f->FindObjectAny("fit");  

  TString cuts0 = "";
  if( inputCuts ) cuts0 = cuts0 + inputCuts + "&&";

  TString cut;    

  
  TCanvas *cPhiT = new TCanvas("cPhiT","cuts Phi T",1600,1200);
  cPhiT->Divide(4,3);  
  int ipad=1;
  for( int irad=0; irad<3; irad++ ){
    TCanvas *c = cPhiT;
    cut = cuts0 + "w>0&&vol==" + irad*3;    
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("dPL",cut+"&&fabs(dPL)<0.1");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("dTL",cut+"&&fabs(dTL)<3.");

    cut = cuts0 + "w>0&&(vol==" + (irad*3+1) + "||vol=="+ (irad*3+2)+")";
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("dPL",cut+"&&fabs(dPL)<0.1");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("dTL",cut+"&&fabs(dTL)<3.");    
  }

  TCanvas *cUV = new TCanvas("cUV","cuts V UZ",1600,1200);
  cUV->Divide(4,3);
  ipad=1;
  for( int irad=0; irad<3; irad++ ){
    TCanvas *c = cUV;    
    cut = cuts0 + "w>0&&vol==" + irad*3;    
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("dv",cut+"&&fabs(normV)<10.");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("duz",cut+"&&fabs(normUZ)<10.");

    cut = cuts0 + "w>0&&(vol==" + (irad*3+1) + "||vol=="+ (irad*3+2)+")";
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("dv",cut+"&&fabs(normV)<10.");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("duz",cut+"&&fabs(normUZ)<10.");
   }
 
  TCanvas *cPick = new TCanvas("cPick","cuts PickV PickUZ",1600,1200);
  cPick->Divide(4,3);
  ipad=1;
  for( int irad=0; irad<3; irad++ ){
    TCanvas *c = cPick;
    cut = cuts0 + "w>0&&vol==" + irad*3;    
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("pickDv",cut+"&&fabs(picknormV)<10.&&dupl>0");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("pickDuz",cut+"&&fabs(picknormUZ)<10.&&dupl>0");
    
    cut = cuts0 + "w>0&&(vol==" + (irad*3+1) + "||vol=="+ (irad*3+2)+")";
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("pickDv",cut+"&&fabs(picknormV)<10.&&dupl>0");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("pickDuz",cut+"&&fabs(picknormUZ)<10.&&dupl>0");
  }

  return 0;
}
