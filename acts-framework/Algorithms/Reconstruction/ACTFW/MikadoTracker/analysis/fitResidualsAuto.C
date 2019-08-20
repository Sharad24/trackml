/*
  root -l analysis/fitResidualsAuto.C'("bV0&&pt>=0.5")'
 */

int fitResidualsAuto( const char *inputCuts=0 )
{
  TFile *f = new TFile("statFit.root","READ");
  
  TNtuple *fit = (TNtuple*) f->FindObjectAny("fit");  
  
  TCanvas *crad = new TCanvas("crad","c radial",1600,1200);
  crad->Divide(6,3);
  TCanvas *cvert = new TCanvas("cvert","c vert",1600,1200);
  cvert->Divide(6,3);

  TString cuts0 = "";
  if( inputCuts ) cuts0 = cuts0 + inputCuts + "&&";

  int ipad=1;
  for( int irad=0; irad<3; irad++ ){
    TCanvas *c = crad;
    TString cut;    
    
    cut = cuts0 + "w>0&&vol==" + irad*3;    
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("dPL",cut+"&&fabs(dPL)<0.1");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("dTL",cut+"&&fabs(dTL)<3.");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("normUZ",cut+"&&fabs(normUZ)<10.");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("normV",cut+"&&fabs(normV)<10.");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("picknormUZ",cut+"&&fabs(picknormUZ)<10.&&dupl>0");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("picknormV",cut+"&&fabs(picknormV)<10.&&dupl>0");
   }
  
  ipad=1;
  for( int irad=0; irad<3; irad++ ){
    TCanvas *c = cvert;
    TString cut;
    cut = cuts0 + "w>0&&(vol==" + (irad*3+1) + "||vol=="+ (irad*3+2)+")";
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("dPL",cut+"&&fabs(dPL)<0.1");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("dTL",cut+"&&fabs(dTL)<3.");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("normUZ",cut+"&&fabs(normUZ)<10.");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("normV",cut+"&&fabs(normV)<10.");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("picknormUZ",cut+"&&fabs(picknormUZ)<10.&&dupl>0");
    c->cd(ipad++);gPad->SetLogy();
    fit->Draw("picknormV",cut+"&&fabs(picknormV)<10.&&dupl>0");
  }

  return 0;
}
