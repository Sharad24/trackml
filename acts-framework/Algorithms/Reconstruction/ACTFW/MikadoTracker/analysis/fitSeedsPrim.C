/*
  root -l analysis/fitSeedsPrim.C'("bV0&&pt>=0.5")'
  root -l analysis/fitSeedsPrim.C'("(bV1||bV2)&&pt>=0.5")'

  root -l analysis/fitSeedsPrim.C'("ptV0>=0.5")'
*/


int fitSeedsPrim( const char *inputCuts=0 )
{
  TFile *f = new TFile("statSeeds.root","READ");
  
  TNtuple *nt = (TNtuple*) f->FindObjectAny("seeds");  

  TString cuts = "";
  if( inputCuts ) cuts = cuts + inputCuts + "&&";

  cuts+="w>0&&layer==2";
  
  TCanvas *c = new TCanvas("cSeeds","cuts seeds",1200,1300);
 
  c->Divide(3,2);  
  
  int ipad=1;
  
  
  c->cd(ipad++);gPad->SetLogy();
  nt->Draw("dPL",cuts+"&&fabs(dPL)<0.1");
  
  c->cd(ipad++);gPad->SetLogy();
  nt->Draw("dTL",cuts+"&&fabs(dTL)<50.");
 
  c->cd(ipad++);gPad->SetLogy();
  nt->Draw("duz",cuts+"&&fabs(duz)<30.");

  c->cd(ipad++);gPad->SetLogy();
  nt->Draw("pt",cuts+"&&pt<5.1");

  c->cd(ipad++);gPad->SetLogy();
  nt->Draw("fitpt",cuts+"&&fitpt<5.");

  return 0;
}
