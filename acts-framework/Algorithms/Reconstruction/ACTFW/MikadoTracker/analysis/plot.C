TCanvas *c = 0;
TNtuple *eff = 0;

void drawXY( const char *cut="part==-2")
{
  if( !c || !eff ) return;
  c->cd(1);
  eff->SetMarkerSize(0.2);
  eff->SetMarkerColor(kBlack);
  eff->Draw("y:x","part%20==0","");
  eff->SetMarkerSize(0.7);
  eff->SetMarkerColor(kRed);
  eff->Draw("y:x",cut,"same");
  c->Update();
}

void drawRZ( const char *cut = "part==-2")
{
  if( !c || !eff ) return;
  c->cd(2);
  eff->SetMarkerSize(0.2);
  eff->SetMarkerColor(kBlack);
  eff->Draw("r:z","part%1==0","");
  eff->SetMarkerSize(0.7);
  eff->SetMarkerColor(kRed);
  eff->Draw("r:z",cut,"same");
  c->Update();
}

void draw( const char *cut = "part==-2" )
{
  drawXY( cut );
  drawRZ( cut );
}

void drawXY(int part)
{
  TString cut = "part==";
  cut+=part;
  drawXY(cut.Data());
}

void drawRZ(int part)
{
  TString cut = "part==";
  cut+=part;
  drawRZ(cut.Data());
}

void draw(int part)
{
  drawXY( part );
  drawRZ( part );
}

void plot(){
  TFile *f = new TFile("plot.root","READ");  
  if( !f ) return;
  c = new TCanvas("c","c",0,0,2000,1000);
  c->Divide(2,1);   
  eff = (TNtuple*) f->FindObjectAny("eff");
  draw();
}
