{
   gROOT->Reset();
   gROOT->SetStyle("Plain");
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);
   gStyle->SetFillColor(0);
   gStyle->SetPadBorderMode(0);
   m = (TH1F*)gROOT->FindObject("h");
   if (m) m->Delete();
   TCanvas *c1 = new TCanvas("c1","Normalization of level density",600,600);
   TH2F *h = new TH2F("h"," ",10,-0.836000,10.684000,50,0.034348,214200.000000);
   ifstream rholev("rholev.cnt"), rhopaw("rhopaw.cnt"), fermi("fermigas.cnt");
   float levels[96],rho[96],rhoerr[96],energy[783],energyerr[783],fermigas[783];
   float Bn[1]={10.198000};
   float Bnerr[1]={0.001};
   float rho_Bn[1]={21420.000000};
   float rho_Bnerr[1]={3381.000000};
   int i = 0;
   float a0 =  -0.8360;
   float a1 =   0.1280;
   float x,y,z;
   while(fermi){
   	fermi >> x;
   	fermigas[i]=x;
   	energy[i]=a0+(a1*i);
   	energyerr[i]=0.0;
      i++;
   }
   i=0;
   while(rhopaw){
   	rhopaw >> y;
   	if(i<95){
   		rho[i]=y;
   	}
   	else{rhoerr[i-95]=y;}
   	i++;
   }
  	i=0;
	while(rholev){
		rholev >> z;
		levels[i]=z;
		i++;
  }
   TGraphErrors *rhoexp = new TGraphErrors(95,energy,rho,energyerr,rhoerr);
   TGraphErrors *rhoBn = new TGraphErrors(1,Bn,rho_Bn,Bnerr,rho_Bnerr);
   TGraph *fermicalc = new TGraph(782,energy,fermigas);
   TGraph *level = new TGraph(95,energy,levels);
   c1->SetLogy();
   c1->SetLeftMargin(0.14);
   h->GetXaxis()->CenterTitle();
   h->GetXaxis()->SetTitle("Excitation energy E (MeV)");
   h->GetYaxis()->CenterTitle();
   h->GetYaxis()->SetTitleOffset(1.4);
   h->GetYaxis()->SetTitle("Level density #rho (E) (MeV^{-1})");
   h->Draw();
   rhoexp->SetMarkerStyle(21);   rhoexp->SetMarkerSize(0.8);
   rhoexp->Draw("P");
   fermicalc->SetLineStyle(2);
   fermicalc->DrawGraph(39,&fermicalc->GetX()[52],&fermicalc->GetY()[52],"L");
   level->SetLineStyle(1);
   level->Draw("L");
   rhoBn->SetMarkerStyle(25);
   rhoBn->SetMarkerSize(0.8);
   rhoBn->Draw("P");
   TLegend *leg = new TLegend(0.15,0.70,0.6,0.85);
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->AddEntry(rhoexp," Oslo data ","P");
   leg->AddEntry(level," Known levels ","L");
   leg->AddEntry(fermicalc," CT or FG model ","L");	
   leg->AddEntry(rhoBn," #rho from neutron res. data ","P");
   leg->Draw();
   TLatex t;
   t.SetTextSize(0.05);
   t.DrawLatex(    8.547,1.071e+05,"^{xx}Yy");
   TArrow *arrow1 = new TArrow(1.084000,57.565828,1.084000,16.868746,0.02,">");
   arrow1->Draw();
   TArrow *arrow2 = new TArrow(3.644000,728.394760,3.644000,213.444443,0.02,">");
   arrow2->Draw();
   TArrow *arrow3 = new TArrow(7.100000,16192.492640,7.100000,4744.951172,0.02,">");
   arrow3->Draw();
   TArrow *arrow4 = new TArrow(8.380000,55213.506215,8.380000,16179.435547,0.02,">");
   arrow4->Draw();
   c1->Update();
   c1->Print("counting.pdf");
   c1->Print("counting.eps");
   c1->Print("counting.ps");
}
