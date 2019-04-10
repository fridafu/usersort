{
   gROOT->Reset();
   gROOT->SetStyle("Plain");
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);
   gStyle->SetFillColor(0);
   gStyle->SetPadBorderMode(0);
   m = (TH1F*)gROOT->FindObject("h");
   if (m) m->Delete();
   TCanvas *c1 = new TCanvas("c1","Normalization of gamma-transmission coefficient",600,600);
   TH2F *h = new TH2F("h"," ",10,-0.752000,  11.048,50,1.615e+00,1.731e+05);
   ifstream sigfile("sigpaw.cnt");
   float sig[92],sigerr[92];
   float energy[167],energyerr[167];
   float extL[168],extH[168];
   int i;
   float a0 = -0.7520;
   float a1 =  0.1200;
   for(i = 0; i < 96; i++){
   	energy[i] = a0 + (a1*i);
   	energyerr[i] = 0.0;
   	extL[i] = 0.0;
   	extH[i] = 0.0;
   }
   float x, y;
   i = 0;
   while(sigfile){
   	sigfile >> x;
   	if(i<91){
   		sig[i]=x;
   	}
   	else{sigerr[i-91]=x;}
   	i++;
   }
   ifstream extendfile("extendLH.cnt");
   i = 0;
   while(extendfile){
   	extendfile >> x >> y ;
   	extL[i]=x;
   	extH[i]=y;
   	i++;
   }
   TGraph *extLgraph = new TGraph(96,energy,extL);
   TGraph *extHgraph = new TGraph(96,energy,extH);
   TGraphErrors *sigexp = new TGraphErrors(91,energy,sig,energyerr,sigerr);
   c1->SetLogy();
   c1->SetLeftMargin(0.14);
   h->GetXaxis()->CenterTitle();
   h->GetXaxis()->SetTitle("#gamma-ray energy E_{#gamma} (MeV)");
   h->GetYaxis()->CenterTitle();
   h->GetYaxis()->SetTitleOffset(1.4);
   h->GetYaxis()->SetTitle("Transmission coeff. (arb. units)");
   h->Draw();
   sigexp->SetMarkerStyle(21);
   sigexp->SetMarkerSize(0.8);
   sigexp->Draw("P");
   extLgraph->SetLineStyle(1);
   extLgraph->DrawGraph(31,&extLgraph->GetX()[0],&extLgraph->GetY()[0],"L");
   extHgraph->SetLineStyle(1);
   extHgraph->DrawGraph(18,&extHgraph->GetX()[78],&extHgraph->GetY()[78],"L");
   TArrow *arrow1 = new TArrow(2.008e+00,1.802e+02,2.008e+00,5.393e+01,0.02,">");
   arrow1->Draw();
   TArrow *arrow2 = new TArrow(2.848e+00,4.084e+02,2.848e+00,1.223e+02,0.02,">");
   arrow2->Draw();
   TArrow *arrow3 = new TArrow(8.608e+00,1.571e+04,8.608e+00,4.703e+03,0.02,">");
   arrow3->Draw();
   TArrow *arrow4 = new TArrow(9.808e+00,1.538e+04,9.808e+00,4.603e+03,0.02,">");
   arrow4->Draw();
   c1->Update();
   c1->Print("sigext.pdf");
   c1->Print("sigext.eps");
   c1->Print("sigext.ps");
}
